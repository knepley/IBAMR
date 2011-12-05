// Filename: ConstraintIBMethod.C
// Written by amneet bhalla on meghnad@mech.northwestern.edu
// Created on 12/01/2011.

// This is a concrete class which implements fast and efficient distributed Lagrange multiplier (fictitious domain) method
// for immersed structures.
     
// References 
//  *  Patankar et al. A new formulation of the distributed Lagrange multiplier/fictitious domain method 
//     for particulate flows. Int. Journal of Multiphase flows, 26, 1509-1524 (2000).
//  *  Shirgaonkar et al. A new mathematical formulation and fast algorithm for fully resolved simulation of self-propulsion.
//     JCP, 228 , 2366-2390 (2009).


#include "ConstraintIBMethod.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMARAI INCLUDES
#include <HierarchyDataOpsManager.h>
#include <VariableDatabase.h>
#include <tbox/Utilities.h>
#include <tbox/SAMRAI_MPI.h>
#include <CartesianGridGeometry.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenPatchStrategy.h>
#include <CoarsenSchedule.h>
#include <RefineAlgorithm.h>
#include <RefinePatchStrategy.h>
#include <RefineSchedule.h>
#include <CoarsenOperator.h>
#include <RefineOperator.h>

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES


// C++ INCLUDES
#include <limits>
#include <sstream>


/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
  
  // Number of ghost cells used for each variable quantity.
  static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1 );
  static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);


  // Type of coarsening to perform prior to setting coarse-fine boundary and
  // physical boundary ghost cell values.
  static const std::string CELL_DATA_COARSEN_TYPE = "CUBIC_COARSEN";
  static const std::string SIDE_DATA_COARSEN_TYPE = "CUBIC_COARSEN";

  // Type of extrapolation to use at physical boundaries.
  static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

  // Whether to enforce consistent interpolated values at Type 2 coarse-fine
  // interface ghost cells.
  static const bool CONSISTENT_TYPE_2_BDRY = false;
  
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

ConstraintIBMethod::ConstraintIBMethod(
    const std::string& object_name,  
    Pointer< Database> input_db,
    bool register_for_restart,
    Pointer< INSHierarchyIntegrator > ins_hier_integrator,
    Pointer< IBHierarchyIntegrator  > ib_hier_integrator) 
    : IBMethod(object_name, input_db, register_for_restart),
    d_ins_hier_integrator(ins_hier_integrator),
    d_ib_hier_integrator(ib_hier_integrator),
    d_hier_math_ops(new HierarchyMathOps(object_name+ "HierarchyMathOps",d_hierarchy)),
    d_no_structures(d_l_data_manager->getLagrangianStructureNames(d_hierarchy->getFinestLevelNumber()).size()),
    d_ib_kinematics_ptr(d_no_structures,Pointer<IBKinematics>(NULL)),
    d_needs_div_free_projection(false),
    d_rigid_trans_vel                       (d_no_structures, std::vector<double>(3,0.0)),
    d_rigid_trans_vel_old                   (d_no_structures, std::vector<double>(3,0.0)),
    d_rigid_rot_vel                         (d_no_structures, std::vector<double>(3,0.0)),
    d_rigid_rot_vel_old                     (d_no_structures, std::vector<double>(3,0.0)),
    d_incremented_angle_from_reference_axis (d_no_structures, std::vector<double>(3,0.0)),
    d_vel_com_def                           (d_no_structures, std::vector<double>(3,0.0)),
    d_vel_com_def_old                       (d_no_structures, std::vector<double>(3,0.0)),
    d_omega_com_def                         (d_no_structures, std::vector<double>(3,0.0)),
    d_omega_com_def_old                     (d_no_structures, std::vector<double>(3,0.0)),
    d_center_of_mass                        (d_no_structures, std::vector<double>(3,0.0)),
    d_tagged_pt_lag_idx                     (d_no_structures, 0),
    d_tagged_pt_position                    (d_no_structures, std::vector<double>(3,0.0)),
    d_rho_fluid(std::numeric_limits<double>::quiet_NaN()), 
    d_mu_fluid (std::numeric_limits<double>::quiet_NaN()),
    d_timestep_counter     (-1),
    d_output_interval      (1),
    d_INS_num_cycles       (2),
    d_INS_current_cycle_num(0),
    d_print_output                  (false),
    d_output_drag_and_kinetic_energy(false),
    d_output_power                  (false),
    d_output_trans_vel              (false),
    d_output_rot_vel                (false),
    d_output_COM_coordinates        (false),
    d_output_MOI                    (false),
    d_dir_name            ("./ConstraintIBMethodDump"),
    d_base_output_filename("ImmersedStructrue"),
    d_l_data_U_interp    (d_hierarchy->getNumberOfLevels(),Pointer<LData>(NULL)),
    d_l_data_U_correction(d_hierarchy->getNumberOfLevels(),Pointer<LData>(NULL)),
    d_l_data_U_new       (d_hierarchy->getNumberOfLevels(),Pointer<LData>(NULL))
{
    // NOTE: Parent class constructor registers class with the restart manager, sets object name. 
    
    
    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);
     
    
    // Obtain the type of INS fluid solver and velocity variable.
    bool collocated_solver = false, staggered_solver = false;
    d_u_fluidSolve_var                             = d_ins_hier_integrator->getVelocityVariable();
    Pointer<CellVariable<NDIM,double> > cc_vel_var = d_u_fluidSolve_var;
    Pointer<SideVariable<NDIM,double> > sc_vel_var = d_u_fluidSolve_var;
    if(!cc_vel_var.isNull()) collocated_solver = true;
    if(!sc_vel_var.isNull()) staggered_solver  = true;
    d_new_context            = d_ins_hier_integrator->getNewContext();
    d_u_fluidSolve_idx       = var_db->mapVariableAndContextToIndex(d_u_fluidSolve_var, d_new_context);

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<CellVariable<NDIM,double> > cc_var      = new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops                              = hier_ops_manager->getOperationsDouble(cc_var, d_hierarchy, true);
    Pointer<SideVariable<NDIM,double> > sc_var      = new SideVariable<NDIM,double>("sc_var");
    d_hier_sc_data_ops                              = hier_ops_manager->getOperationsDouble(sc_var, d_hierarchy, true);
   
    // Initialize  variables & variable contexts associated with projection step.
    VariableDatabase<NDIM>* var_db        = VariableDatabase<NDIM>::getDatabase();
    d_scratch_context                     = var_db->getContext(d_object_name + "::SCRATCH");
    
    if(collocated_solver) d_u_var         = new CellVariable<NDIM,double>(d_object_name + "::u"     );
    if(staggered_solver)  d_u_var         = new SideVariable<NDIM,double>(d_object_name + "::u"     );
    d_Div_u_var                           = new CellVariable<NDIM,double>(d_object_name + "::Div_u" );
    d_phi_var                             = new CellVariable<NDIM,double>(d_object_name + "::phi"   );
    const IntVector<NDIM> cell_ghosts     = CELLG;
    const IntVector<NDIM> side_ghosts     = SIDEG; 
    if(collocated_solver) d_u_scratch_idx = var_db->registerVariableAndContext(d_u_var,    d_scratch_context, cell_ghosts);
    if(staggered_solver) d_u_scratch_idx  = var_db->registerVariableAndContext(d_u_var,    d_scratch_context, side_ghosts);
    d_phi_idx                             = var_db->registerVariableAndContext(d_phi_var,  d_scratch_context, cell_ghosts);
    d_Div_u_scratch_idx                   = var_db->registerVariableAndContext(d_Div_u_var,d_scratch_context, cell_ghosts); 



    // Setup the cell centered Poisson Solver needed for projection.
    if (d_needs_div_free_projection)
    {
        const std::string velcorrection_projection_prefix = "ConstraintIBMethodProjection";
        // Setup the various solver components.
        for (int d = 0; d < NDIM; ++d)
        {
            d_velcorrection_projection_bc_coef.setBoundarySlope(2*d  ,0.0);
            d_velcorrection_projection_bc_coef.setBoundarySlope(2*d+1,0.0);
        }

        d_velcorrection_projection_spec    = new PoissonSpecifications(d_object_name + "::ConstraintIBMethodProjection::Spec");
        d_velcorrection_projection_op      = new CCLaplaceOperator(d_object_name + "ConstraintIBMethodProjection::PoissonOperator",
            *d_velcorrection_projection_spec, &d_velcorrection_projection_bc_coef, true);
        d_velcorrection_projection_op->setHierarchyMathOps(d_hier_math_ops);
        d_velcorrection_projection_solver  = new PETScKrylovLinearSolver(d_object_name + "ConstraintIBMethodProjection::PoissonKrylovSolver",
             velcorrection_projection_prefix);
        d_velcorrection_projection_solver->setInitialGuessNonzero(false);
        d_velcorrection_projection_solver->setOperator(d_velcorrection_projection_op);

       	if (d_velcorrection_projection_fac_pc_db.isNull())
        {
            TBOX_WARNING(d_object_name << "::ConstraintIBMethod():\n" <<
                " ConstraintIBMethodProjection:: Poisson FAC PC solver database is null." << std::endl);
	}

        d_velcorrection_projection_fac_op  = new CCPoissonFACOperator(d_object_name + ":: ConstraintIBMethodProjection::PoissonFACOperator", 
             d_velcorrection_projection_fac_pc_db);
        d_velcorrection_projection_fac_op->setPoissonSpecifications(*d_velcorrection_projection_spec);
        d_velcorrection_projection_fac_pc  = new IBTK::FACPreconditioner(d_object_name + "::ConstraintIBMethodProjection::PoissonPreconditioner",
             *d_velcorrection_projection_fac_op, d_velcorrection_projection_fac_pc_db);
        d_velcorrection_projection_solver->setPreconditioner(d_velcorrection_projection_fac_pc);

        // Set some default options.
        d_velcorrection_projection_solver->setKSPType("gmres");
        d_velcorrection_projection_solver->setAbsoluteTolerance(1.0e-12);
        d_velcorrection_projection_solver->setRelativeTolerance(1.0e-08);
        d_velcorrection_projection_solver->setMaxIterations(25);

        // NOTE: We always use homogeneous Neumann boundary conditions for the
        // velocity correction projection Poisson solver.
        d_velcorrection_projection_solver->setNullspace(true, NULL);
    }
    else
    {
        d_velcorrection_projection_spec   = NULL;
        d_velcorrection_projection_op     = NULL;
        d_velcorrection_projection_fac_op = NULL;
        d_velcorrection_projection_fac_pc = NULL;
        d_velcorrection_projection_solver = NULL;
    }
        
    // Get the variables related to power calculation.
    if(d_output_power)
    {
        d_VisDefPower_var     = new CellVariable<NDIM,double>(d_object_name + "::Power"             , NDIM );
        d_EulDefVel_var       = new CellVariable<NDIM,double>(d_object_name + "::EulerianDefVel"    , NDIM );
        d_ConstraintPower_var = new CellVariable<NDIM,double>(d_object_name + "::ConstraintPower"   , NDIM);
        d_EulTagDefVel_var    = new CellVariable<NDIM,int   >(d_object_name + "::EulerianTagDefVel");
      
        // Map the variables.
        d_VisDefPower_scratch_idx     = var_db->registerVariableAndContext(d_VisDefPower_var,    d_scratch_context, 0);
        d_EulDefVel_scratch_idx       = var_db->registerVariableAndContext(d_EulDefVel_var,      d_scratch_context, 1);
        d_ConstraintPower_scratch_idx = var_db->registerVariableAndContext(d_ConstraintPower_var,d_scratch_context, 0);
        d_EulTagDefVel_scratch_idx    = var_db->registerVariableAndContext(d_EulTagDefVel_var,   d_scratch_context, 1);
    
    }

    // Create several communications algorithms, used in filling ghost cell data
    // and synchronizing data on the patch hierarchy.
    Pointer<Geometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineAlgorithm<NDIM> > refine_alg;
    Pointer<RefineOperator<NDIM> > refine_op;
    RefinePatchStrategy<NDIM>* refine_patch_strategy;
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg;
    Pointer<CoarsenOperator<NDIM> > coarsen_op;

    // Create algorithm for synching:
    // a) u_fluidSolve and interpolating velocity on Lagrangian mesh.
    // b) deformational velocity and tagging of grid cells which have lag pts.
    coarsen_alg = new CoarsenAlgorithm<NDIM>();
    coarsen_op = grid_geom->lookupCoarsenOperator(d_u_fluidSolve_var, "CONSERVATIVE_COARSEN");
    coarsen_alg->registerCoarsen(d_u_fluidSolve_idx, d_u_fluidSolve_idx, coarsen_op);
    d_ib_hier_integrator->registerCoarsenAlgorithm(d_object_name+"Synch::u_fluidSolve", coarsen_alg);
    
    if(d_output_power)
    {
        coarsen_alg = new CoarsenAlgorithm<NDIM>();
        coarsen_op = grid_geom->lookupCoarsenOperator(d_EulDefVel_var, "CONSERVATIVE_COARSEN");
        coarsen_alg->registerCoarsen(d_EulDefVel_scratch_idx,    d_EulDefVel_scratch_idx,    coarsen_op);
	coarsen_alg->registerCoarsen(d_EulTagDefVel_scratch_idx, d_EulTagDefVel_scratch_idx, coarsen_op);
        d_ib_hier_integrator->registerCoarsenAlgorithm(d_object_name+"Synch::u_def/tag", coarsen_alg);
    }
    
    // Create algorithm for filling ghost cells needed during 
    // a) interpolating velocity on Lagrangian mesh.
    // b) calculating viscous deformational power.
    refine_alg = new RefineAlgorithm<NDIM>();
    refine_op = NULL;
    refine_alg->registerRefine(d_u_fluidSolve_idx, d_u_fluidSolve_idx, d_u_fluidSolve_idx, refine_op);
    d_ib_hier_integrator->registerGhostfillRefineAlgorithm(d_object_name+"FILL_GHOSTCELL::u_fluidSolve", refine_alg);

    if(d_output_power)
    {
        refine_alg = new RefineAlgorithm<NDIM>();
        refine_op = NULL;
        refine_alg->registerRefine(d_EulDefVel_scratch_idx,   d_EulDefVel_scratch_idx,    d_EulDefVel_scratch_idx, refine_op);
	refine_alg->registerRefine(d_EulTagDefVel_scratch_idx,d_EulTagDefVel_scratch_idx, d_EulTagDefVel_scratch_idx, refine_op);
        d_ib_hier_integrator->registerGhostfillRefineAlgorithm(d_object_name+"FILL_GHOSTCELL::u_def/tag", refine_alg);
    }
    
    // Create algorithm for spreading correction.
    refine_alg = new RefineAlgorithm<NDIM>();
    refine_op = grid_geom->lookupRefineOperator(d_u_fluidSolve_var, "CONSERVATIVE_LINEAR_REFINE");
    refine_alg->registerRefine(d_u_fluidSolve_idx, d_u_fluidSolve_idx, d_u_fluidSolve_idx, refine_op);
    d_ib_hier_integrator->registerProlongRefineAlgorithm(d_object_name+"Prolong::u_fluidSolve", refine_alg);
    
    // Do printing operation for processor 0 only.
    if( !SAMRAI_MPI::getRank() && d_print_output)
    {
        d_trans_vel_stream         .resize(d_no_structures);
	d_rot_vel_stream           .resize(d_no_structures);
	d_drag_force_stream        .resize(d_no_structures);
	d_moment_of_inertia_stream .resize(d_no_structures);
	d_kinetic_energy_stream    .resize(d_no_structures);
	d_position_COM_stream      .resize(d_no_structures);
	d_power_spent_stream       .resize(d_no_structures);
	
        for(int struct_no = 0; struct_no < d_no_structures; ++struct_no)
	{
	    std::ostringstream trans_vel,rot_vel,drag_force,moment_inertia,
            kinetic_energy, position_com, power_spent;
	  
	    trans_vel      << d_base_output_filename + "_Trans_vel"       + "_struct_no_" << struct_no;
	    rot_vel        << d_base_output_filename + "_Rot_vel"         + "_struct_no_" << struct_no;
	    drag_force     << d_base_output_filename + "_Drag_force"      + "_struct_no_" << struct_no;
	    moment_inertia << d_base_output_filename + "_MOI"             + "_struct_no_" << struct_no;
	    kinetic_energy << d_base_output_filename + "_Kinetic_energy"  + "_struct_no_" << struct_no;
	    position_com   << d_base_output_filename + "_COM_coordinates" + "_struct_no_" << struct_no;
	    power_spent    << d_base_output_filename + "_Power_spent"     + "_struct_no_" << struct_no;
	    
	    if(from_restart) d_trans_vel_stream[struct_no]         = new std::ofstream(trans_vel.str().c_str(), std::fstream::app);
	    else             d_trans_vel_stream[struct_no]         = new std::ofstream(trans_vel.str().c_str(), std::fstream::out);
	    
	    if(from_restart) d_rot_vel_stream[struct_no]           = new std::ofstream(rot_vel.str().c_str(), std::fstream::app);
	    else             d_rot_vel_stream[struct_no]           = new std::ofstream(rot_vel.str().c_str(), std::fstream::out);
	    
	    if(from_restart) d_drag_force_stream[struct_no]        = new std::ofstream(drag_force.str().c_str(), std::fstream::app);
	    else             d_drag_force_stream[struct_no]        = new std::ofstream(drag_force.str().c_str(), std::fstream::out);
	    
	    if(from_restart) d_moment_of_inertia_stream[struct_no] = new std::ofstream(moment_inertia.str().c_str(), std::fstream::app);
	    else             d_moment_of_inertia_stream[struct_no] = new std::ofstream(moment_inertia.str().c_str(), std::fstream::out);
	  
	    if(from_restart) d_kinetic_energy_stream[struct_no]    = new std::ofstream(kinetic_energy.str().c_str(), std::fstream::app);
	    else             d_kinetic_energy_stream[struct_no]    = new std::ofstream(kinetic_energy.str().c_str(), std::fstream::out);
	    
	    if(from_restart) d_position_COM_stream[struct_no]      = new std::ofstream(kinetic_energy.str().c_str(), std::fstream::app);
	    else             d_position_COM_stream[struct_no]      = new std::ofstream(kinetic_energy.str().c_str(), std::fstream::out);
	    
	    if(from_restart) d_power_spent_stream[struct_no]       = new std::ofstream(power_spent.str().c_str(), std::fstream::app);
	    else             d_power_spent_stream[struct_no]       = new std::ofstream(power_spent.str().c_str(), std::fstream::out);
	}

    }
   
   
     // set the initial velocity of lag points.
    setInitialLagrangianVelocity();  
       
    // calculate the volume of material point in the non-elastic domain.
    calculateVolOfLagPointInNonElasticDomain();
    
    // calculate the volume of material point in the elastic domain
    if( d_body_is_partly_elastic )
     { calculateVolOfLagPointInElasticDomain(); }

 return;
  
  
} // ConstraintIBMethod

ConstraintIBMethod::~ConstraintIBMethod()
{
  
    delete d_velcorrection_projection_spec;
    if( !SAMRAI_MPI::getRank() && d_print_output)
    {
        for(int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            delete (d_trans_vel_stream[struct_no]);
            delete (d_rot_vel_stream[struct_no]);
            delete (d_drag_force_stream[struct_no]);
            delete (d_moment_of_inertia_stream[struct_no]);
            delete (d_kinetic_energy_stream[struct_no]);
            delete (d_position_COM_stream[struct_no]);
            delete (d_power_spent_stream[struct_no]); 
        }
    }
    
    return;
} //~ConstraintIBMethod

/////////////////////////////// PRIVATE ///////////////////////////////////////
void 
ConstraintIBMethod::getFromInput(
    Pointer< Database > input_db, 
    const bool from_restart)
{

      //Read in control parameters from input database.
    d_INS_num_cycles                   = input_db->getIntegerWithDefault("num_INS_cycles", d_INS_num_cycles);
    d_needs_div_free_projection        = input_db->getBoolWithDefault("needs_divfree_projection", d_needs_div_free_projection);
    d_rho_fluid                        = input_db->getDoubleWithDefault("rho_fluid", d_rho_fluid);
    d_mu_fluid                         = input_db->getDoubleWithDefault("mu_fluid" , d_mu_fluid);
    d_lag_position_update_method       = input_db->getString("lag_position_update_method");
   
    // check if lagrangian update method is consistent.
    if(d_lag_position_update_method != "user_defined_velocity" && d_lag_position_update_method != "user_defined_position" 
       && d_lag_position_update_method != "background_fluid_velocity" )
    {
      TBOX_ERROR( "ERROR::ConstraintIBMethod::getFromInput( ) " << "\n"
                  << "Update methods supported are user_defined_velocity/user_defined_position/background_fluid_velocity \n\n"
		  << std::endl);
    }
    //Printing stuff to files.
    Pointer<Database> d_output_db     = input_db->getDatabase("PrintOutput");
    d_print_output                    = d_output_db->getBoolWithDefault( "print_output"  ,d_print_output);
    d_output_interval                 = d_output_db->getIntegerWithDefault("output_interval", d_output_interval);
    d_output_drag_and_kinetic_energy  = d_output_db->getBoolWithDefault( "output_drag_kinetic_energy", d_output_drag_and_kinetic_energy);
    d_output_power                    = d_output_db->getBoolWithDefault( "output_power" , d_output_power);
    d_output_trans_vel                = d_output_db->getBoolWithDefault( "output_rig_transvel" , d_output_trans_vel);
    d_output_rot_vel                  = d_output_db->getBoolWithDefault( "output_rig_rotvel" , d_output_rot_vel);
    d_output_COM_coordinates          = d_output_db->getBoolWithDefault( "output_com_coords" , d_output_COM_coordinates);
    d_output_MOI                      = d_output_db->getBoolWithDefault( "output_moment_inertia" , d_output_MOI);
    d_dir_name                        = d_output_db->getStringWithDefault("output_dirname", d_dir_name) + "/";
    d_base_output_filename            = d_dir_name + d_output_db->getStringWithDefault("base_filename",d_base_output_filename);
    
    if(!from_restart) tbox::Utilities::recursiveMkdir(d_dir_name);
  
    return;
} //getFromInput





} //IBAMR

