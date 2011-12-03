
/******************************************************************************************************************************
   Filename: RigidityConstraint.C
   Written by amneet bhalla on meghnad@mech.northwestern.edu
   Created on 01/10/2011.


   This is a concrete class which implements fast and efficient distributed Lagrange multiplier (fictitious domain) method.
   
   References 
   1. Patankar et al. A new formulation of the distributed Lagrange multiplier/fictitious domain method 
      for particulate flows. Int. Journal of Multiphase flows, 26, 1509-1524 (2000).
   2. Shirgaonkar et al. A new mathematical formulation and fast algorithm for fully resolved simulation of self-propulsion.
      JCP, 228 , 2366-2390 (2009).
 
*******************************************************************************************************************************/


#include "RigidityConstraint.h"


//////////////////////////// INCLUDES /////////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LNodeIndexData.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CoarsenOperator.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <HierarchyDataOpsManager.h>
#include <Index.h>
#include <IndexData.h>
#include <Patch.h>
#include <VariableDatabase.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/PIO.h>

// BLITZ++ INCLUDES
#include <blitz/array.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <iterator>
#include <limits>
#include <numeric>
#include <cmath>

using namespace blitz;


// FORTRAN ROUTINES 
#if (NDIM == 2 )

#define CALCULATE_VISC_DEF_POWER_FC FC_FUNC(calculateviscdefpower2d, CALCULATEVISCDEFPOWER2D )

#endif //if (NDIM == 2 )

#if (NDIM == 3 )

#define CALCULATE_VISC_DEF_POWER_FC FC_FUNC(calculateviscdefpower3d, CALCULATEVISCDEFPOWER3D )

#endif //if (NDIM == 3 )


extern "C"
{
    
  void
  CALCULATE_VISC_DEF_POWER_FC(
                     const int* Tag,                           // cells containing lagrangian pts.
		     const double* Vel ,                       // Lagrangian velocity in the cells
		     double*       VisDef_Power,               // Power to be calculated.
		     const double& mu,                         // viscosity of the fluid
		     const int& GCW,                           // Ghost cell width used
		     const int& ilower0, const int& iupper0,   // x-index of patch box 
                     const int& ilower1, const int& iupper1,   // y-index of patch box 
#if (NDIM == 3 )
		     const int& ilower2, const int& iupper2,   // z-index of the patch box
#endif
                     const double* dx                          // grid-spacing at the finest level. 
                   );
  

  
}

///////////////////////////////////////// NAMESPACE ////////////////////////////

namespace IBAMR
{


///////////////////////// Call Back Function to be hooked into IBAMR ///////////

 void callRigidityConstraintCallBackFunction(const double current_time,
                                             const double new_time,
                                             const int cycle_num,
                                             void* ctx )
 {

    RigidityConstraint* rigConstraint_ptr = (RigidityConstraint*)ctx;
	 
    // get control parameters from the RigidityConstraint object.
    static const bool body_is_self_translating          = rigConstraint_ptr->bodyIsSelfTranslating();
    static const bool body_is_self_rotating             = rigConstraint_ptr->bodyIsSelfRotating();
    static const bool needs_divfree_projection          = rigConstraint_ptr->needs_DivFreeProjection(); 
    static const int  INS_num_cycles                    = rigConstraint_ptr->getNoOfINSCycles();
    static const std::string lag_position_update_method = rigConstraint_ptr->getNameOfLagrangianPositionUpdateMethod();
	 
	 
    rigConstraint_ptr->setSimulationTime(new_time);
    rigConstraint_ptr->setHierarchyRelatedData(cycle_num);
	 
    rigConstraint_ptr->interpolateFluidSolveVelocity();
    rigConstraint_ptr->calculateCOMandMOIOfNonElasticBody();
    rigConstraint_ptr->calculateKinematicsVelocity(current_time, new_time);
	 
    if( body_is_self_translating)
     {
       rigConstraint_ptr->subtractDeformationVelocity();
       rigConstraint_ptr->solveRigidTranslationalVelocity();
       
       if(body_is_self_rotating)
	{
          rigConstraint_ptr->solveRigidRotationalVelocity();
	}
	    
     }  

    rigConstraint_ptr->correctVelocityOnLagrangianMesh(); 

         /*if ( body_is_freely_swimming )
             rigConstraint_ptr->ensureMomentumConservation();*/

    rigConstraint_ptr->spreadCorrectedLagrangianVelocity();
    rigConstraint_ptr->synchronizeLevels();	

    if(needs_divfree_projection )
    { rigConstraint_ptr->applyProjection(); }	
	
    if (lag_position_update_method == "user_defined_position")	
    { rigConstraint_ptr->updateLagrangianMarkerPosition(current_time, new_time); }
	 
    if(cycle_num == (INS_num_cycles -1) )
     { rigConstraint_ptr->postprocessData(current_time, new_time); }
       

   return;

}// callRigidityConstraintCallBackFunction

//////////////////////////////////////// STATIC ////////////////////////////////

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


  // PI
  static const double PI = 3.141592653589793238462643383279502884197169399375;


  // Routine to solve 3X3 equation to get rigid body rotational velocity.
  std::vector<double>
  solveSystemOfEqns(std::map<std::string,double>& inertiaTensor, 
                  const double d1,
                  const double d2,
                  const double d3)
  {

      std::vector<double> solution(3,std::numeric_limits<double>::signaling_NaN());

      const double a1 = inertiaTensor["Ixx"],
                   a2 = inertiaTensor["Ixy"],
                   a3 = inertiaTensor["Ixz"],
                   b1 = inertiaTensor["Iyx"],
                   b2 = inertiaTensor["Iyy"],
                   b3 = inertiaTensor["Iyz"],
                   c1 = inertiaTensor["Izx"],
                   c2 = inertiaTensor["Izy"],
                   c3 = inertiaTensor["Izz"];

      const double Dnr       =  (a3*b2*c1 -a2*b3*c1  -a3*b1*c2 + a1*b3*c2 + a2*b1*c3 -a1*b2*c3);
  
      solution[0]            =  (b3*c2*d1 - b2*c3*d1 -a3*c2*d2 + a2*c3*d2 + a3*b2*d3 -a2*b3*d3)/Dnr;

      solution[1]            =  -(b3*c1*d1 - b1*c3*d1 -a3*c1*d2 + a1*c3*d2 + a3*b1*d3 -a1*b3*d3)/Dnr;
 
      solution[2]            =  (b2*c1*d1 - b1*c2*d1 -a2*c1*d2 + a1*c2*d2 + a2*b1*d3 -a1*b2*d3)/Dnr;

      return solution;

  }



}// namespace anonymous



/////////////////////////////////////////    PUBLIC MEMBER FUNCTIONS //////////////////////////////////

  RigidityConstraint::RigidityConstraint(
     const std::string& object_name,
     Pointer< RigidIBStaggeredHierarchyIntegrator > ib_hier_integrator,
     Pointer< INSStaggeredHierarchyIntegrator > ins_hier_integrator,
     Pointer<IBKinematics> ib_kinematics_ptr,
     Pointer< Database > input_db)
    :d_object_name(object_name),
     d_ib_hier_integrator(ib_hier_integrator),
     d_ins_hier_integrator(ins_hier_integrator),
     d_patch_hierarchy(d_ib_hier_integrator->getPatchHierarchy()),
     d_finest_ln(d_patch_hierarchy->getFinestLevelNumber()),
     d_coarsest_ln(0),
     d_hier_math_ops(d_ins_hier_integrator->getHierarchyMathOps()),
     d_lag_data_manager(d_ib_hier_integrator->getLDataManager()),
     d_kinematics_ptr(ib_kinematics_ptr),
     d_lagrangian_global_nodes(d_lag_data_manager->getNumberOfNodes(d_finest_ln)),
     d_rigid_global_nodes(d_lagrangian_global_nodes),
     d_incremented_angle_from_reference_axis(3),
     d_calculate_translational_momentum(3),
     d_calculate_rotational_momentum(3),
     d_rigid_trans_vel(3,0.0),
     d_rigid_trans_vel_old(3,0.0),
     d_rigid_rot_vel(3,0.0),
     d_rigid_rot_vel_old(3,0.0),
     d_vel_com_def(3,0.0),
     d_vel_com_def_old(3,0.0),
     d_omega_com_def(3,0.0),
     d_omega_com_def_old(3,0.0),
     d_center_of_mass_non_elastic_body(3,0.0),
     d_tagged_pt_idx(0),
     d_tagged_point_position(3,0.0),
     d_rho_fluid(std::numeric_limits<double>::quiet_NaN()),
     d_mu_fluid(std::numeric_limits<double>::quiet_NaN()),
     d_timestep_counter(-1),
     d_output_interval(1),
     d_INS_num_cycles(2),
     d_INS_cycle_num(0)
  {
    
   
    // set some default values.
    d_needs_div_free_projection   = true;
    d_body_is_self_rotating       = false;
    d_body_is_self_translating    = false;
    d_kinematics_need_filtering   = false;
    d_body_is_partly_elastic      = false;
    
    d_do_log                          = false;
    d_print_output                    = false;
    d_output_drag_and_kinetic_energy  = false;
    d_output_power                    = false;
    d_output_trans_vel                = false;
    d_output_rot_vel                  = false;
    d_output_COM_coordinates          = false;
    d_output_MOI                      = false;
    d_base_output_filename            = "immersed_body";
    d_dir_name                        = "./OUTPUT" ; 
    // get rest of the values from input file.
    getFromInput(input_db);
    

    // Allocate Memory for lagrangian velocity correction, old velocity, new velocity and deformation velocity data.
    d_lag_vel_correction_data     = d_lag_data_manager->createLNodeLevelData(d_object_name + "::lag_vel_correction_data",
							                     d_finest_ln,NDIM,true);
    d_lag_new_vel_data            = d_lag_data_manager->createLNodeLevelData(d_object_name + "::lag_new_vel_data",
							                     d_finest_ln,NDIM,true);
    d_lag_old_vel_data            = d_lag_data_manager->createLNodeLevelData(d_object_name + "::lag_old_vel_data",
                                                                             d_finest_ln,NDIM,true); 
    d_lag_new_position_data       = d_lag_data_manager->createLNodeLevelData(d_object_name + "::lag_new_position_data",
                                                                             d_finest_ln,NDIM,true);
 
     
    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager   = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<CellVariable<NDIM,double> > cc_var = new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, d_patch_hierarchy, true);
    Pointer<SideVariable<NDIM,double> > sc_var  = new SideVariable<NDIM,double>("sc_var");
    d_hier_sc_data_ops                          = hier_ops_manager->getOperationsDouble(sc_var, d_patch_hierarchy, true);
   
    // Initialize  variables & variable contexts associated with projection step.
    VariableDatabase<NDIM>* var_db        = VariableDatabase<NDIM>::getDatabase();
    d_scratch_context                     = var_db->getContext(d_object_name + "::SCRATCH");
    
    d_U_var                               = new SideVariable<NDIM,double>(d_object_name + "::U"         );
    d_Div_U_var                           = new CellVariable<NDIM,double>(d_object_name + "::Div_U"     );
    d_Phi_var                             = new CellVariable<NDIM,double>(d_object_name + "::Phi"       );
    const IntVector<NDIM> cell_ghosts      = CELLG;
    const IntVector<NDIM> side_ghosts      = SIDEG; 
    d_U_scratch_idx                       = var_db->registerVariableAndContext(d_U_var,    d_scratch_context, side_ghosts);
    d_Phi_idx                             = var_db->registerVariableAndContext(d_Phi_var,  d_scratch_context, cell_ghosts);
    d_Div_U_scratch_idx                   = var_db->registerVariableAndContext(d_Div_U_var,d_scratch_context, cell_ghosts); 

 
    // Get the fluid solve variable and context from INSStaggeredHierarchyIntegrator object 
    // and map it to index.
    d_U_fluidSolve_var       = d_ins_hier_integrator->getVelocityVar();
    d_new_context            = d_ins_hier_integrator->getNewContext();
    d_U_fluidSolve_idx       = var_db->mapVariableAndContextToIndex(d_U_fluidSolve_var, d_new_context);
       
    // Setup the cell centered Poisson Solver needed for projection.
    if (d_needs_div_free_projection)
    {
        const std::string velcorrection_projection_prefix = "velcorrection_projection_";

        // Setup the various solver components.
        for (int d = 0; d < NDIM; ++d)
        {
            d_velcorrection_projection_bc_coef.setBoundarySlope(2*d  ,0.0);
            d_velcorrection_projection_bc_coef.setBoundarySlope(2*d+1,0.0);
        }

        d_velcorrection_projection_spec    = new PoissonSpecifications(d_object_name+"::velcorrection_projection_spec");
        d_velcorrection_projection_op      = new CCLaplaceOperator(d_object_name+"::Velocity Correction Projection Poisson Operator",
                                                                   *d_velcorrection_projection_spec, &d_velcorrection_projection_bc_coef, true);
        d_velcorrection_projection_op->setHierarchyMathOps(d_hier_math_ops);

        d_velcorrection_projection_solver  = new PETScKrylovLinearSolver(d_object_name+"::Velocity Correction Projection Poisson Krylov Solver",
                                                                         velcorrection_projection_prefix);
        d_velcorrection_projection_solver->setInitialGuessNonzero(false);
        d_velcorrection_projection_solver->setOperator(d_velcorrection_projection_op);

       	if (d_velcorrection_projection_fac_pc_db.isNull())
        {
            TBOX_WARNING(d_object_name << "::constructor():\n" <<
                         "  velocity correction projection poisson fac pc solver database is null." << std::endl);

	}

        d_velcorrection_projection_fac_op  = new CCPoissonFACOperator(d_object_name+"::Velocity Correction Projection Poisson FAC Operator", 
                                                                      d_velcorrection_projection_fac_pc_db);
        d_velcorrection_projection_fac_op->setPoissonSpecifications(*d_velcorrection_projection_spec);
        d_velcorrection_projection_fac_pc  = new IBTK::FACPreconditioner(d_object_name+"::Velocity Correction Projection Poisson Preconditioner",
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

    // Get the data transfer operator for synching patch levels after correcting velocity at the finest level.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom     = d_patch_hierarchy->getGridGeometry();
    Pointer<CoarsenOperator<NDIM> > coarsen_operator    = grid_geom->lookupCoarsenOperator(d_U_fluidSolve_var, "CONSERVATIVE_COARSEN");

    // Set the Coarsening operator associated with synching the Eulerian Fluid Velocity.
    d_calgs["SYNCH_EULERIAN_VEL_CORRECTION"]               = new CoarsenAlgorithm<NDIM>();
    d_calgs["SYNCH_EULERIAN_VEL_CORRECTION"]->registerCoarsen(d_U_fluidSolve_idx,  // destination
                                                              d_U_fluidSolve_idx,  // source
			                                      coarsen_operator  ); // operator name
    
        
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
	
     // Refine algorithm and schedule associated with the data copying operation.
      d_PowerRefineAlgorithm.registerRefine( d_EulDefVel_scratch_idx,       //dst
	                                     d_EulDefVel_scratch_idx,       //src
	                                     d_EulDefVel_scratch_idx,       //scratch
					     NULL                           //refine operator
	                                    );
	
      d_PowerRefineAlgorithm.registerRefine( d_EulTagDefVel_scratch_idx,       //dst
	                                     d_EulTagDefVel_scratch_idx,       //src
	                                     d_EulTagDefVel_scratch_idx,       //scratch
					     NULL                              //refine operator
	                                   );
      
    }


   // Do printing operation for processor 0 only.
   if( !SAMRAI_MPI::getRank() && d_print_output)
   {

      std::string  trans_vel      = d_base_output_filename + "_trans_vel";
      std::string  rot_vel        = d_base_output_filename + "_rot_vel"; 
      std::string  drag_force     = d_base_output_filename + "_drag_force";
      std::string  moment_inertia = d_base_output_filename + "_MOI";
      std::string  kinetic_energy = d_base_output_filename + "_kinetic_energy";
      std::string  positioncom    = d_base_output_filename + "_COM_coordinates";
      std::string  powerspent     = d_base_output_filename + "_power_spent";
      
      // Attach streams to the files.
      d_trans_vel_stream.open( trans_vel.c_str(), std::fstream::out );
      d_rot_vel_stream.open( rot_vel.c_str(), std::fstream::out );
      d_drag_force_stream.open( drag_force.c_str(), std::fstream::out); 
      d_moment_of_inertia_stream.open(moment_inertia.c_str(), std::fstream::out);
      d_kinetic_energy_stream.open(kinetic_energy.c_str(), std::fstream::out);
      d_position_COM_stream.open(positioncom.c_str(),  std::fstream::out);
      d_power_spent_stream.open(powerspent.c_str(), std::fstream::out);

   }
   
   
     // set the initial velocity of lag points.
    setInitialLagrangianVelocity();  
       
    // calculate the volume of material point in the non-elastic domain.
    calculateVolOfLagPointInNonElasticDomain();
    
    // calculate the volume of material point in the elastic domain
    if( d_body_is_partly_elastic )
     { calculateVolOfLagPointInElasticDomain(); }

 return;

}// RigidityConstraint


  void
  RigidityConstraint::getFromInput(Pointer< Database > db)
  {

    //Read in control parameters from input database.
    d_INS_num_cycles                   = db->getIntegerWithDefault("num_INS_cycles", d_INS_num_cycles);
    d_needs_div_free_projection        = db->getBoolWithDefault("needs_divfree_projection"       , d_needs_div_free_projection);
    d_kinematics_need_filtering        = db->getBoolWithDefault("filter_deformation_kinematics"  , d_kinematics_need_filtering);
    d_do_log                           = db->getBoolWithDefault("maintain_log"                   , d_do_log);
    d_rho_fluid                        = db->getDoubleWithDefault("rho_fluid", d_rho_fluid);
    d_mu_fluid                         = db->getDoubleWithDefault("mu_fluid" , d_mu_fluid);
    d_calculate_translational_momentum = db->getIntegerArray("calculate_translational_momentum");
    d_calculate_rotational_momentum    = db->getIntegerArray("calculate_rotational_momentum");
    d_lag_position_update_method       = db->getString("lag_position_update_method");
    d_tagged_pt_idx                    = db->getIntegerWithDefault("tag_global_lag_pt", d_tagged_pt_idx);
    
    //CHECK CONSISTENCIES IN INPUT FILE.
    //check if body is self rotating or self translating.
    for(int i = 0; i < 3 ; ++i)
    {
      if(d_calculate_rotational_momentum[i] )
	 d_body_is_self_rotating = true;
      
      if(d_calculate_translational_momentum[i] )
	 d_body_is_self_translating = true;
      
      d_incremented_angle_from_reference_axis[i] = 0.0;
    }
      
    // Only self translating bodies can be self rotating.
    if(d_body_is_self_rotating && !d_body_is_self_translating)
     {
       TBOX_ERROR("ERROR:: RigidityConstraint::getFromInput( ) " << "\n"
              << "Only self translating/propelling bodies can be self rotating." << std::endl );
     }
    // Filtering kinematics velocity makes sense only for cases which are self propelling.
    if(d_kinematics_need_filtering && !d_body_is_self_translating)
     {
       TBOX_ERROR("ERROR:: RigidityConstraint::getFromInput( ) " << "\n"
              << "Deformational kinematics will be filtered only for self translating/propelling bodies." << std::endl );
     }
     
    // check if lagrangian update method is consistent.
    if(d_lag_position_update_method != "user_defined_velocity" && d_lag_position_update_method != "user_defined_position" 
       && d_lag_position_update_method != "background_fluid_velocity" )
    {
      TBOX_ERROR( "ERROR:: RigidityConstraint::getFromInput( ) " << "\n"
                  << "Update methods supported are user_defined_velocity/user_defined_position/background_fluid_velocity \n\n"
		  << std::endl);
    }
    //Printing stuff to files.
    Pointer<Database> d_output_db     = db->getDatabase("PrintOutput");
    d_print_output                    = d_output_db->getBoolWithDefault( "print_output"  ,d_print_output       );
    d_output_interval                 = d_output_db->getIntegerWithDefault("output_interval", d_output_interval);
    d_output_drag_and_kinetic_energy  = d_output_db->getBoolWithDefault( "output_drag_kinetic_energy", d_output_drag_and_kinetic_energy);
    d_output_power                    = d_output_db->getBoolWithDefault( "output_power" , d_output_power);
    d_output_trans_vel                = d_output_db->getBoolWithDefault( "output_rig_transvel" , d_output_trans_vel);
    d_output_rot_vel                  = d_output_db->getBoolWithDefault( "output_rig_rotvel" , d_output_rot_vel);
    d_output_COM_coordinates          = d_output_db->getBoolWithDefault( "output_com_coords" , d_output_COM_coordinates);
    d_output_MOI                      = d_output_db->getBoolWithDefault( "output_moment_inertia" , d_output_MOI);
    d_dir_name                        = d_output_db->getStringWithDefault("output_dirname", d_dir_name) + "/";
    d_base_output_filename            = d_dir_name + d_output_db->getStringWithDefault("base_filename",d_base_output_filename);
    Utilities::recursiveMkdir(d_dir_name);
    
    //check if elastic structure is also there with rigid structure
    d_body_is_partly_elastic          = db->getBoolWithDefault("body_is_partly_elastic", d_body_is_partly_elastic);
    if( d_body_is_partly_elastic)
    {
       Pointer<Database> elasticpart_db    = db->getDatabase("ElasticPartOfBody");
       d_rigid_global_nodes    = elasticpart_db->getIntegerWithDefault( "starting_elastic_lag_idx", d_rigid_global_nodes );
    }
    
  return;

  }//getFromInput


  void
  RigidityConstraint::calculateCOMandMOIOfNonElasticBody()
  {
    
    d_map_lag_position.clear();
    std::vector<std::vector<double> > rig_lag_position_vec;
    std::vector<double> rig_lag_position(NDIM,0.0);   
    std::vector<double> tagged_lag_pt_position(3,0.0);
    
    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
      Pointer<Patch<NDIM> > patch = level->getPatch(p());
      const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
      const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
	
      for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
      {
        const LNodeIndex& node_idx    = *it;
        const int lag_index           = node_idx.getLagrangianIndex();
	const double* const X         = node_idx.getNodeLocation();    
	//store the coordinates of the tagged in the vector.
	if(lag_index == d_tagged_pt_idx)
	{
          for(int dim = 0; dim < NDIM; ++dim)
	   { tagged_lag_pt_position[dim] = X[dim]; }
		
        }    
	if ( lag_index < d_rigid_global_nodes )
	{
	 //get the position of lagrangian nodes             
	  for(int depth = 0; depth < NDIM ; depth++)
	   { rig_lag_position[depth] = X[depth];}
	  
	  rig_lag_position_vec.push_back(rig_lag_position);
          d_map_lag_position.insert(std::make_pair(lag_index,rig_lag_position));
	}

      }// for a patch
      
    }// for all patches.
     
   
   d_rigid_local_nodes = rig_lag_position_vec.size();
   if(!d_body_is_partly_elastic)
   { TBOX_ASSERT(d_rigid_local_nodes == d_lag_data_manager->getNumberOfLocalNodes(d_finest_ln) );}
   
   // store the data in blitz arrays
   d_x_array.resize(d_rigid_local_nodes,NDIM);
   for(int i = 0 ; i < d_rigid_local_nodes; ++i )
    {
      for(int j = 0; j < NDIM; ++j)
       {
          d_x_array(i,j)      = rig_lag_position_vec[i][j];
       }
     }
	 
    calculateCenterOfMassOfNonElasticBody();   
    if(d_body_is_self_rotating)
     { calculateInertiaTensor();}
    
    // store the coordinates of the tagged point on the body.
    if(d_INS_cycle_num == 0)
    {
      SAMRAI_MPI::sumReduction(&tagged_lag_pt_position[0],3);
      d_tagged_point_position = tagged_lag_pt_position;
    }

     return;
    
  }//calculateCOMandMOIOfNonElasticBody
  
    
  void
  RigidityConstraint::calculateMomentumOfKinematicsVelocity()
  {
    
    //calculate velocity of center of mass introduced due to deformation velocities.
    for(int dim = 0; dim < NDIM; ++dim)
    {
      d_vel_com_def[dim]   = 0.0;
      d_omega_com_def[dim] = 0.0;     
    }

    for(std::map<int, std::vector<double> >::const_iterator mitr = (d_kinematics_ptr->d_map_kinematics_vel).begin(); 
	mitr != (d_kinematics_ptr->d_map_kinematics_vel).end(); ++mitr)
    {
       for(int dim = 0; dim < NDIM; ++dim)
       {
	 d_vel_com_def[dim] += (mitr->second)[dim];
       }  
      
    }
    
    for( int dim = 0; dim < NDIM; ++dim)
    {
       d_vel_com_def[dim] /= (d_kinematics_ptr->d_map_kinematics_vel).size(); 
    }
    

    // calculate the rotational momentum 
    if(d_body_is_self_rotating)
    {
            
#if (NDIM == 2)

      double angularMomentum_z = 0.0;
      for(std::map<int, std::vector<double> >::const_iterator mitr = d_map_lag_position.begin(); 
	  mitr != d_map_lag_position.end(); ++mitr)
       {      
	  int lag_idx = mitr->first; 
	  double x_cm = d_map_lag_position[lag_idx][0] - d_center_of_mass_non_elastic_body[0];
	  double y_cm = d_map_lag_position[lag_idx][1] - d_center_of_mass_non_elastic_body[1];
	  angularMomentum_z +=  ( x_cm*(d_kinematics_ptr->d_map_kinematics_vel)[lag_idx][1] 
	                         -y_cm*(d_kinematics_ptr->d_map_kinematics_vel)[lag_idx][0] );
	      
	}// angularMomentumOnThisProcessor.
	
      //calculate omega of body introduced by the deformational kinematics.
      d_omega_com_def[2] = SAMRAI_MPI::sumReduction(angularMomentum_z)/d_map_inertiaTensor["Izz"];

#endif
      
#if (NDIM == 3)
	
     double angularMomentum_x =0.0, angularMomentum_y = 0.0, angularMomentum_z = 0.0;
     for(std::map<int, std::vector<double> >::const_iterator mitr = d_map_lag_position.begin(); 
	    mitr != d_map_lag_position.end(); ++mitr)
      {      
	 int lag_idx = mitr->first; 
	 double x_cm = d_map_lag_position[lag_idx][0] - d_center_of_mass_non_elastic_body[0];
	 double y_cm = d_map_lag_position[lag_idx][1] - d_center_of_mass_non_elastic_body[1];
	 double z_cm = d_map_lag_position[lag_idx][2] - d_center_of_mass_non_elastic_body[2];
	    
	 angularMomentum_x +=  ( y_cm*(d_kinematics_ptr->d_map_kinematics_vel)[lag_idx][2] 
	                        -z_cm*(d_kinematics_ptr->d_map_kinematics_vel)[lag_idx][1] );
	    
	 angularMomentum_y +=  ( -x_cm*(d_kinematics_ptr->d_map_kinematics_vel)[lag_idx][2] 
	                         +z_cm*(d_kinematics_ptr->d_map_kinematics_vel)[lag_idx][0] );
	    
	 angularMomentum_z +=  ( x_cm*(d_kinematics_ptr->d_map_kinematics_vel)[lag_idx][1] 
	                        -y_cm*(d_kinematics_ptr->d_map_kinematics_vel)[lag_idx][0] );
	      
      }// angularMomentumOnThisProcessor.
	
    //add angular momentum on all processors.
     angularMomentum_x = SAMRAI_MPI::sumReduction(angularMomentum_x);
     angularMomentum_y = SAMRAI_MPI::sumReduction(angularMomentum_y);
     angularMomentum_z = SAMRAI_MPI::sumReduction(angularMomentum_z);

   //solve the system of equation i.e Iw = sum_i (r_i X u_def_i )
     d_omega_com_def                  =    solveSystemOfEqns(d_map_inertiaTensor,
                                                              angularMomentum_x, 
                                                              angularMomentum_y,
                                                              angularMomentum_z);
	
#endif
      
      
    }//d_body_is_self_rotating
    
    return;
    
  }//calculateMomentumOfKinematicsVelocity

  void
  RigidityConstraint::filterGivenKinematicsVelocity()
  {
    
    (d_kinematics_ptr->d_map_filtered_kinematics_vel).clear();
    std::vector<double> radiusFromCoM(3,0.0), OmegaCrossR(3,0.0), vec_filtered_def_vel(NDIM);
    // filter the translational momentum and angular momentum due to deformation velocity
    for(std::map<int, std::vector<double> >::const_iterator mitr = d_map_lag_position.begin(); 
	    mitr != d_map_lag_position.end(); ++mitr)
    {
     
       const int lag_idx = mitr->first;
       if(d_body_is_self_rotating)
       {
	  for(int depth = 0; depth < NDIM; ++depth)
	   {
             radiusFromCoM[depth] = d_map_lag_position[lag_idx][depth] - d_center_of_mass_non_elastic_body[depth];
	   }
	    
	  OmegaCrossR[0] =  radiusFromCoM[2]*d_omega_com_def[1] - radiusFromCoM[1]*d_omega_com_def[2];
          OmegaCrossR[1] = -radiusFromCoM[2]*d_omega_com_def[0] + radiusFromCoM[0]*d_omega_com_def[2];
          OmegaCrossR[2] =  radiusFromCoM[1]*d_omega_com_def[0] - radiusFromCoM[0]*d_omega_com_def[1];
	    
	  for(int dim = 0; dim < NDIM; ++dim)
	   {
	      vec_filtered_def_vel[dim] = (d_kinematics_ptr->d_map_kinematics_vel)[lag_idx][dim] 
	                                  -( OmegaCrossR[dim] + d_vel_com_def[dim] );
	   }
	    
	  (d_kinematics_ptr->d_map_filtered_kinematics_vel).insert(std::make_pair(lag_idx,vec_filtered_def_vel));
       }
       else
       {
	 
	  for(int dim = 0; dim < NDIM; ++dim)
	   {
	       vec_filtered_def_vel[dim] = (d_kinematics_ptr->d_map_kinematics_vel)[lag_idx][dim] - d_vel_com_def[dim] ;
	   }
	   
	  (d_kinematics_ptr->d_map_filtered_kinematics_vel).insert(std::make_pair(lag_idx,vec_filtered_def_vel)); 
       }
       
       
    }
    
  
#ifdef DEBUG_CHECK_ASSERTIONS

   /*!
    * check to see if this transformation has resulted in no translational and 
    * rotational momentum in the kinematics velocity.
    */

   std::vector<double> vel_com_def(3,0.0), omega_com_def(3,0.0);	
   //calculate velocity of center of mass introduced due to deformation velocities.

   for(std::map<int, std::vector<double> >::const_iterator mitr = d_map_lag_position.begin(); 
	    mitr != d_map_lag_position.end(); ++mitr) 
    {
      
       const int lag_idx = mitr->first;
       for(int dim = 0; dim < NDIM; ++dim)
       {
	 vel_com_def[dim] += (d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_idx][dim];
       }  
      
    }
    
    SAMRAI_MPI::sumReduction(&vel_com_def[0],NDIM);
    
    // calculate the rotational momentum 
    if(d_body_is_self_rotating)
    {
            
#if (NDIM == 2)

       double angularMomentum_z = 0.0;
       for(std::map<int, std::vector<double> >::const_iterator mitr = d_map_lag_position.begin(); 
	    mitr != d_map_lag_position.end(); ++mitr)
        {
	      
	  int lag_idx = mitr->first; 
	  double x_cm = d_map_lag_position[lag_idx][0] - d_center_of_mass_non_elastic_body[0];
	  double y_cm = d_map_lag_position[lag_idx][1] - d_center_of_mass_non_elastic_body[1];
	  angularMomentum_z +=  ( x_cm*(d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_idx][1] 
	                         -y_cm*(d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_idx][0] );
	      
	}// angularMomentumOnThisProcessor.
	
	//calculate omega of body introduced by the deformational kinematics.
	omega_com_def[2] = SAMRAI_MPI::sumReduction(angularMomentum_z)/d_map_inertiaTensor["Izz"];

#endif
      
#if (NDIM == 3)
	
     double angularMomentum_x =0.0, angularMomentum_y = 0.0, angularMomentum_z = 0.0;
     for(std::map<int, std::vector<double> >::const_iterator mitr = d_map_lag_position.begin(); 
	    mitr != d_map_lag_position.end(); ++mitr)
      {      
	 int lag_idx = mitr->first; 
	 double x_cm = d_map_lag_position[lag_idx][0] - d_center_of_mass_non_elastic_body[0];
	 double y_cm = d_map_lag_position[lag_idx][1] - d_center_of_mass_non_elastic_body[1];
	 double z_cm = d_map_lag_position[lag_idx][2] - d_center_of_mass_non_elastic_body[2];
	    
	 angularMomentum_x +=  ( y_cm*(d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_idx][2] 
	                        -z_cm*(d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_idx][1] );
	    
	 angularMomentum_y +=  ( -x_cm*(d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_idx][2] 
	                         +z_cm*(d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_idx][0] );
	    
	 angularMomentum_z +=  ( x_cm*(d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_idx][1] 
	                        -y_cm*(d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_idx][0] );
	      
     }// angularMomentumOnThisProcessor.
	
     //add angular momentum on all processors.
     angularMomentum_x = SAMRAI_MPI::sumReduction(angularMomentum_x);
     angularMomentum_y = SAMRAI_MPI::sumReduction(angularMomentum_y);
     angularMomentum_z = SAMRAI_MPI::sumReduction(angularMomentum_z);

    //solve the system of equation i.e Iw = sum_i (r_i X u_def_i )
     omega_com_def     =    solveSystemOfEqns(d_map_inertiaTensor,
                                              angularMomentum_x, 
                                              angularMomentum_y,
                                              angularMomentum_z);
	
#endif
           
    }//bodyIsRotating
    
    for( int depth = 0 ; depth < 3; ++depth)
    {
      TBOX_ASSERT( MathUtilities<double>::equalEps(vel_com_def[depth], 0.0)   );
      TBOX_ASSERT( MathUtilities<double>::equalEps(omega_com_def[depth], 0.0) );
    }
    
#endif

    return;
    
  }//filterGivenKinematicsVelocity



  void
  RigidityConstraint::setInitialLagrangianVelocity()
  {
      
    calculateCOMandMOIOfNonElasticBody();     
    d_kinematics_ptr->calculateGivenKinematicsVelocity(0.0,d_incremented_angle_from_reference_axis,
						       d_center_of_mass_non_elastic_body, d_tagged_point_position);
    if( d_body_is_self_translating )
     { calculateMomentumOfKinematicsVelocity(); }
    if(d_kinematics_need_filtering)
     { filterGivenKinematicsVelocity();}
       
   // Get the background fluid velocity at the lagrangian points.
   getBackgroundFluidVelocity(d_lag_old_vel_data, 0.0);
   
    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
     {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
         {
            const LNodeIndex& node_idx    = *it;
            const int lag_index           = node_idx.getLagrangianIndex();
            const int local_petsc_index   = node_idx.getLocalPETScIndex();

            // fill in the kinematic velocities at the LNodeLevelData
            if( (lag_index < d_rigid_global_nodes))
            {
	       if(d_kinematics_need_filtering)
	       {
		  for(int depth = 0; depth < NDIM ; depth++)
		   {
		     (*d_lag_old_vel_data)(local_petsc_index,depth) = (d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_index][depth];
		   }	  
	       }
	      else
	      {
	          for(int depth = 0; depth < NDIM ; depth++)
		   {
		      (*d_lag_old_vel_data)(local_petsc_index,depth) = (d_kinematics_ptr->d_map_kinematics_vel)[lag_index][depth];
		   }
	      }
		
            }// for non elastic domain.

	 }// over a patch

     }// over patches


    return;
  }//setInitialLagrangianVelocity

  
  void
  RigidityConstraint::calculateVolOfLagPointInNonElasticDomain()
  {

    double volume_body_onthisproc = 0.0;
    Pointer< PatchLevel<NDIM> > finest_level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
    const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
	
     // Initialize variables and variable contexts associated with Eulerian tracking of the lagrangian points.
    VariableDatabase<NDIM>* var_db        = VariableDatabase<NDIM>::getDatabase();
    const IntVector<NDIM> cell_ghosts   = 0; 
    Pointer< CellVariable< NDIM,int > > LagPtsInEulerianCell_var = new CellVariable<NDIM, int>(d_object_name + "::LagPtsInEulerianCellNonElasticDomain");
    const int LagPtsInEulerianCell_scratch_idx                   = var_db->registerVariableAndContext(LagPtsInEulerianCell_var, 
												       d_scratch_context, cell_ghosts);	
    finest_level->allocatePatchData(LagPtsInEulerianCell_scratch_idx, 0.0);
 
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
       Pointer<CellData<NDIM,int   > > LagPtsInEulerianCell_scratch_idx_Data  = patch->getPatchData(LagPtsInEulerianCell_scratch_idx);
       LagPtsInEulerianCell_scratch_idx_Data->fill(0,0);
     }


    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
       Pointer<CellData<NDIM,int   > > LagPtsInEulerianCell_scratch_idx_Data  = patch->getPatchData(LagPtsInEulerianCell_scratch_idx);
       Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
       const double* const dx                       = pgeom->getDx(); 
       const double* const XLower                   = pgeom->getXLower();
       const Box<NDIM>& patch_box                   = patch->getBox();
       const Index<NDIM>& ilower                    = patch_box.lower();
	    
       d_finest_cell_vol                            = dx[0]*dx[1]
#if (NDIM > 2)
                                                      *dx[2]
#endif
                                                            ;
							   
       const Pointer<LNodeIndexData> idx_data       = patch->getPatchData(lag_node_index_idx);	
       for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
                 it != idx_data->lnode_index_end(); ++it)
        {
	
	   const LNodeIndex& node_idx    = *it;
           const int lag_index           = node_idx.getLagrangianIndex();
	   const double* const X         = node_idx.getNodeLocation();
		    
	   if(lag_index < d_rigid_global_nodes )
	    {
               const CellIndex<NDIM> Lag2Eul_cellindex(  Index<NDIM> (int( floor((X[0]-XLower[0])/dx[0]))+ilower(0)
#if (NDIM > 1)
                                                         ,int(floor((X[1]-XLower[1])/dx[1]))+ilower(1)
#if (NDIM > 2)
                                                         ,int(floor((X[2]-XLower[2])/dx[2]))+ilower(2)
#endif
#endif
                                                                     )
					               );
		    
		(*LagPtsInEulerianCell_scratch_idx_Data)(Lag2Eul_cellindex)++;
		    
             } // only in the non elastic domain.
		    
	 }// on a patch
	    
	// Iterate over the same patch to find the fraction of volume of the body on this patch.
	 for(CellData<NDIM,double>::Iterator it(patch_box); it; it++)
          {
	    if( (*LagPtsInEulerianCell_scratch_idx_Data)(*it) )
	     { volume_body_onthisproc += d_finest_cell_vol; }
	     
	  }// on the same patch
	    
     }// on all the patches.
	

   // Do the parallel reduction
    d_volume_of_body = SAMRAI_MPI::sumReduction(volume_body_onthisproc);
    d_volume_of_non_elastic_lag_pt = d_volume_of_body/d_rigid_global_nodes;
   
   
   if(d_do_log)
   {
     tbox::plog << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n\n\n" 
                << " VOLUME OF THE CELL AT FINEST LEVEL       = " << d_finest_cell_vol << "\n"
                << " VOLUME OF THE NON ELASTIC MATERIAL POINT = " << d_volume_of_non_elastic_lag_pt << "\n"
                << " VOLUME OF THE NON ELASTIC BODY           = " << d_volume_of_body << "\n"
                << " MASS   OF THE NON ELASTIC BODY           = " << d_rho_fluid*d_volume_of_body << "\n\n\n"
                << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" << std::endl;
   }            
     
   // Deallocate non elastic tag variable data.
   if( finest_level->checkAllocated(LagPtsInEulerianCell_scratch_idx) )
    { 
      finest_level->deallocatePatchData(LagPtsInEulerianCell_scratch_idx); 
      var_db->removePatchDataIndex(LagPtsInEulerianCell_scratch_idx);
    }
    
   return;
   
  } //calculateVolOfLagPointInNonElasticDomain
  

  void
  RigidityConstraint::calculateVolOfLagPointInElasticDomain()
  {

    double volume_body_onthisproc = 0.0;
    Pointer< PatchLevel<NDIM> > finest_level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
    const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
	
   // Initialize variables and variable contexts associated with Eulerian tracking of the lagrangian points.
    VariableDatabase<NDIM>* var_db        = VariableDatabase<NDIM>::getDatabase();
    const IntVector<NDIM> cell_ghosts   = 0; 
    Pointer< CellVariable< NDIM,int > > LagPtsInEulerianCell_var = new CellVariable<NDIM, int>(d_object_name + "::LagPtsInEulerianCellElasticDomain");
    const int LagPtsInEulerianCell_scratch_idx                   = var_db->registerVariableAndContext(LagPtsInEulerianCell_var, 
													  d_scratch_context, cell_ghosts);
    finest_level->allocatePatchData(LagPtsInEulerianCell_scratch_idx, 0.0);
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
       Pointer<CellData<NDIM,int   > > LagPtsInEulerianCell_scratch_idx_Data  = patch->getPatchData(LagPtsInEulerianCell_scratch_idx);
       LagPtsInEulerianCell_scratch_idx_Data->fill(0,0);
     }

    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
       Pointer<CellData<NDIM,int   > > LagPtsInEulerianCell_scratch_idx_Data   = patch->getPatchData(LagPtsInEulerianCell_scratch_idx);
       Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
       const double* const dx                       = pgeom->getDx(); 
       const double* const XLower                   = pgeom->getXLower();
       const Box<NDIM>& patch_box                   = patch->getBox();
       const Index<NDIM>& ilower                    = patch_box.lower();
	    
       d_finest_cell_vol                            = dx[0]*dx[1]
#if (NDIM > 2)
                                                     *dx[2]
#endif
                                                           ;
							   
       const Pointer<LNodeIndexData> idx_data       = patch->getPatchData(lag_node_index_idx);	
       for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
              it != idx_data->lnode_index_end(); ++it)
        {
	  const LNodeIndex& node_idx    = *it;
	  const int lag_index           = node_idx.getLagrangianIndex();
	  const double* const X         = node_idx.getNodeLocation();
		    
          if( lag_index >= d_rigid_global_nodes )
           {
             const CellIndex<NDIM> Lag2Eul_cellindex(  Index<NDIM> (int( floor((X[0]-XLower[0])/dx[0]))+ilower(0)
#if (NDIM > 1)
                                                       ,int(floor((X[1]-XLower[1])/dx[1]))+ilower(1)
#if (NDIM > 2)
                                                       ,int(floor((X[2]-XLower[2])/dx[2]))+ilower(2)
#endif
#endif
                                                                   )
					             );
		    
	     (*LagPtsInEulerianCell_scratch_idx_Data)(Lag2Eul_cellindex)++;
		    
	   } // only in the non elastic domain.
		    
	}// on a patch
	    
       // Iterate over the same patch to find the fraction of volume of the body on this patch.
       for(CellData<NDIM,double>::Iterator it(patch_box); it; it++)
        {
	  if( (*LagPtsInEulerianCell_scratch_idx_Data)(*it) )
	   { volume_body_onthisproc += d_finest_cell_vol;}
	     
	}// on the same patch
	
     }// on all the patches.
	

   // Do the parallel reduction
    const double elastic_volume = SAMRAI_MPI::sumReduction(volume_body_onthisproc);
    d_volume_of_body += elastic_volume;
    d_volume_of_elastic_lag_pt = d_volume_of_body/(d_lagrangian_global_nodes - d_rigid_global_nodes);
    
    if(d_do_log)
     {
       tbox::plog << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n\n\n" 
                  << " VOLUME OF THE ELASTIC MATERIAL POINT = " << d_volume_of_elastic_lag_pt << "\n"
                  << " VOLUME OF THE ELASTIC BODY   = " << elastic_volume << std::endl
                  << " MASS   OF THE ELASTIC BODY   = " << d_rho_fluid*elastic_volume << "\n\n\n"
                  << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
     }            	    
          	
    if( finest_level->checkAllocated(LagPtsInEulerianCell_scratch_idx) )
    { 
      finest_level->deallocatePatchData(LagPtsInEulerianCell_scratch_idx);
      var_db->removePatchDataIndex(LagPtsInEulerianCell_scratch_idx);
    }
    
    return;
  } //calculateVolOfLagPointInElasticDomain



  const std::string
  RigidityConstraint::getNameOfLagrangianPositionUpdateMethod()
  {
     return d_lag_position_update_method;
    
  }// getNameOfLagrangianPositionUpdateMethod
   
  const std::string
  RigidityConstraint::getNameOfLagrangianPositionUpdateDataStructure()
  {
      std::string data_structure_name;
      if (d_lag_position_update_method == "user_defined_velocity" )
      {
	 data_structure_name = d_object_name + "::lag_new_vel_data";
      }
      else if (d_lag_position_update_method == "user_defined_position" )
      {
	 data_structure_name = d_object_name + "::lag_new_position_data";
      }
      else if(d_lag_position_update_method == "background_fluid_velocity" )
      {
	 data_structure_name = "";
      }

    return data_structure_name;
  }// getNameOfLagrangianPositionUpdateDataStructure
 
 
 
  void
  RigidityConstraint::setSimulationTime(const double new_time)
  {
    d_integrator_time = new_time;
    return;
  }//setSimulationTime

  double
  RigidityConstraint::getSimulationTime() const
  {

    return d_integrator_time;
  }//getIntegratorTime


  
 void
 RigidityConstraint::setHierarchyRelatedData(const int cycle_num)
 {
    d_INS_cycle_num = cycle_num;
    if( d_INS_cycle_num ==  (d_INS_num_cycles -1) )
    { ++d_timestep_counter; }
       
       
    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(d_patch_hierarchy);
    d_hier_cc_data_ops->resetLevels(0, d_finest_ln);

    d_hier_sc_data_ops->setPatchHierarchy(d_patch_hierarchy);
    d_hier_sc_data_ops->resetLevels(0, d_finest_ln);

    // Get the control volume weight variables and patch data descriptor
    // indices.
    d_wgt_cc_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Get the volume of the physical domain.
    d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    return;

 } //setHierarchyRelatedData


 void
 RigidityConstraint::calculateKinematicsVelocity(const double current_time, const double new_time)
 {

    const double dt = new_time - current_time;
   // Update the orientation only at the first cycle. Rest of the cycles will have the same orientation 
   // Theta_new = Theta_old + Omega_old*dt
    if (d_INS_cycle_num == 0)
    {
      if(d_kinematics_need_filtering)
       { 
         d_omega_com_def_old               = d_omega_com_def;
	 d_vel_com_def_old                 = d_vel_com_def;
	 for(int i = 0 ; i < d_incremented_angle_from_reference_axis.size(); ++i)
	  { d_incremented_angle_from_reference_axis[i] +=  (d_rigid_rot_vel[i] - d_omega_com_def_old[i])*dt;}
       }
       else
       {
	 for(int i = 0 ; i < d_incremented_angle_from_reference_axis.size(); ++i)
	  { d_incremented_angle_from_reference_axis[i] +=  d_rigid_rot_vel[i]*dt;}
       }
       
       d_kinematics_ptr->calculateGivenKinematicsVelocity(new_time, d_incremented_angle_from_reference_axis, 
							  d_center_of_mass_non_elastic_body, d_tagged_point_position);
    }
    
    if(d_body_is_self_translating)
    { calculateMomentumOfKinematicsVelocity(); }
    
    if(d_kinematics_need_filtering)
    { filterGivenKinematicsVelocity();}


   return;
 }// calculateKinematicsVelocity



  void
  RigidityConstraint::interpolateFluidSolveVelocity()
  {
   

 /*!
   * The following data members are used in spreading/interpolating the deformation velocities 
   * from the lagrangian grid to the Eulerian grid and vice versa.The vector should have all the
   * components for each level of AMR grid.
   */
    std::vector< Pointer< LNodeLevelData > > F_data(d_finest_ln+1,  Pointer<LNodeLevelData>(NULL) );
    std::vector< Pointer< LNodeLevelData > > X_data(d_finest_ln+1,  Pointer<LNodeLevelData>(NULL) );

   // Fill in the above vectors at the finest level.
    F_data[d_finest_ln] = d_lag_vel_correction_data;   
    X_data[d_finest_ln] = d_lag_data_manager->getLNodeLevelData("X",d_finest_ln);

   // Interpolate the fluid solve velocity.
    d_lag_data_manager->interp(d_U_fluidSolve_idx,
			       F_data,
			       X_data,
                               std::vector< Pointer< CoarsenSchedule<NDIM> > > (),
                               std::vector< Pointer< RefineSchedule<NDIM>  > > (),
                               getSimulationTime(),
			       d_finest_ln,
			       d_finest_ln);
    

    // Cache the interpolated data to d_lag_vel_data_cached
     d_lag_vel_correction_data->restoreLocalFormVec();
     (*d_lag_new_vel_data) = (*d_lag_vel_correction_data);

    return;

  }// interpolateFluidSolveVelocity
  
  
  void
  RigidityConstraint::getBackgroundFluidVelocity(Pointer< LNodeLevelData > lag_data, 
						 const double time)
  {
     
     std::vector< Pointer< LNodeLevelData > > F_data(d_finest_ln+1,  Pointer<LNodeLevelData>(NULL) );
     std::vector< Pointer< LNodeLevelData > > X_data(d_finest_ln+1,  Pointer<LNodeLevelData>(NULL) );


   // Fill in the above vectors at the finest level.
    F_data[d_finest_ln] = lag_data;   
    X_data[d_finest_ln] = d_lag_data_manager->getLNodeLevelData("X",d_finest_ln);

   // Interpolate the fluid solve velocity.
    d_lag_data_manager->interp(d_U_fluidSolve_idx,
			       F_data,
			       X_data,
                               std::vector< Pointer< CoarsenSchedule<NDIM> > > (),
                               std::vector< Pointer< RefineSchedule<NDIM>  > > (),
                               time,
			       d_finest_ln,
			       d_finest_ln);
    return;
    
  }//getBackgroundFluidVelocity



 void
 RigidityConstraint::subtractDeformationVelocity()
 {

    std::vector< std::vector<double> > rig_lag_vel_vec;
    std::vector<double> rig_lag_vel(NDIM,0.0);
    
    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = level->getPatch(p());
       const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
       const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
	
       for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
          const LNodeIndex& node_idx    = *it;
          const int lag_index           = node_idx.getLagrangianIndex();
          const int local_petsc_index   = node_idx.getLocalPETScIndex();
	    
	  if ( lag_index < d_rigid_global_nodes )
	   {
	    // fill in the deformation velocities at the LNodeLevelData
	     for(int depth = 0; depth < NDIM ; depth++)
	      {
	        if(d_kinematics_need_filtering)
		 {
		   rig_lag_vel[depth] = (*d_lag_vel_correction_data)(local_petsc_index,depth);
		 }
		else
		 {
		   rig_lag_vel[depth] = (*d_lag_vel_correction_data)(local_petsc_index,depth) 
		                         -(d_kinematics_ptr->d_map_kinematics_vel)[lag_index][depth];
		 }  
		    
	      }
	    rig_lag_vel_vec.push_back(rig_lag_vel);
	   }

	 }// for a patch
     
     }// for all patches.
     
     TBOX_ASSERT(d_rigid_local_nodes == static_cast<int>(rig_lag_vel_vec.size()) );
     
     // store the data in blitz arrays
     d_u_diff_array.resize(d_rigid_local_nodes,NDIM);   
     for(int i = 0; i < d_rigid_local_nodes; ++i)
     {
       for(int j = 0; j < NDIM; ++j)
       {
	 d_u_diff_array(i,j) = rig_lag_vel_vec[i][j];
       }
     }
	 
     return;
 }//subtractDeformationVelocity

 
  void
  RigidityConstraint::solveRigidTranslationalVelocity()
  {

    if(d_INS_cycle_num == 0)
     { d_rigid_trans_vel_old = d_rigid_trans_vel;}
    
#if (NDIM == 2)

    double u_x = 0.0, u_y = 0.0;
    for(int i = 0; i < d_rigid_local_nodes; ++i)
     { 
         u_x += d_u_diff_array(i,0);
         u_y += d_u_diff_array(i,1);
     }
  
    if(d_calculate_translational_momentum[0])
     { d_rigid_trans_vel[0] = SAMRAI_MPI::sumReduction(u_x)/d_rigid_global_nodes; }
    else
     { d_rigid_trans_vel[0] = 0.0; }
    
    if(d_calculate_translational_momentum[1])
     { d_rigid_trans_vel[1] = SAMRAI_MPI::sumReduction(u_y)/d_rigid_global_nodes; }
    else
     { d_rigid_trans_vel[1] = 0.0; }
         
    d_rigid_trans_vel[2] = 0.0;

#endif

#if (NDIM == 3)

    double u_x = 0.0, u_y = 0.0, u_z = 0.0;
    for(int i = 0; i < d_rigid_local_nodes; ++i)
     { 
       u_x += d_u_diff_array(i,0);
       u_y += d_u_diff_array(i,1);
       u_z += d_u_diff_array(i,2);
     }
     
   if(d_calculate_translational_momentum[0])
    { d_rigid_trans_vel[0] = SAMRAI_MPI::sumReduction(u_x)/d_rigid_global_nodes; }
   else
    { d_rigid_trans_vel[0] = 0.0; }
   
   if(d_calculate_translational_momentum[1]) 
    { d_rigid_trans_vel[1] = SAMRAI_MPI::sumReduction(u_y)/d_rigid_global_nodes;}
   else
    { d_rigid_trans_vel[1] = 0.0; }
   
   if(d_calculate_translational_momentum[2]) 
    { d_rigid_trans_vel[2] = SAMRAI_MPI::sumReduction(u_z)/d_rigid_global_nodes;}
   else
    { d_rigid_trans_vel[2] = 0.0; }

#endif


   // write rigid translational velocity to the file
    if( !SAMRAI_MPI::getRank() && d_print_output && d_output_trans_vel && d_INS_cycle_num == (d_INS_num_cycles -1) 
       && (d_timestep_counter % d_output_interval ) == 0)
     {
       d_trans_vel_stream << getSimulationTime()  << '\t' << d_rigid_trans_vel[0] << '\t'
                          << d_rigid_trans_vel[1] << '\t' << d_rigid_trans_vel[2] << '\t'
                          << d_vel_com_def[0]     << '\t' << d_vel_com_def[1]     << '\t'
                          << d_vel_com_def[2]     << std::endl;

     }

 
    return;

  }//solveRigidTranslationalVelocity
  
  void
  RigidityConstraint::calculateCenterOfMassOfNonElasticBody()
  {
      
    for(int depth = 0; depth < NDIM; ++depth)
     { d_center_of_mass_non_elastic_body[depth] = 0.0;}
	     
    for(int i =0; i < d_rigid_local_nodes; ++i)
     {
       for(int depth = 0; depth < NDIM; ++depth)
       {
	 d_center_of_mass_non_elastic_body[depth] += d_x_array(i,depth);
	 
       }       
     }
     
    //do the parallel reductions for center of mass.
    SAMRAI_MPI::sumReduction(&d_center_of_mass_non_elastic_body[0],NDIM);
    for(int depth = 0; depth < NDIM; ++depth)
      { d_center_of_mass_non_elastic_body[depth] /= d_rigid_global_nodes;}
      
    if( !SAMRAI_MPI::getRank() && d_print_output && d_output_COM_coordinates && d_INS_cycle_num == (d_INS_num_cycles -1) 
       && (d_timestep_counter % d_output_interval ) == 0)
     {
       d_position_COM_stream << getSimulationTime() << '\t' << d_center_of_mass_non_elastic_body[0] << '\t' 
                             << d_center_of_mass_non_elastic_body[1] << '\t' << d_center_of_mass_non_elastic_body[2] 
                             << std::endl;
     }  
    
    return;
    
  }//calculateCenterOfMassOfNonElasticBody



   void
   RigidityConstraint::calculateInertiaTensor()
   {
  
    /*!
     * Map the global result. The components are Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz.
     */

#if (NDIM == 2)

    blitz::Range ALL = blitz::Range::all();
    blitz::Range COLUMN1(0,0);
    blitz::Range COLUMN2(1,1);
   
    //Get the position vector from center of mass of the lagrangian structure.  
    blitz::Array<double,2>  xcol_offset(d_rigid_local_nodes,1), ycol_offset(d_rigid_local_nodes,1);
    xcol_offset   = d_x_array(ALL, COLUMN1) - d_center_of_mass_non_elastic_body[0];
    ycol_offset   = d_x_array(ALL, COLUMN2) - d_center_of_mass_non_elastic_body[1];
    
   // Calculate local moment of inertias.
    double Ixx_local = 0.0, Ixy_local = 0.0, Iyy_local = 0.0, Izz_local = 0.0;
    blitz::Array<double,2> i_xx( pow2(ycol_offset) ), i_xy( -1 * xcol_offset*ycol_offset ), i_yy( pow2(xcol_offset) ), 
                           i_zz( pow2(xcol_offset)  + pow2(ycol_offset) );

    for(int i = 0; i < d_rigid_local_nodes; ++i)
    {
           Ixx_local +=  i_xx(i,0);          // sum(m_i *y_i^2)
           Ixy_local +=  i_xy(i,0);          // -sum(m_i*x_i*y_i)
           Iyy_local +=  i_yy(i,0);          // sum(m_i*x_i^2)
           Izz_local +=  i_zz(i,0);          // sum(m_i*[x_i^2 + y_i^2])

    }

   // Perform the necessary reductions.
    d_map_inertiaTensor["Ixx"] = SAMRAI_MPI::sumReduction(Ixx_local);
    d_map_inertiaTensor["Ixy"] = SAMRAI_MPI::sumReduction(Ixy_local);
    d_map_inertiaTensor["Ixz"] = 0.0;
    d_map_inertiaTensor["Iyx"] = d_map_inertiaTensor["Ixy"];
    d_map_inertiaTensor["Iyy"] = SAMRAI_MPI::sumReduction(Iyy_local);
    d_map_inertiaTensor["Iyz"] = 0.0;
    d_map_inertiaTensor["Izx"] = 0.0;
    d_map_inertiaTensor["Izy"] = 0.0;
    d_map_inertiaTensor["Izz"] = SAMRAI_MPI::sumReduction(Izz_local);

#endif

#if (NDIM == 3)
   
    blitz::Range ALL = blitz::Range::all();
    blitz::Range COLUMN1(0,0);
    blitz::Range COLUMN2(1,1);
    blitz::Range COLUMN3(2,2);

    // Get the position vector from center of mass of the lagrangian structure.  
    blitz::Array<double,2>  xcol_offset(d_rigid_local_nodes,1), ycol_offset(d_rigid_local_nodes,1), zcol_offset(d_rigid_local_nodes,1);
    xcol_offset   = d_x_array(ALL, COLUMN1) - d_center_of_mass_non_elastic_body[0];
    ycol_offset   = d_x_array(ALL, COLUMN2) - d_center_of_mass_non_elastic_body[1];
    zcol_offset   = d_x_array(ALL, COLUMN3) - d_center_of_mass_non_elastic_body[2];


   // Calculate local moment of inertias.
    double Ixx_local = 0.0, Ixy_local = 0.0, Ixz_local = 0.0,  Iyy_local = 0.0, Iyz_local = 0.0, Izz_local = 0.0;

    blitz::Array<double, 2> i_xx( pow2(ycol_offset) + pow2(zcol_offset) ), i_xy( -1*xcol_offset*ycol_offset ), i_xz( -1*xcol_offset*zcol_offset ),
                 i_yy( pow2(xcol_offset) + pow2(zcol_offset) ), i_yz( -1*ycol_offset*zcol_offset ), i_zz( pow2(xcol_offset) + pow2(ycol_offset) );

    for(int i = 0 ; i < d_rigid_local_nodes; ++i)
    {
        Ixx_local += i_xx(i,0);
        Ixy_local += i_xy(i,0);
        Ixz_local += i_xz(i,0);
        Iyy_local += i_yy(i,0);
        Iyz_local += i_yz(i,0);
        Izz_local += i_zz(i,0);

    }

   // Perform the necessary reductions.
    d_map_inertiaTensor["Ixx"] = SAMRAI_MPI::sumReduction(Ixx_local);
    d_map_inertiaTensor["Ixy"] = SAMRAI_MPI::sumReduction(Ixy_local);
    d_map_inertiaTensor["Ixz"] = SAMRAI_MPI::sumReduction(Ixz_local);
    d_map_inertiaTensor["Iyx"] = d_map_inertiaTensor["Ixy"];
    d_map_inertiaTensor["Iyy"] = SAMRAI_MPI::sumReduction(Iyy_local);
    d_map_inertiaTensor["Iyz"] = SAMRAI_MPI::sumReduction(Iyz_local);
    d_map_inertiaTensor["Izx"] = d_map_inertiaTensor["Ixz"];
    d_map_inertiaTensor["Izy"] = d_map_inertiaTensor["Iyz"];
    d_map_inertiaTensor["Izz"] = SAMRAI_MPI::sumReduction(Izz_local);


#endif
    
    // write the Inertia Tensor to the file
   if( !SAMRAI_MPI::getRank() && d_print_output && d_output_MOI && d_INS_cycle_num == (d_INS_num_cycles -1) 
       && (d_timestep_counter % d_output_interval ) == 0 )
   {
     d_moment_of_inertia_stream << getSimulationTime()        << '\t' << d_map_inertiaTensor["Ixx"] << '\t'
                                << d_map_inertiaTensor["Ixy"] << '\t' << d_map_inertiaTensor["Ixz"] << '\t'
                                << d_map_inertiaTensor["Iyx"] << '\t' << d_map_inertiaTensor["Iyy"] << '\t'
                                << d_map_inertiaTensor["Iyz"] << '\t' << d_map_inertiaTensor["Izx"] << '\t'
                                << d_map_inertiaTensor["Izy"] << '\t' << d_map_inertiaTensor["Izz"] << std::endl;

   }
    

 return;

  }// calculateInertiaTensor



  void
  RigidityConstraint::solveRigidRotationalVelocity()
  {
    // Get the data related to position of lagrangian points from center of mass.
    blitz::Array<double,2> x_array_offset;
    x_array_offset.resize(d_rigid_local_nodes, NDIM);

#if (NDIM == 2)

    blitz::Range ALL = blitz::Range::all();
    blitz::Range COLUMN1(0,0);
    blitz::Range COLUMN2(1,1);

    x_array_offset(ALL, COLUMN1)  =  d_x_array(ALL, COLUMN1) - d_center_of_mass_non_elastic_body[0];
    x_array_offset(ALL, COLUMN2)  =  d_x_array(ALL, COLUMN2) - d_center_of_mass_non_elastic_body[1];

    double  r_cross_u_udef_local_z = 0.0;
    blitz::Array<double, 2> ang_momentum_z( x_array_offset(ALL, COLUMN1) * d_u_diff_array(ALL, COLUMN2) -
                                            x_array_offset(ALL, COLUMN2) * d_u_diff_array(ALL, COLUMN1) );

    for(int i = 0; i < d_rigid_local_nodes; ++i)
     { r_cross_u_udef_local_z += ang_momentum_z(i,0);}

    // Perform the necessary reductions.
    double r_cross_u_udef_global_z  =   SAMRAI_MPI::sumReduction(r_cross_u_udef_local_z);

    // Store the result.
    d_rigid_rot_vel[0]              =  0.0;
    d_rigid_rot_vel[1]              =  0.0;
    d_rigid_rot_vel[2]              =  r_cross_u_udef_global_z/d_map_inertiaTensor["Izz"];
  
#endif

#if (NDIM == 3)    

    blitz::Range ALL = blitz::Range::all();
    blitz::Range COLUMN1(0,0);
    blitz::Range COLUMN2(1,1);
    blitz::Range COLUMN3(2,2);

    x_array_offset(ALL, COLUMN1)  =  d_x_array(ALL, COLUMN1) - d_center_of_mass_non_elastic_body[0];
    x_array_offset(ALL, COLUMN2)  =  d_x_array(ALL, COLUMN2) - d_center_of_mass_non_elastic_body[1];
    x_array_offset(ALL, COLUMN3)  =  d_x_array(ALL, COLUMN3) - d_center_of_mass_non_elastic_body[2];

    double r_cross_u_udef_local_x = 0.0, r_cross_u_udef_local_y = 0.0, r_cross_u_udef_local_z = 0.0;

    blitz::Array<double,2> ang_momentum_x(  x_array_offset(ALL, COLUMN2 )  * d_u_diff_array(ALL,COLUMN3) 
                                          - x_array_offset(ALL, COLUMN3)   * d_u_diff_array(ALL, COLUMN2) ) ,
                                          
                           ang_momentum_y( -1*x_array_offset(ALL, COLUMN1) * d_u_diff_array(ALL, COLUMN3)
                                            + x_array_offset(ALL, COLUMN3) * d_u_diff_array(ALL, COLUMN1) ),
                                            
                           ang_momentum_z( x_array_offset(ALL, COLUMN1)    * d_u_diff_array(ALL,COLUMN2 )
                                           - x_array_offset(ALL, COLUMN2)  * d_u_diff_array(ALL, COLUMN1) ); 

    for(int i = 0; i < d_rigid_local_nodes; ++ i)
    {
      r_cross_u_udef_local_x += ang_momentum_x(i,0);
      r_cross_u_udef_local_y += ang_momentum_y(i,0);
      r_cross_u_udef_local_z += ang_momentum_z(i,0);
    }

    // Perform the necessary reductions.
    double r_cross_u_udef_global_x    =   SAMRAI_MPI::sumReduction(r_cross_u_udef_local_x);
    double r_cross_u_udef_global_y    =   SAMRAI_MPI::sumReduction(r_cross_u_udef_local_y);
    double r_cross_u_udef_global_z    =   SAMRAI_MPI::sumReduction(r_cross_u_udef_local_z);

  
    // Solve the 3X3 system of equations to get Omega_x, Omega_y & Omega_z.
    d_rigid_rot_vel                  =    solveSystemOfEqns(d_map_inertiaTensor,
                                                            r_cross_u_udef_global_x, 
                                                            r_cross_u_udef_global_y,
                                                            r_cross_u_udef_global_z);
    
    // Lock those components which are not needed.
    for(int dim = 0; dim < 3; ++dim)
     {
       if(!d_calculate_rotational_momentum[dim])
        { d_rigid_rot_vel[dim] = 0.0; }
     }
     
#endif

   // write rigid rotational velocity to the file
    if( !SAMRAI_MPI::getRank() && d_print_output && d_output_rot_vel && d_INS_cycle_num == (d_INS_num_cycles -1) 
         && (d_timestep_counter % d_output_interval ) == 0 )
     {
       d_rot_vel_stream << getSimulationTime() << '\t' << d_rigid_rot_vel[0]  << '\t'
                        << d_rigid_rot_vel[1]  << '\t' << d_rigid_rot_vel[2]  << '\t'
                        << d_omega_com_def[0]  << '\t' << d_omega_com_def[1]  << '\t'
                        << d_omega_com_def[2]  << std::endl;
                
     }

    return;

  }//solveRigidRotationalVelocity



  void
  RigidityConstraint::correctVelocityOnLagrangianMesh()
  {

    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
    const Pointer< LNodeLevelData > x_data      = d_lag_data_manager->getLNodeLevelData("X",d_finest_ln);
    std::vector<double> omegaCrossR(3,0.0), radiusFromCoM(3,0.0);
   
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = level->getPatch(p());
       Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
       const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
       const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
	
       for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
          const LNodeIndex& node_idx = *it;
          const int lag_index = node_idx.getLagrangianIndex();
          const int local_petsc_index   = node_idx.getLocalPETScIndex();

	  if(d_body_is_self_rotating)
	   {   
            // add rigidity constraint force in non elastic part of the body and zero out the constraint force
	    // in the elastic part of the body.
	     if ( lag_index < d_rigid_global_nodes ) // work in non-elastic domain begin
	      {
                for(int depth = 0; depth < NDIM ; depth++)
		 { radiusFromCoM[depth] = (*x_data)(local_petsc_index,depth) - d_center_of_mass_non_elastic_body[depth]; }
	            
		omegaCrossR[0] =  radiusFromCoM[2]*d_rigid_rot_vel[1] - radiusFromCoM[1]*d_rigid_rot_vel[2];
                omegaCrossR[1] = -radiusFromCoM[2]*d_rigid_rot_vel[0] + radiusFromCoM[0]*d_rigid_rot_vel[2];
                omegaCrossR[2] =  radiusFromCoM[1]*d_rigid_rot_vel[0] - radiusFromCoM[0]*d_rigid_rot_vel[1];
            
                //fill in the deformation velocity correction at the LNodeLevelData and scale it up for spreading operation.
	        for(int depth = 0; depth < NDIM; ++depth)
                 {   
		   if(d_kinematics_need_filtering)
		    {  
	              (*d_lag_vel_correction_data)(local_petsc_index,depth) = (d_rigid_trans_vel[depth] + omegaCrossR[depth] 
                                   + (d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_index][depth]
                                   - (*d_lag_vel_correction_data)(local_petsc_index,depth) 
									      )*d_volume_of_non_elastic_lag_pt;
										      									      
                      (*d_lag_new_vel_data)(local_petsc_index,depth)        = 0.5*(d_rigid_trans_vel[depth] + omegaCrossR[depth] 
                                   + (d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_index][depth] 
                                   + (*d_lag_old_vel_data)(local_petsc_index,depth)  );
		    }
		   else
		    { 
		      (*d_lag_vel_correction_data)(local_petsc_index,depth)  = (d_rigid_trans_vel[depth] + omegaCrossR[depth] 
                                   + (d_kinematics_ptr->d_map_kinematics_vel)[lag_index][depth]
                                   - (*d_lag_vel_correction_data)(local_petsc_index,depth) 
									       )*d_volume_of_non_elastic_lag_pt;
										     								      
                      (*d_lag_new_vel_data)(local_petsc_index,depth)         = 0.5*( d_rigid_trans_vel[depth] + omegaCrossR[depth] 
                                   + (d_kinematics_ptr->d_map_kinematics_vel)[lag_index][depth]  
                                   + (*d_lag_old_vel_data)(local_petsc_index,depth)   );
			      
		   }
                 }
	      } // work in non-elastic domain end
	     else // work in elastic domain.
              {   	                 
	        for(int depth = 0; depth < NDIM; ++depth)
                 {
	           (*d_lag_vel_correction_data)(local_petsc_index,depth) =  0.0;
                   (*d_lag_new_vel_data)(local_petsc_index,depth)        = 0.5*( (*d_lag_new_vel_data)(local_petsc_index,depth) 
                                                                                     + (*d_lag_old_vel_data)(local_petsc_index,depth) );
                 }
		       
	      }

	   } // body isrotating end...
	  else // body is not rotating begin...
	   {
	     if ( lag_index < d_rigid_global_nodes ) // work in non-elastic domain
              {
                for(int depth = 0; depth < NDIM; ++depth)
                 {	  
		   if(d_kinematics_need_filtering)
	            {
	              (*d_lag_vel_correction_data)(local_petsc_index,depth)  = (d_rigid_trans_vel[depth] 
	                        + (d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_index][depth]
                                - (*d_lag_vel_correction_data)(local_petsc_index,depth) 
									       )*d_volume_of_non_elastic_lag_pt;
										       										       
		      (*d_lag_new_vel_data)(local_petsc_index,depth)         = 0.5*(d_rigid_trans_vel[depth] 
                                + (d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_index][depth] 
                                + (*d_lag_old_vel_data)(local_petsc_index,depth)   );
	             }
	            else
	             {    
	               (*d_lag_vel_correction_data)(local_petsc_index,depth)  = (d_rigid_trans_vel[depth] 
	                       + (d_kinematics_ptr->d_map_kinematics_vel)[lag_index][depth]
	                       - (*d_lag_vel_correction_data)(local_petsc_index,depth) 
										 )*d_volume_of_non_elastic_lag_pt;
										       							       
		       (*d_lag_new_vel_data)(local_petsc_index,depth)         = 0.5*(d_rigid_trans_vel[depth] 
                               +  (d_kinematics_ptr->d_map_kinematics_vel)[lag_index][depth] 
                               +  (*d_lag_old_vel_data)(local_petsc_index,depth)       );
			     
		     }
                 }
	      }
             else //work in elastic domain.
	      {   
	        for(int depth = 0; depth < NDIM; ++depth)
                 {
	           (*d_lag_vel_correction_data)(local_petsc_index,depth)     =  0.0;
                   (*d_lag_new_vel_data)(local_petsc_index,depth) = 0.5*( (*d_lag_new_vel_data)(local_petsc_index,depth)  
                                                                        + (*d_lag_old_vel_data)(local_petsc_index,depth) 
									);
                 }
	      }

           }// body is not rotating end...
	       
        }// over a patch

      }//over patches

  return;

  }// correctVelocityOnLagrangianMesh


  void
  RigidityConstraint::calculateDragFromInertia(const double current_time, const double new_time)
  {
    
    const double dt = new_time - current_time;
    double kinetic_energy_local = 0.0;
    
    const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = level->getPatch(p());
       Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
       const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
       const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
       
       for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
          const LNodeIndex& node_idx = *it;
          const int local_petsc_index   = node_idx.getLocalPETScIndex();
	  const int lag_index           = node_idx.getLagrangianIndex();

	  if ( lag_index < d_rigid_global_nodes ) // in non-elastic domain
	   {
             for(int depth = 0; depth < NDIM; ++depth)
	      {
	        kinetic_energy_local += 0.5*d_volume_of_non_elastic_lag_pt*d_rho_fluid*pow( (*d_lag_new_vel_data)(local_petsc_index,depth),2 );
                (*d_lag_old_vel_data)(local_petsc_index,depth) = ( (*d_lag_new_vel_data)(local_petsc_index,depth)
                                                                  -(*d_lag_old_vel_data)(local_petsc_index,depth) ) * d_volume_of_non_elastic_lag_pt;
	      }
	   }
	  else // in elastic domain
	   {
	     for(int depth = 0; depth < NDIM; ++depth)
	      {
	        kinetic_energy_local += 0.5*d_volume_of_elastic_lag_pt*d_rho_fluid*pow( (*d_lag_new_vel_data)(local_petsc_index,depth),2 );
                (*d_lag_old_vel_data)(local_petsc_index,depth) = ( (*d_lag_new_vel_data)(local_petsc_index,depth)
                                                                -(*d_lag_old_vel_data)(local_petsc_index,depth) ) * d_volume_of_elastic_lag_pt;
	      }
	   }
	        
	}// over a patch
     }// over patches
      
    const double kinetic_energy_global = SAMRAI_MPI::sumReduction(kinetic_energy_local);

    // Get the data related to (U_new - U_old)dv_i.
    const Pointer< LNodeLevelData > U_diff      = d_lag_data_manager->getLNodeLevelData(d_object_name + "::lag_old_vel_data",d_finest_ln);
    const int  local_node_count_udiff           = U_diff->getLocalNodeCount();
    const Vec& petsc_vec_udiff                  = U_diff->getGlobalVec();
    const int  depth_udiff                      = U_diff->getDepth();
    blitz::Array<double,2>                        u_diff_array(local_node_count_udiff,depth_udiff);
    PetscScalar*                                  underlyingArray_udiff;

    // Get the underlying array from PETSc vec and store it in u_diff_array.
    VecGetArray(petsc_vec_udiff,&underlyingArray_udiff);
    int p = -1;
    for(int i = 0; i < local_node_count_udiff; ++i)
      {
        for(int j =0 ; j < depth_udiff ; ++j)
	  {
	    u_diff_array(i,j) = underlyingArray_udiff[++p];
	  }

      }   
    // Restore the array.
    VecRestoreArray(petsc_vec_udiff,&underlyingArray_udiff);

    std::vector<double> inertiaForce(NDIM,0.0);
    for(int j = 0; j < NDIM ; ++j)
    {
        for (int i = 0; i < u_diff_array.rows(); ++i)
        {
        	inertiaForce[j] +=  u_diff_array( i, j ) ;
        }
    }

    // Perform the sum reduction.
    SAMRAI_MPI::sumReduction(&inertiaForce[0],NDIM);    
    for(int i = 0; i < NDIM ; ++i)
     { inertiaForce[i] *= ( d_rho_fluid/dt );}



   //Get the data related to (U_cm + Omega X r + U_def - U_interp)dv_i for the solid domain.   
    std::vector<double> rigConstraintForce(NDIM,0.0);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = level->getPatch(p());
       const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
       const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
	
       for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
          const LNodeIndex& node_idx    = *it;
          const int lag_index           = node_idx.getLagrangianIndex();
          const int local_petsc_index   = node_idx.getLocalPETScIndex();
          if ( lag_index < d_rigid_global_nodes )
	   {     
           // fill in the deformation velocities at the LNodeLevelData
	     for(int depth = 0; depth < NDIM ; depth++)
	      { rigConstraintForce[depth] += (*d_lag_vel_correction_data)(local_petsc_index,depth); }
		  
	   }
	}// for a patch       
     }// for all patches.
     
    // Perform the sum reduction.
    SAMRAI_MPI::sumReduction(&rigConstraintForce[0],NDIM);    
    for(int i = 0; i < NDIM ; ++i)
    { rigConstraintForce[i] *= ( d_rho_fluid/dt );}

   std::vector<double> dragForce(3,0.0);
   for(int i = 0; i < NDIM; ++i)
   { dragForce[i] = inertiaForce[i] - rigConstraintForce[i]; }

   // write drag force to the file
   if( !SAMRAI_MPI::getRank() && d_print_output && (d_timestep_counter % d_output_interval ) == 0 )
   {
     d_drag_force_stream     << getSimulationTime() << '\t' << inertiaForce[0] << '\t'
                             << inertiaForce[1] << '\t' << inertiaForce[2] << '\t' << rigConstraintForce[0] << '\t'
                             << rigConstraintForce[1] << '\t'<< rigConstraintForce[2] << std::endl;
		       
     d_kinetic_energy_stream << getSimulationTime() << '\t' << kinetic_energy_global
                             << std::endl;
   }

   return;

  }// calculateDragFromInertia




  void
  RigidityConstraint::ensureMomentumConservation()
  {

    std::vector<double> u_sum(NDIM,0.0);
       
    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
     Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
     
     for (PatchLevel<NDIM>::Iterator p(level); p; p++)
     {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
	
        for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
         {
            const LNodeIndex& node_idx    = *it;
            const int lag_index           = node_idx.getLagrangianIndex();
            const int local_petsc_index   = node_idx.getLocalPETScIndex();
    
	    if ( lag_index < d_rigid_global_nodes )
	    {
	         
              // fill in the deformation velocities at the LNodeLevelData
	         for(int depth = 0; depth < NDIM ; depth++)
		  {
		     u_sum[depth] += (*d_lag_vel_correction_data)(local_petsc_index,depth);
		  }

	    }

	 }// for a patch

       
     }// for all patches.
     


    // Perform the reduction.
    SAMRAI_MPI::sumReduction(&u_sum[0],NDIM);
    for(int i = 0; i < NDIM ; ++i)
	u_sum[i] /= d_rigid_global_nodes;

  

    // Subtract the u_sum from the lagrangian velocity in the non-elastic domain.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
     {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
         {
              const LNodeIndex& node_idx = *it;            
	      const int lag_index           = node_idx.getLagrangianIndex();
              const int local_petsc_index   = node_idx.getLocalPETScIndex();
	      
	    if ( lag_index < d_rigid_global_nodes )
	    {
	      
              // Fill in the deformation velocity correction at the LNodeLevelData.
	      for(int depth = 0; depth < NDIM ; depth++)
		  (*d_lag_vel_correction_data)(local_petsc_index,depth) -= u_sum[depth];
	      
	    }

	 }// over a patch


     }// over patches

   return;
 
  }// ensureMomentumConservation


  void
  RigidityConstraint::spreadCorrectedLagrangianVelocity()
  {
       
 /*!
   * The following data members are used in spreading/interpolating the deformation velocities 
   * from the lagrangian grid to the Eulerian grid and vice versa.The vector should have all the
     components for each level of AMR grid.
   */
   std::vector< Pointer< LNodeLevelData > > F_data(d_finest_ln+1,  Pointer<LNodeLevelData>(NULL) );
   std::vector< Pointer< LNodeLevelData > > X_data(d_finest_ln+1,  Pointer<LNodeLevelData>(NULL) );


   // Fill in the above vectors at the finest level.
    F_data[d_finest_ln] = d_lag_vel_correction_data;   
    X_data[d_finest_ln] = d_lag_data_manager->getLNodeLevelData("X",d_finest_ln);
 
  // Spread the deformation velocities.
   d_lag_data_manager->spread(d_U_fluidSolve_idx ,
                              F_data,
			      X_data,
                              std::vector< Pointer< RefineSchedule< NDIM > > > (),
			      true,
			      true,
			      d_finest_ln,
			      d_finest_ln);

   return;

  }// spreadCorrectedLagrangianVelocity

 void 
 RigidityConstraint::synchronizeLevels()
 {

    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it!= d_calgs.end(); ++it)
     {
       d_cscheds[(*it).first].resize(d_finest_ln+1);
     }

    if(d_finest_ln > 0)
     {
       // (Re)build coarsen communication schedules.  These are set only for levels >= 1.
       for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
                 it != d_calgs.end(); ++it)
        {
          for (int ln = std::max(d_coarsest_ln,1); ln <= d_finest_ln; ++ln)
           {
             Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(ln);
             Pointer<PatchLevel<NDIM> > coarser_level = d_patch_hierarchy->getPatchLevel(ln-1);
             d_cscheds[(*it).first][ln] = (*it).second->createSchedule(coarser_level, level, d_cstrategies[(*it).first]);
           }
       }

     }
     
    // Do the Coarsening Operation.
    for(int ln = d_finest_ln; ln > d_coarsest_ln; --ln)
    {
      d_cscheds["SYNCH_EULERIAN_VEL_CORRECTION"][ln]->coarsenData();
    }

    return;

 }// synchronizeLevels
 
 
 
 void
 RigidityConstraint::updateLagrangianMarkerPosition(const double current_time, const double new_time)
 {
   
   const double dt = new_time - current_time;
   static std::vector<double> x_center_of_mass_old(3);
   
   if(d_INS_cycle_num == 0)
    {
      x_center_of_mass_old = d_center_of_mass_non_elastic_body;
      d_kinematics_ptr->updateLagrangianMarkerShapeAndOrientation(d_incremented_angle_from_reference_axis, new_time);
    }
   
   // Get the patch data descriptor index for the LNodeIndexData.
   const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
   Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
   for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
      Pointer<Patch<NDIM> > patch = level->getPatch(p());
      const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
      const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
      for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
       {
         const LNodeIndex& node_idx    = *it;
         const int lag_index           = node_idx.getLagrangianIndex();
         const int local_petsc_index   = node_idx.getLocalPETScIndex();

         // fill in the new shape at the LNodeLevelData
         if( (lag_index < d_rigid_global_nodes))
          {
	    for(int depth = 0; depth < NDIM; ++depth)
	     { 
               if(d_kinematics_need_filtering)
		{
	          (*d_lag_new_position_data)(local_petsc_index,depth) = (d_kinematics_ptr->d_map_new_shape)[lag_index][depth]
	                      + x_center_of_mass_old[depth] - d_vel_com_def_old[depth]*dt
	                      + 0.5*(d_rigid_trans_vel[depth] + d_rigid_trans_vel_old[depth] )*dt;
		}
	       else
		{
		  (*d_lag_new_position_data)(local_petsc_index,depth) = (d_kinematics_ptr->d_map_new_shape)[lag_index][depth]
	                     + x_center_of_mass_old[depth] 
	                     + 0.5*(d_rigid_trans_vel[depth] + d_rigid_trans_vel_old[depth] )*dt;
			     										   
		} 
	      }							       
	  }
       }// over a patch
    }// over all the patches on the finest level.
   
   
 }//updateLagrangianMarkerPosition

 
  void
  RigidityConstraint::applyProjection()
  {

    const int coarsest_ln = 0;
    const int finest_ln   = d_patch_hierarchy->getFinestLevelNumber();

    // Allocate temporary data.
    ComponentSelector scratch_idxs;
    scratch_idxs.setFlag(d_U_scratch_idx);
    scratch_idxs.setFlag(d_Phi_idx);
    scratch_idxs.setFlag(d_Div_U_scratch_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, getSimulationTime());
    }

    // Compute div U before applying the projection operator.
    const bool U_current_cf_bdry_synch = true;
    d_hier_math_ops->div(
        d_Div_U_scratch_idx, d_Div_U_var,                // dst
        +1.0,                                            // alpha
        d_U_fluidSolve_idx, d_U_fluidSolve_var,          // src
        d_no_fill_op,                                    // src_bdry_fill
        getSimulationTime(),                               // src_bdry_fill_time
        U_current_cf_bdry_synch);                        // src_cf_bdry_synch


    if (d_do_log)
    {
        const double Div_U_norm_1  = d_hier_cc_data_ops->L1Norm( d_Div_U_scratch_idx, d_wgt_cc_idx);
        const double Div_U_norm_2  = d_hier_cc_data_ops->L2Norm( d_Div_U_scratch_idx, d_wgt_cc_idx);
        const double Div_U_norm_oo = d_hier_cc_data_ops->maxNorm(d_Div_U_scratch_idx, d_wgt_cc_idx);
        plog << d_object_name << "::applyProjection():\n"
             << "  performing velocity correction projection\n"
             << "  before projection:\n"
             << "    ||Div U||_1  = " << Div_U_norm_1  << "\n"
             << "    ||Div U||_2  = " << Div_U_norm_2  << "\n"
             << "    ||Div U||_oo = " << Div_U_norm_oo << "\n";
    }

    // Setup the solver vectors.
    d_hier_cc_data_ops->setToScalar(d_Phi_idx, 0.0, false);
    d_hier_cc_data_ops->scale(d_Div_U_scratch_idx, -1.0, d_Div_U_scratch_idx);
    const double Div_U_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(d_Div_U_scratch_idx, d_wgt_cc_idx);
    d_hier_cc_data_ops->addScalar(d_Div_U_scratch_idx, d_Div_U_scratch_idx, -Div_U_mean);

    SAMRAIVectorReal<NDIM,double> sol_vec(d_object_name+"::sol_vec", d_patch_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(d_Phi_var, d_Phi_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM,double> rhs_vec(d_object_name+"::rhs_vec", d_patch_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_Div_U_var, d_Div_U_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    // Setup the Poisson solver.
    d_velcorrection_projection_spec->setCZero();
    d_velcorrection_projection_spec->setDConstant(-1.0);

    d_velcorrection_projection_op->setPoissonSpecifications(*d_velcorrection_projection_spec);
    d_velcorrection_projection_op->setPhysicalBcCoef(&d_velcorrection_projection_bc_coef);
    d_velcorrection_projection_op->setHomogeneousBc(true);
    d_velcorrection_projection_op->setTime(getSimulationTime());
    d_velcorrection_projection_op->setHierarchyMathOps(d_hier_math_ops);

    d_velcorrection_projection_fac_op->setPoissonSpecifications(*d_velcorrection_projection_spec);
    d_velcorrection_projection_fac_op->setPhysicalBcCoef(&d_velcorrection_projection_bc_coef);
    d_velcorrection_projection_fac_op->setTime(getSimulationTime());

    d_velcorrection_projection_solver->setInitialGuessNonzero(false);
    d_velcorrection_projection_solver->setOperator(d_velcorrection_projection_op);

    // Solve the projection Poisson problem.
    d_velcorrection_projection_solver->initializeSolverState(sol_vec,rhs_vec);
    d_velcorrection_projection_solver->solveSystem(sol_vec,rhs_vec);
    d_velcorrection_projection_solver->deallocateSolverState();

    // NOTE: We always use homogeneous Neumann boundary conditions for the
    // velocity correction projection Poisson solver.
    d_velcorrection_projection_solver->setNullspace(true, NULL);

    // Setup the interpolation transaction information.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_bc_component(d_Phi_idx, CELL_DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY,  
                                                      &d_velcorrection_projection_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> Phi_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    Phi_bdry_bc_fill_op->initializeOperatorState(Phi_bc_component, d_patch_hierarchy);

    // Fill the physical boundary conditions for Phi.
    Phi_bdry_bc_fill_op->setHomogeneousBc(true);
    Phi_bdry_bc_fill_op->fillData(getSimulationTime());

    // Set U := U - grad Phi.
    const bool U_scratch_cf_bdry_synch = true;
    d_hier_math_ops->grad(
        d_U_scratch_idx, d_U_var,  // dst
        U_scratch_cf_bdry_synch,   // dst_cf_bdry_synch
        1.0,                       // alpha
        d_Phi_idx, d_Phi_var,      // src
        d_no_fill_op,              // src_bdry_fill
        getSimulationTime());        // src_bdry_fill_time
    d_hier_sc_data_ops->axpy(d_U_fluidSolve_idx, -1.0, d_U_scratch_idx, d_U_fluidSolve_idx);

    // Compute div U after applying the projection operator
    if (d_do_log)
    {
        const bool U_current_cf_bdry_synch = true;
        d_hier_math_ops->div(
            d_Div_U_scratch_idx, d_Div_U_var,                // dst
            +1.0,                                            // alpha
            d_U_fluidSolve_idx, d_U_fluidSolve_var,          // src
            d_no_fill_op,                                    // src_bdry_fill
            getSimulationTime(),                               // src_bdry_fill_time
            U_current_cf_bdry_synch);                        // src_cf_bdry_synch
        const double Div_U_norm_1  = d_hier_cc_data_ops->L1Norm( d_Div_U_scratch_idx, d_wgt_cc_idx);
        const double Div_U_norm_2  = d_hier_cc_data_ops->L2Norm( d_Div_U_scratch_idx, d_wgt_cc_idx);
        const double Div_U_norm_oo = d_hier_cc_data_ops->maxNorm(d_Div_U_scratch_idx, d_wgt_cc_idx);
        plog << "  after projection:\n"
             << "    ||Div U||_1  = " << Div_U_norm_1  << "\n"
             << "    ||Div U||_2  = " << Div_U_norm_2  << "\n"
             << "    ||Div U||_oo = " << Div_U_norm_oo << "\n";
    }

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }
    return;

  }// applyProjection


  
  RigidityConstraint::~RigidityConstraint()
  {

     delete d_velcorrection_projection_spec;
     
  }// ~RigidityConstraint


  void 
  RigidityConstraint::postprocessData(const double current_time, const double new_time )
  {
    
    if(d_output_drag_and_kinetic_energy)
    { calculateDragFromInertia(current_time, new_time); }
    
    if(d_output_power)
    { calculatePower(current_time,new_time); }

   // U_old = U_new for both elatic and non-elastic domain.
   const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
   Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
   for (PatchLevel<NDIM>::Iterator p(level); p; p++)
   {
     Pointer<Patch<NDIM> > patch = level->getPatch(p());        
     const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
     const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
     for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
      {
        const LNodeIndex& node_idx = *it;
        const int local_petsc_index   = node_idx.getLocalPETScIndex();

        for(int depth = 0; depth < NDIM; ++depth)
	 { (*d_lag_old_vel_data)(local_petsc_index,depth) =  (*d_lag_new_vel_data)(local_petsc_index,depth); }

      }
   }
    
   return;
    
  }//postprocessData


  void
  RigidityConstraint::calculatePower(const double current_time, const double new_time)
  {
    
    const double dt  = new_time - current_time;
    Pointer< PatchLevel<NDIM> > finest_level = d_patch_hierarchy->getPatchLevel(d_finest_ln);
    const int lag_node_index_idx = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
      
    // allocate data for tagging cells and calculating power.
    if( finest_level->checkAllocated(d_EulTagDefVel_scratch_idx) )
     { finest_level->deallocatePatchData(d_EulTagDefVel_scratch_idx);}
      
    if( finest_level->checkAllocated(d_EulDefVel_scratch_idx) )
     { finest_level->deallocatePatchData(d_EulDefVel_scratch_idx); }
      
    if( finest_level->checkAllocated(d_VisDefPower_scratch_idx) )
     { finest_level->deallocatePatchData(d_VisDefPower_scratch_idx); }
      
    if( finest_level->checkAllocated(d_ConstraintPower_scratch_idx) )
     { finest_level->deallocatePatchData(d_ConstraintPower_scratch_idx); }
	
    finest_level->allocatePatchData(d_EulTagDefVel_scratch_idx    , new_time);
    finest_level->allocatePatchData(d_EulDefVel_scratch_idx       , new_time);
    finest_level->allocatePatchData(d_VisDefPower_scratch_idx     , new_time);
    finest_level->allocatePatchData(d_ConstraintPower_scratch_idx , new_time);
     
    // initialize variables. 
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
       Pointer<CellData<NDIM,int   > > EulTagDefVel_scratch_Data      = patch->getPatchData(d_EulTagDefVel_scratch_idx);
       Pointer<CellData<NDIM,double> > EulDefVel_scratch_Data         = patch->getPatchData(d_EulDefVel_scratch_idx);
       Pointer<CellData<NDIM,double> > VisDefPower_scratch_Data       = patch->getPatchData(d_VisDefPower_scratch_idx);
       Pointer<CellData<NDIM,double> > ConstraintPower_scratch_Data   = patch->getPatchData(d_ConstraintPower_scratch_idx);
	 
       EulTagDefVel_scratch_Data->fill(0,EulTagDefVel_scratch_Data->getGhostBox());
       EulDefVel_scratch_Data->fillAll(0.0,EulDefVel_scratch_Data->getGhostBox());
       VisDefPower_scratch_Data->fillAll(0.0,VisDefPower_scratch_Data->getBox() );
       ConstraintPower_scratch_Data->fillAll(0.0, ConstraintPower_scratch_Data->getBox());
     }

    // Tag cells and set velocity in Eulerian cells.
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
       Pointer<CellData<NDIM,int   > > EulTagDefVel_scratch_Data      = patch->getPatchData(d_EulTagDefVel_scratch_idx);
       Pointer<CellData<NDIM,double> > EulDefVel_scratch_Data         = patch->getPatchData(d_EulDefVel_scratch_idx);
       Pointer<CellData<NDIM,double> > ConstraintPower_scratch_Data   = patch->getPatchData(d_ConstraintPower_scratch_idx);
	    
       Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
       const double* const dx                       = pgeom->getDx(); 
       const double* const XLower                   = pgeom->getXLower();
       const Box<NDIM>& patch_box                   = patch->getBox();
       const Index<NDIM>& ilower                    = patch_box.lower();					   
       const Pointer<LNodeIndexData> idx_data       = patch->getPatchData(lag_node_index_idx);	
	    
       for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
	  const LNodeIndex& node_idx    = *it;
	  const int lag_index           = node_idx.getLagrangianIndex();
	  const int local_petsc_index   = node_idx.getLocalPETScIndex();
	  const double* const X         = node_idx.getNodeLocation();
		    
	  if(lag_index < d_rigid_global_nodes )
	   {
             const CellIndex<NDIM> Lag2Eul_cellindex( Index<NDIM> ( int( floor((X[0]-XLower[0])/dx[0]))+ilower(0)
#if (NDIM > 1)
                                                 ,int(floor((X[1]-XLower[1])/dx[1]))+ilower(1)
#if (NDIM > 2)
                                                 ,int(floor((X[2]-XLower[2])/dx[2]))+ilower(2)
#endif
#endif
                                                                  )
					           );
		    
            (*EulTagDefVel_scratch_Data)(Lag2Eul_cellindex)++;
		   
	    for(int depth = 0; depth < NDIM; ++depth)
	     {
	       if(d_kinematics_need_filtering)
	        { (*EulDefVel_scratch_Data)(Lag2Eul_cellindex,depth) +=  (d_kinematics_ptr->d_map_filtered_kinematics_vel)[lag_index][depth];}
	       else
	        { (*EulDefVel_scratch_Data)(Lag2Eul_cellindex,depth) +=  (d_kinematics_ptr->d_map_kinematics_vel)[lag_index][depth]; }
		      
		(*ConstraintPower_scratch_Data)(Lag2Eul_cellindex,depth) +=  (d_rho_fluid/dt)*(*d_lag_vel_correction_data)(local_petsc_index,depth)
		                                                           * (*d_lag_new_vel_data)(local_petsc_index,depth);
	     }
		    
	   } // only in the non elastic domain.
		    
	}// on a patch
	    
       // Now iterate over the patch box and take the average of velocity in each cell
       for(CellData<NDIM,double>::Iterator cell_itr(patch_box); cell_itr; cell_itr++)
        {
	  const int num_pts_in_cell = (*EulTagDefVel_scratch_Data)(*cell_itr);
	  if( num_pts_in_cell )
	   {
	     for(int depth = 0; depth < NDIM; ++depth)
	      { (*EulDefVel_scratch_Data)(*cell_itr,depth)       /= num_pts_in_cell;}
	   }
	}// on the same patch
	    
     }// on all the patches.
      
    
    // copy the data from the neighboring pathes into the ghost region of the patch.
    // use refine schedule to do the job.
    Pointer<RefineSchedule<NDIM> > fill_ghost_schedule = d_PowerRefineAlgorithm.createSchedule(finest_level,finest_level);
    fill_ghost_schedule->fillData(new_time);
    
    // calculate power for each patch on the finest level.
    std::vector<double> rigConstraintPower(3,0.0), viscdefPower(3,0.0);
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
     {
       Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
       Pointer<CellData<NDIM,int   > > EulTagDefVel_scratch_Data    = patch->getPatchData(d_EulTagDefVel_scratch_idx);
       Pointer<CellData<NDIM,double> > EulDefVel_scratch_Data       = patch->getPatchData(d_EulDefVel_scratch_idx);
       Pointer<CellData<NDIM,double> > VisDefPower_scratch_Data     = patch->getPatchData(d_VisDefPower_scratch_idx);
       Pointer<CellData<NDIM,double> > ConstraintPower_scratch_Data = patch->getPatchData(d_ConstraintPower_scratch_idx);
	    
       Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
       const double* const dx                       = pgeom->getDx(); 
       const Box<NDIM>& patch_box                   = patch->getBox();
       const Index<NDIM>& ilower                    = patch_box.lower(); 
       const Index<NDIM>& iupper                    = patch_box.upper();
    
       CALCULATE_VISC_DEF_POWER_FC(EulTagDefVel_scratch_Data->getPointer(),
			  EulDefVel_scratch_Data->getPointer(),
	                  VisDefPower_scratch_Data->getPointer(),
			  d_mu_fluid,	    
			  1,
			  ilower(0), iupper(0),
			  ilower(1), iupper(1),
#if (NDIM == 3)
			  ilower(2), iupper(2),
#endif
			  dx);	
	
       //Iterate over this patch box to compute power.
      for(CellData<NDIM, double>::Iterator citr(patch_box); citr ; citr++)
       {
	 const int num_pts_in_cell = (*EulTagDefVel_scratch_Data)(*citr);
	 if( num_pts_in_cell)
	  {
            for(int depth = 0; depth < NDIM; ++depth)
	     {
	       rigConstraintPower[depth] += (*ConstraintPower_scratch_Data)(*citr, depth);
	       viscdefPower[depth]       += (*VisDefPower_scratch_Data)(*citr,depth) * num_pts_in_cell;
	     }
	      
	  }
       }
	
     }// on all patches.
     
     // Do the parallel reductions to sum the power on all patches.
     SAMRAI_MPI::sumReduction(&rigConstraintPower[0],3);
     SAMRAI_MPI::sumReduction(&viscdefPower[0], 3);
     
     for(int depth = 0 ; depth < NDIM; ++depth)
     { viscdefPower[depth]     *= d_finest_cell_vol;}
	
     
     // write to file.
     if(!SAMRAI_MPI::getRank() && d_print_output && (d_timestep_counter % d_output_interval ) == 0 )
     {
       
        d_power_spent_stream << new_time << "\t\t" << rigConstraintPower[0] << "\t\t" 
	                     << rigConstraintPower[1] << "\t\t" << rigConstraintPower[2] << "\t\t"
			     << viscdefPower[0] << "\t\t" << viscdefPower[1] << "\t\t" << viscdefPower[2] << std::endl;
       
     }


    return;
    
    
  }//calculatePower




} // namespace IBAMR
