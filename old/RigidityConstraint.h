
/******************************************************************************************************************************
   Filename: RigidityConstraint.h
   Written by amneet bhalla on Taylor@mech.northwestern.edu
   Created on 09/30/2011.


   This is a concrete class which implements fast and efficient distributed Lagrange multiplier (fictitious domain) method
   for immersed structure having density as that of surrounding fluid.
     
   References 
   1. Patankar et al. A new formulation of the distributed Lagrange multiplier/fictitious domain method 
      for particulate flows. Int. Journal of Multiphase flows, 26, 1509-1524 (2000).
   2. Shirgaonkar et al. A new mathematical formulation and fast algorithm for fully resolved simulation of self-propulsion.
      JCP, 228 , 2366-2390 (2009).
 
*******************************************************************************************************************************/

/////////////////////// INCLUDE GUARD //////////////////////////////////////////////

#ifndef included_RigidityConstraint
#define included_RigidityConstraint

///////////////////////////// INCLUDES /////////////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Array.h>


// IBAMR INCLUDES
#include "IBKinematics.h"
#include "RigidIBStaggeredHierarchyIntegrator.h"
#include <tbox/Database.h>

// IBTK INCLUDES
#include <ibtk/HierarchyGhostCellInterpolation.h>

// C++ STDLIB INCLUDES
#include <map>
#include <vector>
#include <string>
#include <fstream>

////////////////////////////// CLASS DEFINITION ////////////////////////////////////


namespace IBAMR
{



  /*!
   * Post processing call back function to be hooked into IBAMR::INSStaggeredHierachyIntegrator class.
   * 
   * \param current_time is the time at t_n.
   * \param new_time is the time at t_n+1 = t_n + dt.
   * \param cycle_num is the cycle of predictor-corrector scheme.
   * \param ctx is the pointer to IBAMR::RigidityConstraint class object.
   */
  
  

  void callRigidityConstraintCallBackFunction(const double current_time,
					              const double new_time,
                                                      const int cycle_num, 
						      void* ctx); 

  ///////////////////////////////////////////////////////////////////////////////////////////////

  /*!
   * \brief RigidityConstraint class.
   * 
   * This class provides implementation of distributed Lagrange multiplier force field
   * theory to model the transport of solids(rigid + deforming) in IBAMR framework.
   */

class RigidityConstraint
{

  /////////////////////////////// PUBLIC MEMBER FUNCTIONS ///////////////////////////////

  public:

  /*! 
   * Constructor.
   * Initialize the lagrangian deformation velocity data, hierarchy operation object, 
   * variable contexts and variables.
   */
   
  
  RigidityConstraint( const std::string& object_name,
                      SAMRAI::tbox::Pointer< IBAMR::RigidIBStaggeredHierarchyIntegrator >  ib_hier_integrator,
                      SAMRAI::tbox::Pointer< IBAMR::INSStaggeredHierarchyIntegrator > ins_hier_integrator,
		      SAMRAI::tbox::Pointer<IBAMR::IBKinematics> ib_kinematics_ptr,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

  
  /*!
   * Return the method used to update position of Lagrangian nodes.
   */
  const std::string
  getNameOfLagrangianPositionUpdateMethod();
   
  /*!
   * Return the name of LNodeLevelData object which stores the user defined velocity or position
   * of updated Lagrangian point.
   */
  const std::string
  getNameOfLagrangianPositionUpdateDataStructure();
  
   /*!
   * set the time at which Rigidity Constraint is to be applied in RigidityConstraint object.
   */
  void
  setSimulationTime(const double new_time);

  /*!
   * get the simulation time at which Rigidity Constraint was applied.
   */
  double
  getSimulationTime() const;

  /*!
   * Set hierarchy related stuff in RigidityConstraint object. This has to done after
   * every fluid solve.This method should be called at each time the Rigidity Constraint 
   * is invoked as the hierarchy might have changed and the hierarchy object might have been
   * destroyed during course of simulation.
   */
  void
  setHierarchyRelatedData(const int cycle_num); 
  
   /*!
    * calculate center of mass and moment of inertia of the Non Elastic domain
    * when the body is rotating.
    */
   void
   calculateCOMandMOIOfNonElasticBody();
   
   /*!
    * get volume associated with non elastic material point.
    */
    double
    getVolumeOfNonElasticLagrangianPoint() const
    {
      return d_volume_of_non_elastic_lag_pt;
    }//getVolumeOfNonElasticLagrangianPoint
    
   /*!
    * get volume associated with elastic material point.
    */
    double
    getVolumeOfElasticLagrangianPoint() const
    {
      return d_volume_of_elastic_lag_pt;
    }//getVolumeOfElasticLagrangianPoint   

  /*!
   * calculate & store kinematics velocities on different lagrangian points in a map object.
   */
  void 
  calculateKinematicsVelocity(const double current_time, const double new_time);

  /*!
   * Interpolate fluid solve velocity from Eulerian grid onto the lagrangian mesh.
   */
  void
  interpolateFluidSolveVelocity();

  /*!
   * Subtract the calculated deformation velocities from the interpolated fluidSolve 
   * velocity on the LNodeLevelData.
   */
  void
  subtractDeformationVelocity();

  /*!
   * Solve for rigid body translation velocity.
   */
  void 
  solveRigidTranslationalVelocity();

  /*!
   * Check if rigid body is self rotating or not?
   */
  bool
  bodyIsSelfRotating() const
  {
    return d_body_is_self_rotating;
    
  }//bodyIsSelfRotating

 /*!
  * Check if body is self translating or not?
  */
  bool
  bodyIsSelfTranslating() const
  {
    return d_body_is_self_translating;
    
  }// bodyIsSelfTranslating
  
  /*!
   * get no of INS cycles.
   */
  int
  getNoOfINSCycles() const
  {
    return d_INS_num_cycles;
    
  }//getNoOfINSCycles

  /*!
   * Solve for rigid rotational velocity.
   */
  void
  solveRigidRotationalVelocity();

  /*!
   * Correct solid velocity on Lagrangian mesh. Set the velocity on Lagrangian 
   * mesh as U_lag_corr = U_trans + Omega X r + U_def - U_interpolated.
   */
  void
  correctVelocityOnLagrangianMesh();

 
  /*!
   * Total momentum on the Eulerian mesh should be conserved. To make sure of this we sprinkle 
   * (-)U_lag_corr_avg on each material point.
   */
   void
   ensureMomentumConservation();
 
  /*!
   *  Spread the momentum conserved corrected velocities(U_lag_corr_mom_conserved) at the
   *  Lagrangian mesh to the Eulerian Grid. The lagrangian deformation velocities are stored 
   *  in workable d_lag_vel_data and it needs to be spread to the Eulerian Grid on variable 
   *  d_U_fluidSolve_var.
   */
  void 
  spreadCorrectedLagrangianVelocity();

  /*!
   * if div free projection step is needed?
   */
  bool
  needs_DivFreeProjection() const 
  {
    return d_needs_div_free_projection;
    
  }//needs_DivFreeProjection
    

  /*!
   * Synch the invalid regions of the levels which are smaller than finest_level in the hierarchy.
   */
  void 
  synchronizeLevels();

  /*!
   * The correction on Eulerian Grid can lead to kinkiness in the velocity field. To remove such kinkiness
   * we project the corrected velocity field onto a divergence free field.
   */
  void
  applyProjection();
  
  void
  updateLagrangianMarkerPosition(const double current_time, const double new_time);

   
   /*!
    * Calculate Drag, Kinetic Energy, Power.
    */
   void 
   postprocessData(const double current_time, const double new_time);

  /*!
   * Close the files in the destructor. Return the resources
   */ 
   ~RigidityConstraint();

  

  private:

  /////////////////////////////// PRIVATE MEMBER FUNCTIONS ///////////////////////////////

  /*!
   * The default constructor is not implemented and should not be used.
   */
  RigidityConstraint();

  /*!
   * Copy Constructor should not be used.
   */
  RigidityConstraint(const RigidityConstraint& from);
     
  /*!
   * Assignment operator should not be used.
   */
   RigidityConstraint&
   operator = (const RigidityConstraint& that);

  /*!
   * Get values from input database.
   */
   void
   getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);
   

  /*!
   * Set the velocity at the d_lag_vel_data_old to zero. The initial velocity at lag points is zero.
   */
   void
   setInitialLagrangianVelocity();
   
   /*!
    * Calculate center of mass for the non-elastic part of the body.
    */
    void
    calculateCenterOfMassOfNonElasticBody();
    
  /*!
   * calculate the Inertia Tensor of the rigid body at the current configuration.
   * The components are Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz.
   */
  void
  calculateInertiaTensor();
   
  /*!
    * calculate translational and angular momentum introduced in the system by the 
    * kinematic velocity.
    */
   void
   calculateMomentumOfKinematicsVelocity();
   
   /*!
    * filter the given deformational kinematics.
    */
   void
   filterGivenKinematicsVelocity();

   /*!
    * Calculate the volume of a Lagrangian material point in elastic and non-elastic domains.
    */
   void
   calculateVolOfLagPointInNonElasticDomain();
   
   void
   calculateVolOfLagPointInElasticDomain();
   
   /*!
    * Get the background mesh velocity on the Lagrangian points.
    */
   void
   getBackgroundFluidVelocity(SAMRAI::tbox::Pointer< IBTK::LNodeLevelData > lag_data, 
			      const double time);
   


  /*!
   * Calculate the drag Force on the body as
   * F_drag = rho*dU/dt  - rho*(U_rig - U_interp)/dt
   */
  void
  calculateDragFromInertia(const double current_time, const double new_time);
  
  /*!
   *  calculate power spent during the swimming motion. The power spent comprises of two parts:
   * (a) RigidityConstraint power := rho*(U_b - U_interp)/dt 
   * (b) Viscous Deformational power := 2*mu D(U_def):D(U_def)
   * 
   * Power spent during swimming is := RigidityConstraint power - Viscous Deformational power
   */
  void
  calculatePower(const double current_time, const double new_time);
    
    
  ////////////////////////////////  PRIVATE DATA ////////////////////////////////////////  
    
  /*!
   *  Object Name
   */
  std::string d_object_name;
  
  /*!
   * Pointer to RigidIBStaggeredHierarchyIntegrator and INSStaggeredHierarchyIntegrator,
   * which is the time integrators associated with immersed body and navier stokes
   * solver respectively.
   */
   SAMRAI::tbox::Pointer< IBAMR::RigidIBStaggeredHierarchyIntegrator >  d_ib_hier_integrator;
   SAMRAI::tbox::Pointer< IBAMR::INSStaggeredHierarchyIntegrator     >  d_ins_hier_integrator;
   
  /*!
   * Pointers to the patch hierarchy object associated with the time integration object.
   */
   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_patch_hierarchy;
   int d_finest_ln , d_coarsest_ln;
   
  /*!
   * Pointer to hierarchy math operator.
   */ 
   SAMRAI::tbox::Pointer< IBTK::HierarchyMathOps > d_hier_math_ops;

  /*!
   * The LDataManager is used to coordinate the distribution of Lagrangian
   * data on the patch hierarchy.
   */
  IBTK::LDataManager* d_lag_data_manager;
  
  /*!
   * Name of method and name of LNodeLevelData object which stores information
   * of how lagrangian points are updated in space.
   */
  std::string d_lag_position_update_method;
  
  /*!
   * Pointer to kinematics object. The kinematics object is used to calculate deformational
   * kinematics or imposed kinematics at different lagrangian points.
   */
   SAMRAI::tbox::Pointer<IBAMR::IBKinematics> d_kinematics_ptr;
   
  /*!
   * Total no. of Lagrangian nodes.
   * Total No. of Lagrangian points at global and local level which are to follow rigidity 
   * constraint.
   */
  int d_lagrangian_global_nodes, d_rigid_global_nodes, d_rigid_local_nodes;
 
  /*!
   * Incremented angle from x, y and z axis.
   */
  SAMRAI::tbox::Array<double> d_incremented_angle_from_reference_axis;
  
   
  /*!
   *Bools to control the variants of rigidity constraint method for class of different applications.
   */
   bool d_needs_div_free_projection, d_kinematics_need_filtering, d_body_is_partly_elastic, 
        d_body_is_self_rotating, d_body_is_self_translating;
   SAMRAI::tbox::Array<int> d_calculate_translational_momentum, d_calculate_rotational_momentum;
   
  /*!
   * Vector of translational & rotational velocity components of the body.
   * Vectors of translational and rotational velocity introduced due to deformational kinematics.
   */
  std::vector<double>  d_rigid_trans_vel, d_rigid_trans_vel_old,
                       d_rigid_rot_vel, d_rigid_rot_vel_old,
                       d_vel_com_def, d_vel_com_def_old, 
		       d_omega_com_def, d_omega_com_def_old;
  
  /*!
   * Vector of center of mass of the non elastic body.
   */
  std::vector<double> d_center_of_mass_non_elastic_body;
  
  /*!
   * Vector of coordinates of the tagged point on the body of the structure.
   */
  int d_tagged_pt_idx;
  std::vector<double> d_tagged_point_position;
   
  /*!
   * Fluid and material properties. 
   */
   double d_rho_fluid, d_mu_fluid;

 // Printing out the information as the code proceeds   
   
  /*!
   * Iteration_counter for printing stuff.
   */
   int d_timestep_counter, d_output_interval, d_INS_num_cycles, d_INS_cycle_num;
   
  /*!
   * Bools for outputing stuff which is calculated on the fly.
   */
   bool d_do_log, d_output_drag_and_kinetic_energy, d_output_power, d_print_output,d_output_trans_vel, 
        d_output_rot_vel, d_output_COM_coordinates, d_output_MOI;
   
  /*!
   * output file name string.
   */
   std::string d_dir_name, d_base_output_filename;
   
  /*!
   * Pointer to Lagrangian velocity datas. The d_lag_vel_correction_data is workable data on which 
   * velocity correction is done. 
   * vel_correction = U_RBM - U_interpolated
   * vel_old        = U_RBM,n
   * vel_new        = 0.5( U_RBM,n+1 + U_RBM,n)
   * lag_new        = contains the information related to new position of the lagrangian markers.
   */
  SAMRAI::tbox::Pointer< IBTK::LNodeLevelData > d_lag_vel_correction_data, d_lag_old_vel_data, d_lag_new_vel_data, d_lag_new_position_data;
  
  /*!
   *  Variables and variable contexts associated with calculating divergence free projection.
   */
  SAMRAI::tbox::Pointer< SAMRAI::pdat::SideVariable< NDIM,double > > d_U_var;
  SAMRAI::tbox::Pointer< SAMRAI::pdat::CellVariable< NDIM,double > > d_U_cc_var;
  SAMRAI::tbox::Pointer< SAMRAI::pdat::SideVariable< NDIM,double > > d_U_fluidSolve_var;
  SAMRAI::tbox::Pointer< SAMRAI::pdat::CellVariable< NDIM,double > > d_Phi_var;
  SAMRAI::tbox::Pointer< SAMRAI::pdat::CellVariable< NDIM,double > > d_Div_U_var;
   
  SAMRAI::tbox::Pointer< SAMRAI::hier::VariableContext > d_scratch_context, d_new_context;
  int d_U_scratch_idx, d_U_fluidSolve_idx , d_Phi_idx, d_Div_U_scratch_idx, d_U_cc_scratch_idx;
  
  /*!
   * Hierarchy operations object. Needed for projection step.
   */
  SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > d_hier_sc_data_ops;
  SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
  SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >             d_wgt_cc_var;
  SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>                d_no_fill_op;
  int                                                                         d_wgt_cc_idx;
  double                                                                      d_volume;
  
  
  /*!
   * The following variables are needed to solve cell centered poison equation for \phi ,which is
   * used to project the corrected background fluid velocity on divergence free field to remove kinkiness
   * introduced via rigidbody velocity constraint.
   */
  SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>        d_velcorrection_projection_bc_coef;
  SAMRAI::solv::PoissonSpecifications*                 d_velcorrection_projection_spec;
  SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator>       d_velcorrection_projection_op;
  SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> d_velcorrection_projection_solver;
  SAMRAI::tbox::Pointer<IBTK::CCPoissonFACOperator>    d_velcorrection_projection_fac_op;
  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>        d_velcorrection_projection_fac_pc_db;
  SAMRAI::tbox::Pointer<IBTK::FACPreconditioner>       d_velcorrection_projection_fac_pc;


 
  /*!
   * Coarsening Algorithm object associated with fluidSolve_var object. Coarsening operation
   * synchronizes the data in all the invalid regions of the hierarchy.
   */
  typedef std::map<std::string,SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > >              CoarsenAlgMap;
  typedef std::map<std::string,SAMRAI::xfer::CoarsenPatchStrategy<NDIM>* >                                 CoarsenPatchStrategyMap;
  typedef std::map<std::string,std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > > CoarsenSchedMap;

  CoarsenAlgMap           d_calgs;
  CoarsenPatchStrategyMap d_cstrategies;
  CoarsenSchedMap         d_cscheds;
   
  /*!
   * Variables associated with Power Calculation.
   */
  SAMRAI::tbox::Pointer< SAMRAI::pdat::CellVariable< NDIM,double> > d_VisDefPower_var;
  SAMRAI::tbox::Pointer< SAMRAI::pdat::CellVariable< NDIM,double> > d_EulDefVel_var;
  SAMRAI::tbox::Pointer< SAMRAI::pdat::CellVariable< NDIM,double> > d_ConstraintPower_var;
  SAMRAI::tbox::Pointer< SAMRAI::pdat::CellVariable< NDIM,int   > > d_EulTagDefVel_var;
  int d_VisDefPower_scratch_idx, d_EulDefVel_scratch_idx, d_ConstraintPower_scratch_idx, d_EulTagDefVel_scratch_idx;
  
  /*!
   * Refine algorithm and schedule needed to fill the ghost box needed to calculate gradient of 
   * deformation velocities.
   */
  SAMRAI::xfer::RefineAlgorithm<NDIM> d_PowerRefineAlgorithm;

  /*!
   * File streams associated for the output.
   */
   std::fstream d_trans_vel_stream, d_rot_vel_stream, d_drag_force_stream, d_moment_of_inertia_stream, d_kinetic_energy_stream,
                d_position_COM_stream, d_power_spent_stream;
   		
  /*!
   * Volume of the cell at finest level, total volume of the body and volume of each material point.
   */
   double d_finest_cell_vol, d_volume_of_body, d_volume_of_non_elastic_lag_pt, d_volume_of_elastic_lag_pt;
		
  /*! 
   * Time after Fluid Solve. Its the time at which Rigidity Constraint has to be applied.
   */
  double d_integrator_time;
 
 
  /*!
   * Map containing lagrangian index number as key and vector of 
   * lagrangian node position as mapped_value. It contains local lag pts, but with global index as the key.
   */
  std::map< int, std::vector< double > >  d_map_lag_position; 

  /*!
   * Map containing Inertia component as key and Inertia component
   * value as mapped_value.
   */
  std::map< std::string, double > d_map_inertiaTensor;

  /*!
   * Blitz array of position vectors of lagrangian nodes.
   * Blitz array containing (u_interp - u_def)
   */
  blitz::Array<double,2> d_x_array, d_u_diff_array;
  
}; // class RigidityConstraint

}// namespace IBAMR

#endif //#define included_RigidityConstraint
