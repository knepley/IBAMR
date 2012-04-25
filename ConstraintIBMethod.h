// Filename: ConstraintIBMethod.h
// Written by amneet bhalla on meghnad@mech.northwestern.edu
// Created on 12/01/2011.

// This is a concrete class which implements fast and efficient distributed Lagrange multiplier (fictitious domain) method
// for immersed structures. The particular algorithm implemented is FuRMoRP.

// References
//  *  Patankar et al. A new formulation of the distributed Lagrange multiplier/fictitious domain method
//     for particulate flows. Int. Journal of Multiphase flows, 26, 1509-1524 (2000).
//  *  Shirgaonkar et al. A new mathematical formulation and fast algorithm for fully resolved simulation of self-propulsion.
//     JCP, 228 , 2366-2390 (2009).

#ifndef included_ConstraintIBMethod
#define included_ConstraintIBMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSC INCLUDES
#include <petscsys.h>

// SAMRAI INCLUDES
#include <tbox/Pointer.h>
#include <tbox/Database.h>
#include <Variable.h>
#include <CellVariable.h>
#include <SideVariable.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <VariableContext.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>

// IBAMR INCLUDES
#include <ibamr/IBMethod.h>
#include <ibamr/INSHierarchyIntegrator.h>
#include <ibamr/IBHierarchyIntegrator.h>
#include "ConstraintIBKinematics.h"

// IBTK INCLUDES
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/CCPoissonFACOperator.h>
#include <ibtk/FACPreconditioner.h>

// BLITZ INCLUDES
#include <blitz/array.h>

// C++ STDLIB INCLUDES
#include <vector>
#include <string>
#include <fstream>

namespace IBAMR
{
/*!
 * \brief Class ConstraintIBMethod implements the FuRMoRP algorithm
 * of the IB method.
 */
class ConstraintIBMethod
    : public IBAMR::IBMethod
{

public:

    /*!
     * \brief Constructor
     */
    ConstraintIBMethod(
        const std::string& object_name,
	SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
	const int no_structures,
	bool register_for_restart = true);

    /*!
     * \brief Destructor
     */
    ~ConstraintIBMethod();

    /*!
     *  Initialize Hierarchy operators
     */
   void
   initializeHierarchyOperatorsandData();

    /*!
     * Register Eulerian variables with base IBStrategy class.
     */
    virtual void
    registerEulerianVariables();

    /*!
     *  Register Eulerian communication algorithms with base IBStrategy class.
     */
    virtual void
    registerEulerianCommunicationAlgorithms();
    
    /*!
     *  Create Lagrangian workspace.
     */
    virtual void
    preprocessIntegrateData(
        double current_time,
	double new_time,
	int num_cycles);
    
    /*!
     *  Destroy Lagrangian workspace.
     */
    virtual void
    postprocessIntegrateData(
        double current_time,
	double new_time,
	int num_cycles);

    /*!
     * \brief Register kinematics of the immersed structure(s) with this class.
     */
    void
    registerConstraintIBKinematics(
        const std::vector<SAMRAI::tbox::Pointer<IBAMR::ConstraintIBKinematics> >& ib_kinematics_op);

    /*!
     * \brief Register any preprocess fluid solve callback functions.
     */
    inline void
    registerPreProcessSolveFluidEquationsCallBackFunction(
	void (*ptr_preprocess_callbackfnc)(const double, const double, const int, void* ), 
        void* ctx)
    {

        d_prefluidsolve_callback_fns.push_back(ptr_preprocess_callbackfnc);
        d_prefluidsolve_callback_fns_ctx.push_back(ctx);

        return;
    }//registerPreProcessSolveFluidEquationsCallBackFunction

    /*!
     * \brief Calculate any body forces for INS solver over here.
     */
    virtual void
    preprocessSolveFluidEquations(
        double current_time,
	double new_time,
	int cycle_num);

    /*!
     * \brief Register any postprocess fluid solve callback functions.
     */
    inline void
    registerPostProcessSolveFluidEquationsCallBackFunction(
	void (*ptr_postprocess_callbackfnc)(const double, const double, const int, void* ) , 
        void* ctx)
    {

        d_postfluidsolve_callback_fns.push_back(ptr_postprocess_callbackfnc);
        d_postfluidsolve_callback_fns_ctx.push_back(ctx);

        return;
    }//registerPostProcessSolveFluidEquationsCallBackFunction

    /*!
     * \brief Apply the FuRMoRP algorithm in IBAMR::IBStrategy::postprocessSolveFluidEquations method.
     */
    virtual void
    postprocessSolveFluidEquations(
        double current_time,
	double new_time,
	int cycle_num);

    /*!
     * \brief Override the eulerStep method of the base IBMethod class.
     */
    virtual void
    eulerStep(
        double current_time,
	double new_time);

    /*!
     * \brief Override the midpointStep method of the base IBMethod class.
     */
    virtual void
    midpointStep(
        double current_time,
	double new_time);

    /*!
     * \brief Override the putToDatabase method of the base Serializable class.
     */
    virtual void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Get the volume element associated with material points.
     */
    inline const std::vector<double>&
    getVolumeElement()
    {
      return d_vol_element;
      
    }//getVolumeElement
    
    /*!
     * Get LData associated with Lagrange multiplier force field.
     */
    inline const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >&
    getLagrangeMultiplier()
    {
        return d_l_data_U_correction; 
    }
    
private:
 
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    ConstraintIBMethod();

    /*!
     * \brief Default copy constructor.
     *
     * \note This copy constructor is not implemented and should not be used.
     */
    ConstraintIBMethod(
        const ConstraintIBMethod& from);

    /*!
     * \brief Default assignment operator.
     *
     * \note This assignment operator is not implemented and should not be used.
     */
    ConstraintIBMethod&
    operator=(
        const ConstraintIBMethod& that);

    /*!
     * \brief Get values from input file.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
	const bool from_restart);

    /*!
     * \brief Get values from restart file.
     */
    void
    getFromRestart();

    /*!
     * Set initial Lagrangian velocity on material points.
     */
    void
    setInitialLagrangianVelocity();

    /*!
     * Calculate center of mass and moment of inertia of immersed
     * structures.
     */
    void
    calculateCOMandMOIOfStructures();

    /*!
     * Calculate & store kinematics velocity on different lagrangian points.
     */
    void
    calculateNewKinematicsVelocity();

    /*!
     * Calculate momentum of Kinematics velocity.
     */
    void
    calculateMomentumOfNewKinematicsVelocity(const int position_handle);

    /*!
     * Calculate volume element associated with material points.
     */
    void
    calculateVolumeElement();

    /*!
     * \brief Set the INS cycle number for this object.
     */
    inline void
    setINSCycleNumberAndCounter(
        const int cycle_num)
    {
        d_INS_current_cycle_num = cycle_num;
        if (cycle_num == (d_INS_num_cycles -1))
        {
            ++d_timestep_counter;
        }
        return;
    }// setINSCycleNumberAndCounter

    /*!
     * \brief Set the time at which FuRMoRP is applied.
     */
    void
    setFuRMoRPTime(
        const double current_time,
        const double new_time)
    {
        d_FuRMoRP_current_time = current_time;
        d_FuRMoRP_new_time     = new_time;
        return;
    }// setFuRMoRPTime

    /*!
     * Interpolate fluid solve velocity from Eulerian grid onto the Lagrangian mesh.
     */
    void
    interpolateFluidSolveVelocity();

    /*!
     * Calculate the rigid translational velocity.
     */
    void
    calculateRigidTranslationalMomentum();
    
    /*!
     * Calculate the rigid rotational velocity.
     */
    void
    calculateRigidRotationalMomentum();

    /*!
     * Calculate current velocity on the material points.
     */
    void
    calculateCurrentLagrangianVelocity();

    /*!
     * Correct velocity on Lagrangian mesh. Set the velocity on Lagrangian
     * mesh as U_lag_corr = U_trans + Omega X r + U_def - U_interpolated.
     */
    void
    correctVelocityOnLagrangianMesh();

    /*!
     *  Spread the corrected velocity at the Lagrangian mesh to the Eulerian Grid.
     */
    void
    spreadCorrectedLagrangianVelocity();

    /*!
     * Synch the patch hierarchy.
     */
    void
    synchronizeLevels();

    /*!
     * The correction on Eulerian grid can lead to kinkiness in the velocity field. To remove such kinkiness
     * we project the corrected velocity field onto a divergence free field.
     */
    void
    applyProjection();

    /*!
     * Update the position of structures according to forward Euler step method.
     */
    void
    updateStructurePositionEulerStep();

    /*!
     * Update the position of structures according to mid point step method.
     */
    void
    updateStructurePositionMidPointStep();

    /*!
     * Compute U_half = 0.5(U_current + U_new);
     */
    void
    calculateMidPointVelocity();

    /*!
     * Calculate Drag on the immersed structures.
     */
    void
    calculateDrag();
    
    /*!
     * Calculate kinetic energy of the immersed structures.
     */
    void
    calculateKE();
 
    /*!
     * Calculate power spent during swimming.
     */
    void
    calculatePower();
       
    /*!
     * Hierarchy math operator.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;

    /*!
     * Type of fluid solver.
     */
    bool d_collocated_solver, d_staggered_solver;

    /*!
     * No of immersed structures.
     */
    const int d_no_structures;

    /*!
     * Pointer to the kinematics of the immersed structures.
     */
    std::vector<SAMRAI::tbox::Pointer<IBAMR::ConstraintIBKinematics> > d_ib_kinematics;

    /*!
     * FuRMoRP apply time.
     */
    double d_FuRMoRP_current_time, d_FuRMoRP_new_time;

    /*!
     * Volume element associated with material points.
     */
    std::vector<double> d_vol_element;

    /*!
     * If divergence free projection is needed after FuRMoRP algorithm?
     */
    bool d_needs_div_free_projection;

    /*!
     * Rigid translational velocity of the structures.
     */
    std::vector<std::vector<double> > d_rigid_trans_vel_current, d_rigid_trans_vel_new;

    /*!
     * Rigid rotational velocity of the structures.
     */
    std::vector<std::vector<double> > d_rigid_rot_vel_current, d_rigid_rot_vel_new;

    /*!
     * Incremented angle from x, y and z axis when the body is rotating.
     */
    std::vector<std::vector<double> > d_incremented_angle_from_reference_axis;

    /*!
     * Translational velocity of the structures due to deformational kinematics.
     */
    std::vector<std::vector<double> > d_vel_com_def_current, d_vel_com_def_new;

    /*!
     * Rotational velocity of the structures due to deformational kinematics.
     */
    std::vector<std::vector<double> > d_omega_com_def_current, d_omega_com_def_new;

    /*!
     * Center of mass of the immersed structures.
     */
    std::vector<std::vector<double> > d_center_of_mass_current, d_center_of_mass_new;

    /*!
     * Moment of inertia of the structures.
     */
    std::vector<blitz::Array<double,2> > d_moment_of_inertia_current, d_moment_of_inertia_new;

    /*!
     * Tag a Lagrangian point (generally eye of the fish) of the immersed structures.
     */
    std::vector<int> d_tagged_pt_lag_idx;

    /*!
     * Coordinates of the tagged points of different structures.
     */
    std::vector<std::vector<double> > d_tagged_pt_position;

    /*!
     * Density and viscosity of the fluid.
     */
    double d_rho_fluid, d_mu_fluid;

    /*!
     * Iteration_counter for printing stuff.
     */
    int d_timestep_counter, d_output_interval, d_INS_num_cycles, d_INS_current_cycle_num;

    /*!
     * Bools for outputing stuff which is calculated on the fly.
     */
    bool d_print_output, d_output_drag, d_output_kinetic_energy, d_output_power,d_output_trans_vel,
         d_output_rot_vel, d_output_COM_coordinates, d_output_MOI;

    /*!
     * output file name string.
     */
    std::string d_dir_name, d_base_output_filename;

    /*!
     * Store LData for only those levels which contain immersed structures.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_l_data_U_interp, d_l_data_U_correction,
        d_l_data_U_new, d_l_data_U_current, d_l_data_U_half, d_l_data_X_new_Euler, d_l_data_X_new_MidPoint;

    /*!
     * Hierarchy operations object. Needed for projection step.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>                d_no_fill_op;
    int                                                                         d_wgt_cc_idx;
    double                                                                      d_volume;

    /*!
     *  Variables and variable contexts associated with calculating divergence free projection.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_u_var;
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_u_fluidSolve_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_phi_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_Div_u_var;

    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_scratch_context;
    int d_u_scratch_idx, d_u_fluidSolve_idx , d_phi_idx, d_Div_u_scratch_idx;

    /*!
     * The following variables are needed to solve cell centered poison equation for \f$ \phi \f$ ,which is
     * used to project the corrected background fluid velocity on divergence free field to remove kinkiness
     * introduced via FuRMoRP algorithm.
     */
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>        d_velcorrection_projection_bc_coef;
    SAMRAI::solv::PoissonSpecifications*                 d_velcorrection_projection_spec;
    SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator>       d_velcorrection_projection_op;
    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> d_velcorrection_projection_solver;
    SAMRAI::tbox::Pointer<IBTK::CCPoissonFACOperator>    d_velcorrection_projection_fac_op;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>        d_velcorrection_projection_fac_pc_db;
    SAMRAI::tbox::Pointer<IBTK::FACPreconditioner>       d_velcorrection_projection_fac_pc;

    /*!
     * Variables associated with Power calculation.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_VisDefPower_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_EulDefVel_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_ConstraintPower_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,int> >    d_EulTagDefVel_var;
    int d_VisDefPower_scratch_idx, d_EulDefVel_scratch_idx, d_ConstraintPower_scratch_idx, d_EulTagDefVel_scratch_idx;

    /*!
     * File streams associated for the output.
     */
    std::vector<std::ofstream*>  d_trans_vel_stream, d_rot_vel_stream, d_drag_force_stream, d_moment_of_inertia_stream,
        d_kinetic_energy_stream,d_position_COM_stream, d_power_spent_stream;

    /*!
     * Pre and post fluid solve call back functions and contexts.
     */
    std::vector<void(*)(const double, const double, const int, void*)> d_prefluidsolve_callback_fns, d_postfluidsolve_callback_fns;
    std::vector<void*> d_prefluidsolve_callback_fns_ctx, d_postfluidsolve_callback_fns_ctx;
};
}// namespace IBAMR

///////////////////////////////////////// INLINE ////////////////////////

//#include "ConstraintIBMethod.I"

/////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ConstraintIBMethod
