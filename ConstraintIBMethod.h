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
 
#ifndef included_constraintibmethod
#define included_constraintibmethod

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Pointer.h>
#include <tbox/Database.h>

// IBAMR INCLUDES
#include <ibamr/IBMethod.h>
#include <ibamr/INSHierarchyIntegrator.h>
#include "../../Examples/IBKinematics.h"

// C++ STDLIB INCLUDES
#include <vector>


namespace IBAMR
{

/*!
 * \brief Class ConstraintIBMethod implements the FuRMoRP algorithm
 * of the IB method.
 */

class ConstraintIBMethod : public IBAMR::IBMethod
{

  
public:
  
    /*!
     * \brief Constructor 
     */
    ConstraintIBMethod(
        const std::string& object_name,  
	SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
	bool register_for_restart = true,
	SAMRAI::tbox::Pointer< IBAMR::INSHierarchyIntegrator > ins_hier_integrator
        );
        
    /*!
     * \brief Destructor 
     */
    ~ConstraintIBMethod();
    
    /*!
     * \brief Register kinematics of the immersed structure(s) with this class.
     */
    void
    registerKinematics(
        std::vector< SAMRAI::tbox::Pointer<IBAMR::IBKinematics> > ib_kinematics_op);
    
    /*!
     * \brief Apply the FuRMoRP algorithm in IBAMR::IBStrategy::postprocessSolveFluidEquations method.
     */    
    virtual void
    postprocessSolveFluidEquations(
        double current_time, 
	double new_time, 
	int cycle_num);
    
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
        const ConstraintIBMethod& from );
    
    /*!
     * \brief Default assignment operator.
     *
     * \note This assignment operator is not implemented and should not be used.
     */
    ConstraintIBMethod&
    operator = (const ConstraintIBMethod& that);
    
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
     * \brief Set the INS cycle number for this object.
     */
    void
    setINSCycleNumber(
        const int cycle_num);
   
    /*!
     * \brief Set the time at which FuRMoRP is applied.
     */
    void
    setFurmorpTime(
        const double new_time);
    
    /*!
     * Calculate & store kinematics velocity on different lagrangian points.
     */
    void 
    calculateKinematicsVelocity(
        const double current_time, 
	const double new_time);

    /*!
     * Interpolate fluid solve velocity from Eulerian grid onto the Lagrangian mesh.
     */
    void
    interpolateFluidSolveVelocity();

    /*!
     * Subtract the calculated deformation velocity from the interpolated fluidSolve 
     * velocity.
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
    getNumberOfINSCycles() const
    {
        return d_INS_num_cycles;
    
    }//getNoOfINSCycles

    /*!
     * Solve for rigid rotational velocity.
     */
    void
    solveRigidRotationalVelocity();

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
     * if div free projection step is needed?
     */
    bool
    needs_DivFreeProjection() const 
    {
        return d_needs_div_free_projection;
    
    }//needs_DivFreeProjection
    

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
     * Translate and rotate the shape with rigid body velocity.
     */
    void
    updateLagrangianPosition(
        const double current_time, 
	const double new_time);

   
    /*!
     * Calculate Drag, Kinetic Energy, Power.
     */
    void 
    postprocessData(
        const double current_time, 
        const double new_time);
   
/////////////////////////    PRIVATE DATA MEMBERS ////////////////////////////////
    
    /*!
     * Cache IBAMR-IBTK pointers for functionality.
     */
    SAMRAI::tbox::Pointer< IBAMR::INSHierarchyIntegrator > d_ins_hier_integrator;
    SAMRAI::tbox::Pointer< IBTK::HierarchyMathOps > d_hier_math_ops;
    
    /*!
     * Pointer to the kinematics of the immersed structures.
     */
    
    
    /*!
     * Name of the immersed Lagrangian structures in the domain.
     */
    std::vector<std::string> d_structure_names; 
    
    /*!
     * ID of immersed Lagrangian structures in the domain.
     */
    std::vector<int> d_structure_ids;
    
    /*!
     * No of immersed structures.
     */
    const int d_no_structures;
    
    /*!
     * Bools to control the variants of FuRMoRP algorithm for class of different applications.
     */
    bool d_needs_div_free_projection, d_kinematics_need_filtering, d_body_is_self_rotating, d_body_is_self_translating;
    std::map<int, std::vector<bool> > d_calculate_translational_momentum, d_calculate_rotational_momentum;
        
    /*!
     * Rigid translational velocity of the structures.
     */
    std::map<int, std::vector<double> > d_rigid_trans_vel, d_rigid_trans_vel_old;
    
    /*!
     * Rigid rotational velocity of the structures.
     */
    std::map<int, std::vector<double> > d_rigid_rot_vel, d_rigid_rot_vel_old;
    
    /*!
     * Incremented angle from x, y and z axis when the body is rotating.
     */
    std::map<int, SAMRAI::tbox::Array<double> > d_incremented_angle_from_reference_axis;
    
    /*!
     * Translational velocity of the structures due to deformational kinematics.
     */
    std::map<int, std::vector<double> > d_vel_com_def, d_vel_com_def_old;
    
    /*!
     * Rotational velocity of the structures due to deformational kinematics.
     */
    std::map<int, std::vector<double> > d_omega_com_def, d_omega_com_def_old;
    
    /*!
     * Center of mass of the immersed structures.
     */
    std::map<int, std::vector<double> > d_center_of_mass;
    
    /*!
     * Tag a Lagrangian point (generally eye of the fish) of the immersed structures.
     */
    std::map<int, std::vector<int> > d_tagged_pt_lag_idx;
    
    /*!
     * Coordinates of the tagged points of different structures.
     */
    std::map<int, std::vector<double> > d_tagged_pt_position;
    
    /*!
     * Density and viscosity of the fluid.
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
    bool d_output_drag_and_kinetic_energy, d_output_power, d_print_output,d_output_trans_vel, 
         d_output_rot_vel, d_output_COM_coordinates, d_output_MOI;
   
    /*!
     * output file name string.
     */
    std::string d_dir_name, d_base_output_filename;
     
};
   
}

























#endif //#ifndef included_constraintibmethod