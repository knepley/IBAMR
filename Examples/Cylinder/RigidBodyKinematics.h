// Filename: RigidBodyKinematics.h
// Written by amneet bhalla on meghnad@mech.northwestern.edu
// Created on 1/10/2012.

// This is a concrete class which provides kinematics of a rigid solid to ConstraintIBMethod class.
     
 
#ifndef included_RigidBodyKinematics
#define included_RigidBodyKinematics

///////////////////////////////////////// INCLUDES //////////////////////////////////////////

//IBAMR INCLUDES
#include "../../ConstraintIBKinematics.h"


//C++ INCLUDES
#include <vector>


/////////////////////////////////////// FORWARD DECLARATION ////////////////////////////////

namespace mu
{
  class Parser;
}// namespace mu

namespace IBAMR
{
  
  /*!
   * \brief Class RigidBodyKinematics provides definition for base ConstraintIBKinematics class.
   */
class RigidBodyKinematics 
    : public ConstraintIBKinematics
{
  
  
public:
    
    /*!
     * \brief Constructor.
     */
    RigidBodyKinematics(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
	IBTK::LDataManager* l_data_manager,
        bool register_for_restart = true);
    
    /*!
     * Destructor.
     */
    virtual 
    ~RigidBodyKinematics();
    
    /*!
     * Set kinematics velocity at new time for rigid body.
     */
    virtual void
    setNewKinematicsVelocity(
        const double new_time,
        const std::vector<double>& incremented_angle_from_reference_axis,
        const std::vector<double>& center_of_mass,
        const std::vector<double>& tagged_pt_position);
    
    /*!
     * Get the kinematics velocity at new time for rigid body on the specified level.
     */
    virtual const std::vector<std::vector<double> >&
    getNewKinematicsVelocity(const int level) const;
    
    /*!
     * Get the kinematics velocity at current time for rigid body on the specified level.
     */
    virtual const std::vector<std::vector<double> >&
    getCurrentKinematicsVelocity(const int level) const;
    
  
    /*!
     * Set the shape of rigid body at new time on all levels.
     */
    virtual void
    setNewShape();
    
    /*!
     * Get the shape of rigid body at new time on the specified level.
     */
    virtual const std::vector<std::vector<double> >&
    getNewShape(const int level) const;

    /*!
     * Override the base Serializable method.
     */
    virtual void
    putToDatabase(
	SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);
    

private:
  
    /*!
     * Deleted default ctor.
     */
    RigidBodyKinematics();
    
    /*!
     * Deleted default copy ctor.
     */
    RigidBodyKinematics(
        const RigidBodyKinematics& from);
  
    /*!
     * Deleted default assignment.
     */
    RigidBodyKinematics&
    operator = (
        const RigidBodyKinematics& that);

    /*!
     * Get necessary data from restart manager for restarted runs.
     */
    void
    getFromRestart();

    /*!
     * Set rigid body velocity.
     */
    void
    setRigidBodySpecificVelocity(
        const double time,
        const std::vector<double>& incremented_angle_from_reference_axis,
        const std::vector<double>& center_of_mass,
        const std::vector<double>& tagged_pt_position);

    /*!
     * The mu::Parser objects which evaluate the data-setting functions.
     */
    std::vector<mu::Parser*> d_kinematicsvel_parsers;
    std::vector<mu::Parser*> d_all_parsers;
  
    /*!
     * Input kiematics velocity functions.
     */
    std::vector<std::string> d_kinematicsvel_function_strings;

    /*!
     * Parser variables.
     */
    double* d_parser_time;    
    double* d_parser_posn;

    /*!
     * Current time (t) and new time (t+dt).
     */
    double d_current_time, d_new_time;

    /*!
     * New kinematics velocity. New shape of the body.
     * 
     * \NOTE Current velocity is always equal to new velocity. Position is 
     * updated via CONSTRAINT_VELOCITY method, so new shape is not filled in.
     */
    std::vector<std::vector<std::vector<double> > >  d_new_kinematics_vel, d_current_kinematics_vel; 
    std::vector<std::vector<double> > d_new_shape;

    /*!
     * Save COM, tagged point position and incremented angle from reference axis for restarted runs.
     */
    std::vector<double> d_center_of_mass, d_incremented_angle_from_reference_axis, d_tagged_pt_position;
 
}; //RigidBodyKinematics 
  
 
} //IBAMR

#endif //included_RgidBodyKinematics
