// Filename: MovingPlateKinematics.h
// Written by amneet bhalla on meghnad@mech.northwestern.edu
// Created on 12/20/2011.

// This is a concrete class which provides kinematics of Moving plate to ConstraintIBMethod class.
     
 
#ifndef included_MovingPlateKinematics
#define included_MovingPlateKinematics

///////////////////////////////////////// INCLUDES //////////////////////////////////////////

//IBAMR INCLUDES
#include "../../ConstraintIBKinematics.h"


//C++ INCLUDES
#include <vector>

namespace IBAMR
{
  
  /*!
   * \brief Class MovingPlateKinematics provides definition for base ConstraintIBKinematics class.
   */
class MovingPlateKinematics 
    : public ConstraintIBKinematics
{
  
  
public:
    
    /*!
     * \brief Constructor.
     */
    MovingPlateKinematics(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
	IBTK::LDataManager* l_data_manager,
        bool register_for_restart = true);
    
    /*!
     * Destructor.
     */
    virtual 
    ~MovingPlateKinematics();
    
    /*!
     * Set kinematics velocity at new time for Moving plate.
     */
    virtual void
    setNewKinematicsVelocity(
        const double Time,
        const std::vector<double>& incremented_angle_from_reference_axis,
        const std::vector<double>& center_of_mass,
        const std::vector<double>& tagged_pt_position);
    
    /*!
     * Get the kinematics velocity at new time for Moving plate on the specified level.
     * 
     */
    virtual const std::vector<std::vector<double> >&
    getNewKinematicsVelocity(const int level) const;
    
    /*!
     * Get the kinematics velocity at current time for Moving plate on the specified level.
     * 
     */
    virtual const std::vector<std::vector<double> >&
    getCurrentKinematicsVelocity(const int level) const;
    
  
    /*!
     * Set the shape of Moving plate at new time for the structure on all levels.
     */
    virtual void
    setNewShape();
    
    /*!
     * Get the shape of structure at new time  on the specified level.
     */
    virtual const std::vector<std::vector<double> >&
    getNewShape(const int level) const;
    

private:
  
    /*!
     * Deleted default ctor.
     */
    MovingPlateKinematics();
    
    /*!
     * Deleted default copy ctor.
     */
    MovingPlateKinematics(
        const MovingPlateKinematics& from);
  
    /*!
     * Deleted default assignment.
     */
    MovingPlateKinematics&
    operator = (
        const MovingPlateKinematics& that);
  
    
    /*!
     * New kinematics velocity. New shape of the body.
     * 
     * \NOTE Current velocity is always equal to new velocity. Position is 
     * updated via CONSTRAINT_VELOCITY method, so new shape is not filled in.
     */
    std::vector<std::vector<std::vector<double> > >  d_new_kinematics_vel; 
    std::vector<std::vector<double> > d_new_shape;
  
  
  
}; // MovingPlateKinematics 
  
 
} //IBAMR

#endif //included_MovingPlateKinematics
