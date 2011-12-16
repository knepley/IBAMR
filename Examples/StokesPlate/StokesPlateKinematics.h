// Filename: StokesPlateKinematics.h
// Written by amneet bhalla on meghnad@mech.northwestern.edu
// Created on 12/14/2011.

// This is a concrete class which provides kinematics of Stokes' plate to ConstraintIBMethod class.
     
 
#ifndef included_stokesplatekinematics
#define included_stokesplatekinematics

///////////////////////////////////////// INCLUDES //////////////////////////////////////////

//IBAMR INCLUDES
#include "../../ConstraintIBKinematics.h"


//C++ INCLUDES
#include <vector>

namespace IBAMR
{
  
  /*!
   * \brief Class StokesPlateKinematics provides definition for base ConstraintIBKinematics class.
   */
class StokesPlateKinematics : public ConstraintIBKinematics
{
  
  
public:
    
    /*!
     * \brief Constructor.
     */
    StokesPlateKinematics(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
	IBTK::LDataManager* l_data_manager,
        bool register_for_restart = true);
    
    /*!
     * Destructor.
     */
    virtual 
    ~StokesPlateKinematics();
    
    /*!
     * Set kinematics velocity at new time for Stokes' plate.
     */
    virtual void
    setNewKinematicsVelocity(
        const double Time,
        const std::vector<double>& incremented_angle_from_reference_axis,
        const std::vector<double>& center_of_mass,
        const std::vector<double>& tagged_pt_position);
    
    /*!
     * Get the kinematics velocity at new time for Stokes' plate on the specified level.
     * 
     */
    virtual const std::vector<std::vector<double> >&
    getNewKinematicsVelocity(const int level) const;
    
    /*!
     * Get the kinematics velocity at current time for Stokes' plate on the specified level.
     * 
     */
    virtual const std::vector<std::vector<double> >&
    getCurrentKinematicsVelocity(const int level) const;
    
  
    /*!
     * Set the shape of Stokes' plate at new time for the structure on all levels.
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
    StokesPlateKinematics();
    
    /*!
     * Deleted default copy ctor.
     */
    StokesPlateKinematics(
        const StokesPlateKinematics& from);
  
    /*!
     * Deleted default assignment.
     */
    StokesPlateKinematics&
    operator = (const StokesPlateKinematics& that);
  
    
    /*!
     * New kinematics velocity. New shape of the body.
     * 
     * \NOTE Current velocity is always equal to new velocity. Position is 
     * updated via CONSTRAINT_VELOCITY method, so new shape is not filled in.
     */
    std::vector<std::vector<double> >  d_new_kinematics_vel, d_new_shape;
  
  
  
}; //class StokesPlateKinematics 
  
 
} //namespace IBAMR

#endif //included_stokesplatekinematics