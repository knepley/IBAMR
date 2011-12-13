// Filename: ConstraintIBKinematics.h
// Written by amneet bhalla on meghnad@mech.northwestern.edu
// Created on 12/10/2011.

// This is an abstract class which encapsulates structure information and provides kinematics and updated shape of the structure to 
// ConstraintIBMethod class.
     
 
#ifndef included_constraintibkinematics
#define included_constraintibkinematics

///////////////////////////////////////// INCLUDES //////////////////////////////////////////


//SAMRAI INCLUDES
#include <tbox/DescribedClass.h>
#include <tbox/Serializable.h>
#include <tbox/Database.h>

//IBAMR INCLUDES


//IBTK INCLUDES
#include <ibtk/LDataManager.h>

//C++ INCLUDES
#include <vector>

namespace IBAMR
{
  
/*!
 * \brief Class ConstraintIBKinematics encapsulates structure information and provides abstraction to get 
 * kinematics(deformational or imposed) of immersed structure to ConstraintIBMethod class.
 */

class ConstraintIBKinematics
    : public virtual SAMRAI::tbox::DescribedClass,
      public SAMRAI::tbox::Serializable
{

  
public:
  
class StructureParameters
{
    
public:  
  
    /*!
     * \brief Constructor.
     */
    StructureParameters(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
	IBTK::LDataManager* l_data_manager);
    
    /*!
     * Lagrangian point to tag on this structure.
     */
    int 
    getTaggedPtIdx() const;

    /*!
     * Check if the structure has translational degree unlocked.
     */
    bool
    getStructureIsSelfTranslating() const;
    
    /*!
     * Check if the structure has rotational degree unlocked.
     */
    bool
    getStructureIsSelfRotating() const;

    /*!
     * The coarsest level on which the structure resides.
     */
    int
    getCoarsestLevelNumber() const;
    
    /*!
     * The finest level on which the structure resides.
     */
    int
    getFinestLevelNumber() const;
    
    /*!
     * Global Lagrangian indices managed for this structure.
     */
    const std::vector<std::pair<int,int> >& 
    getLagIdxRange() const;

    /*!
     * Total number of Lagrangian nodes managed for this structure.
     */
    int 
    getTotalNodes() const;

    /*!
     * Lagrangian nodes update method for this structure.
     */
    std::string
    getPositionUpdateMethod() const;
   
private:
  
     std::string d_lag_position_update_method; 
     int d_coarsest_ln, d_finest_ln;
     std::vector<std::pair<int,int> > d_idx_range;
     int d_total_nodes;
     int d_tagged_pt_idx;
     bool d_struct_is_self_translating, d_struct_is_self_rotating;
    
};
    /*!
     * \brief Constructor.
     */
    ConstraintIBKinematics(
        std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
	IBTK::LDataManager* l_data_manager,
        bool register_for_restart = true);
    
    /*!
     * \brief Destructor.
     */
    virtual
    ~ConstraintIBKinematics();
    
    /*!
     * Get the object enclosing this structure's parameters.
     */
    const StructureParameters&
    getStructureParameters() const;

    /*!
     * Set the kinematics velocity(deformational or imposed) at Lagrangian points managed by this object.
     * 
     * \param Time new time at which kinematics velocity is to be set.
     * 
     * \param incremented_angle_from_reference_axis angle made with x,y & z axis due to rigid rotational velocity.
     * \f$ \theta_n = \theta_{n-1} + \omega_{n-1} \triangle t \f$
     * 
     * \param center_of_mass COM of the structure at current time.
     * 
     * \param tagged_pt_position Coordinates of the tagged point on this structure.
     */
    virtual void
    setNewKinematicsVelocity(
        const double Time,
        const std::vector<double>& incremented_angle_from_reference_axis,
        const std::vector<double>& center_of_mass,
        const std::vector<double>& tagged_pt_position) = 0;

    /*!
     * Get the kinematics velocity at new time for the structure on the specified level.
     * 
     * \param level kinematics velocity of the structure on this level
     */
    virtual const std::vector<std::vector<double> >&
    getNewKinematicsVelocity(const int level) const = 0;
    
    
    /*!
     * Get the kinematics velocity at current time for the structure on the specified level.
     * 
     * \param level kinematics velocity of the structure on this level
     */
    virtual const std::vector<std::vector<double> >&
    getCurrentKinematicsVelocity(const int level) const = 0;
    
    /*!
     * Set the shape of structure at new time for the structure on all levels.
     */
    virtual void
    setNewShape() = 0;

    /*!
     * Get the shape of structure at new time  on the specified level.
     * 
     * \param level new shape of the structure on this level
     */
    virtual const std::vector<std::vector<double> >&
    getNewShape(const int level) const = 0;
    
    /*!
     * Write out object state to the given database.
     *
     * An empty default implementation is provided.
     */
    virtual void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);
    
///////////////////////////////////////PRIVATE////////////////////////////////// 
private:    
  
    /*!
     * Deleted default ctor.
     */
    ConstraintIBKinematics();
    
    /*!
     * Deleted default copy ctor.
     */
    ConstraintIBKinematics(
        const ConstraintIBKinematics& from);
   
    /*!
     * Object enclosing all the parameters of the structure.
     */ 
    StructureParameters d_struct_param;
    
///////////////////////////////////////PROTECTED////////////////////////////////// 
protected:
  
    /*!
     * Name of the object.
     */
    std::string d_object_name;
    
    /*!
     * If the object is registred for restart.
     */
    bool d_registered_for_restart;
    
       
};//ConstraintIBKinematics


} //IBAMR


/////////////////////////////////// INLINE ///////////////////////////////////////

#include "ConstraintIBKinematics.I"

/////////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_constraintibkinematics