/******************************************************************************************************************************
   Filename : IBRigidStructureKinematics.h
   Written by amneet bhalla on Taylor@Courant Institute of Mathematical Sciences
   Created on 8/20/2011.


   This is a concrete class which sets the deformation velocity to zero and implements the shape update method of the abstract class 
   IBAMR::IBKinematics for a rigid structure needed by IBAMR::RigidityConstraint class.
 
*******************************************************************************************************************************/

#ifndef included_ibrigidstructurekinematics
#define included_ibrigidstructurekinematics


/////////////////////////////////////////// INCLUDES ///////////////////////////////////////

// APPLICATION INCLUDES
#include "../IBKinematics.h"

// SAMRAI INCLUDES
#include <tbox/Pointer.h>
#include <tbox/Database.h>
#include <PatchHierarchy.h>

// C++ INCLUDES
#include <string>

namespace IBAMR
{
  /*!
   * \brief  This is a concrete class which sets the deformation velocity to zero and implements the shape update method of the abstract class 
   * IBAMR::IBKinematics for a rigid structure needed by IBAMR::RigidityConstraint class.
   */

  class IBRigidStructureKinematics : public IBKinematics
  {

    /////////////////////// PRIVATE DATAMEMBERS //////////////////////////////
    /*!
     * Name of the object.
     */
    std::string d_object_name;



    ////////////////////////// PRIVATE MEMBER FUNCTIONS ////////////////////////

    /*!
     * \note The default ctor is not implemented and should not be used.
     */
    IBRigidStructureKinematics();

    /*!
     * \note The default copy ctor is not implemented and should not be used.
     */
    IBRigidStructureKinematics(const IBRigidStructureKinematics& from );

    /*!
     * \note The default assignment operator is not implemented and should
     * not be used.
     */
    IBRigidStructureKinematics&
    operator = (const IBRigidStructureKinematics& that );
  


    ////////////////////////// PUBLIC MEMBER FUNCTIONS /////////////////////////

 public: 
      
      /*!
       * ctor. This is the only ctor for this object.
       */
      IBRigidStructureKinematics( SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, 
	                SAMRAI::tbox::Pointer< SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);
      
      /*!
       * Destructor.
       */
      ~IBRigidStructureKinematics();
      
      /*!
       * calculate the deformation velocity at material point.
       */
       virtual void  
       calculateGivenKinematicsVelocity(const double Time, const SAMRAI::tbox::Array<double>& incremented_angle_from_reference_axis, 
				        const std::vector<double>& X_COM, const std::vector<double>& X_tagged );
           
      /*!
       * update the shape and orientation of the body for the next time step.
       */
       virtual void 
       updateLagrangianMarkerShapeAndOrientation(const SAMRAI::tbox::Array<double>& incremented_angle_from_reference_axis,
					     const double new_time);


  }; // class IBRigidStructureKinematics




}//namespace IBAMR



#endif // ifndef included_ibrigidstructurekinematics
