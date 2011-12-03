/******************************************************************************************************************************
   Filename : IBKinematics.h
   Written by amneet bhalla on meghnad@mech.northwestern.edu
   Created on 6/14/2011.
 
*******************************************************************************************************************************/
 
/////////////////////////////////////////////////// INCLUDE GAURDS /////////////////////////////////////////////////////

#ifndef included_IBKinematics
#define included_IBKinematics

// SAMRAI INCLUDES
#include <tbox/DescribedClass.h>
#include <tbox/Array.h>

// C++ INCLUDES
 #include <map>
 #include <vector>



 namespace IBAMR
 {

   /*!
    * \brief Abstract class for specifing kinematics (deformational or otherwise).
    * 
    *  This class provides an abstraction to the IBAMR::RigidityConstraint class to get kinematics and the
    *  updated shape of the body.
    */

 class IBKinematics : public virtual SAMRAI::tbox::DescribedClass
 {
   
   ////////////////////////////////// PRIVATE MEMBERS//////////////////////////////////////////////////////////////////

   // Default copy constructor is not implemented and should not be used.
   IBKinematics(const IBKinematics& from);

   // Default assignment operator is not implemented and should not be used.
   IBKinematics& 
   operator = (const IBKinematics& that);
   

 
   ////////////////////////////////// PUBLIC MEMBERS//////////////////////////////////////////////////////////////////
   public:
     
   /*!
    * Data structures needed to store global lagrangian kinematic velocity and shape information.
    */  
   std::map< int,std::vector<double> > d_map_kinematics_vel, d_map_filtered_kinematics_vel, d_map_new_shape;
  
   /*!
    * Default constructor. It does nothing interesting.
    */
   IBKinematics();

   /*!
    * Virtual destructor.
    */
   virtual
     ~IBKinematics();
   

   /*!
    * \brief The following virtual function is used to calculate kinematic velocity of the material points. 
    * 
    * The kinematic velocity is calculated at t_n+1 = t_n + dt. It is calculated at cycle_num = 0 and 
    * remains the same for the consecutive cycles.
    * 
    * \param Time is the time at which kinematic velocity is calculated. It is at t_n+1 = t_n + dt.
    * 
    * \param incremented_angle_from_reference_axis is the angle made by body axis with the reference axis.
    * 
    * \param X_COM is vector containing coordinates of the center of mass.
    * 
    * \param X_tagged is vector containing coordinates of the tagged point .
    */
   virtual void  
   calculateGivenKinematicsVelocity(const double Time, const SAMRAI::tbox::Array<double>& incremented_angle_from_reference_axis, 
				   const std::vector<double>& X_COM, const std::vector<double>& X_tagged ) = 0;

   /*!
    * \brief The following virtual function is used to update the shape of immersed body. 
    * 
    * The body shape and its orientation is calculated in this method and IBAMR::RigidityConstraint translates 
    * the body to correct position for t_n+1/2 = 0.5*( t_n+1 + t_n).
    * 
    * \param incremented_angle_from_reference_axis is the incremented angle made by body with the reference axis.
    * \param new_time is the time at which body shape has to be updated. It is at t_n+1 = t_n + dt.
    */
   virtual void 
   updateLagrangianMarkerShapeAndOrientation(const SAMRAI::tbox::Array<double>& incremented_angle_from_reference_axis,
					     const double new_time) = 0 ;
  

 };// IBKinematics


 } //namespace IBAMR
 


#endif // #ifndef included_IBKinematics 
