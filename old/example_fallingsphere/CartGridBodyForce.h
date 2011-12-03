
/****************************************************************************************************************

  // Filename: CartGridBodyForce.h
  // Created by amneet bhalla on 08/11/2011 @ Courant Institute Of Mathematical Sciences
 
****************************************************************************************************************/

// INCLUDE GAURDS
#ifndef included_CartGridBodyForce
#define included_CartGridBodyForce


//////////////////////////////////INCLUDES///////////////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Pointer.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <Variable.h>



// IBAMR INCLUDES



// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>



// C++ INCLUDES



namespace IBTK
{
    
  
  
  /*!
   * \brief CartGridBodyForce class.
   * 
   * This class is a concrete implementation of IBTK::CartGridFunction class and is used to set the
   * body force on patches for the fluid solver.
   */
    
    
    class CartGridBodyForce : public IBTK::CartGridFunction
    {

        ///////////////////////////// PRIVATE DATA MEMBERS /////////////////////
        /*!
	 * Patch index for Electric Field body force object.
	 */
        int d_sc_BodyForce_idx;
	
	////////////////////////////// PRIVATE MEMBER FUNCTIONS //////////////////
	/*!
	 * \note The default ctor is not implemented and should not be used.
	 */
	CartGridBodyForce();
	
	/*!
	 * \note The default copy ctor is not implemented and should not be used.
	 */
	CartGridBodyForce(const CartGridBodyForce& from);
	
        /*!
	 * \note The default assignment operator is not implemented and should not be used.
	 */
	CartGridBodyForce&
	operator = (const CartGridBodyForce& that);
	
       //////////////////////////// PUBLIC MEMBER FUNCTIONS ////////////////////////
	
       public:
	 
	/*!
	 * \brief Constructor
	 * 
	 * \param body_force_idx is the array index of the patch where body force due to Electric Field is
	 * calculated by IBAMR::IBElectricField class.
	 */
	CartGridBodyForce(const int body_force_idx);
	
	// The following methods redefine the pure virtual functions of the CartGridFunction abstract base class
	// isTimeDependent ()
	// setDataOnPatch()
	
	/*!
	 * Returns if the body force cartesian grid function is time dependent or not.
	 */
	virtual bool
	isTimeDependent() const;
	
        /*!
	 * This method copies the body force calculated apiori on the  patch interior passed to this method.
	 */
	virtual void
	setDataOnPatch(	const int data_idx,
                        SAMRAI::tbox::Pointer< SAMRAI::hier::Variable< NDIM > > var,
                        SAMRAI::tbox::Pointer< SAMRAI::hier::Patch< NDIM > > patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAI::tbox::Pointer< SAMRAI::hier::PatchLevel< NDIM > > patch_level = SAMRAI::tbox::Pointer< SAMRAI::hier::PatchLevel< NDIM > >(NULL)	 
                      );	




    }; // CartGridBodyForce






}// namespace IBTK




#endif //#ifndef included_CartGridBodyForce
