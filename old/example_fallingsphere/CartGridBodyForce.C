/****************************************************************************************************************

  // Filename: CartGridBodyForce.h
  // Created by amneet bhalla on 08/11/2011 Taylor @ Courant Institute Of Mathematical Sciences
 
****************************************************************************************************************/

#include "CartGridBodyForce.h"

////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <VariableDatabase.h>
#include <tbox/Utilities.h>
#include <SideVariable.h>
#include <SideData.h>
#include <tbox/PIO.h>


// IBAMR INCLUDES
#include <ibamr/namespaces.h>



// IBTK INCLUDES


// C++ INCLUDES


namespace IBTK
{
  
  
  
    namespace 
    {
    
    }// namespace anonymous
  
    CartGridBodyForce::CartGridBodyForce(const int body_force_idx)
                                     : d_sc_BodyForce_idx(body_force_idx)
    {
        // this is intentionally left blank
        return;
      
    }//CartGridBodyForce
    
    bool
    CartGridBodyForce::isTimeDependent() const
    {
        return true;
    }//isTimeDependent


    void
    CartGridBodyForce::setDataOnPatch(const int data_idx, 
                                      Pointer< Variable<NDIM> > var, 
	                              Pointer< Patch<NDIM> > patch, 
	                              const double data_time, 
				      const bool initial_time, 
				      Pointer< PatchLevel<NDIM> > patch_level)
    {

      Pointer< SideVariable<NDIM,double> > copyto_side_var = var;
      
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!patch.isNull());
      TBOX_ASSERT(!copyto_side_var.isNull());
#endif 
      
      Pointer< SideData<NDIM, double> > copyto_sidedata   = patch->getPatchData(data_idx);
      const Pointer< SideData<NDIM, double> > copyfrom_sidedata = patch->getPatchData(d_sc_BodyForce_idx);
      
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!copyto_sidedata.isNull());
      TBOX_ASSERT(!copyfrom_sidedata.isNull());
#endif
      
      copyto_sidedata->copy(*copyfrom_sidedata);
      
      
      
      return;
    }//setDataonPatch

  
  
}// namespace IBTK
