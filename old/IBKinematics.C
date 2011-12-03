/******************************************************************************************************************************
   Filename : IBKinematics.C
   Written by amneet bhalla on meghnad@mech.northwestern.edu
   Created on 10/06/2011.
 
*******************************************************************************************************************************/

// APPLICATION INCLUDES
#include "IBKinematics.h"

//IBAMR INCLUDES
#include <ibamr/namespaces.h>


namespace IBAMR
{
  
  IBKinematics::IBKinematics()
  {
    //intentionally left blank
    return;
  }

  IBKinematics::~IBKinematics()
  {
    //intentionally left blank
    return;
  }
 
  
} //namespace IBAMR


/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

 #include <tbox/Pointer.C>
template class Pointer<IBAMR::IBKinematics>;

//////////////////////////////////////////////////////////////////////////////
