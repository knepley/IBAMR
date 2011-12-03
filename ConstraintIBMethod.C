// Filename: ConstraintIBMethod.h
// Written by amneet bhalla on meghnad@mech.northwestern.edu
// Created on 12/01/2011.

// This is a concrete class which implements fast and efficient distributed Lagrange multiplier (fictitious domain) method
// for immersed structures.
     
// References 
//  *  Patankar et al. A new formulation of the distributed Lagrange multiplier/fictitious domain method 
//     for particulate flows. Int. Journal of Multiphase flows, 26, 1509-1524 (2000).
//  *  Shirgaonkar et al. A new mathematical formulation and fast algorithm for fully resolved simulation of self-propulsion.
//     JCP, 228 , 2366-2390 (2009).


#include "ConstraintIBMethod.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES


// C++ INCLUDES

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{

  
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

ConstraintIBMethod::ConstraintIBMethod(
    const std::string& object_name,  
    Pointer< Database> input_db,
    bool register_for_restart,
    Pointer< INSHierarchyIntegrator > ins_hier_integrator) 
    : IBMethod(object_name, input_db, register_for_restart)
{
    // NOTE: Parent class constructor registers class with the restart manager, sets object name. 
    
    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);
  
  
} // ConstraintIBMethod



/////////////////////////////// PRIVATE ///////////////////////////////////////
void 
ConstraintIBMethod::getFromInput(
    Pointer< Database > input_db, 
    const bool from_restart)
{

} //getFromInput





}

