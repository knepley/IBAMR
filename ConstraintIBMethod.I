// Filename: ConstraintIBMethod.I
// Written by amneet bhalla on meghnad@mech.northwestern.edu
// Created on 12/07/2011.

// This is a concrete class which implements fast and efficient distributed Lagrange multiplier (fictitious domain) method
// for immersed structures.
     
// References 
//  *  Patankar et al. A new formulation of the distributed Lagrange multiplier/fictitious domain method 
//     for particulate flows. Int. Journal of Multiphase flows, 26, 1509-1524 (2000).
//  *  Shirgaonkar et al. A new mathematical formulation and fast algorithm for fully resolved simulation of self-propulsion.
//     JCP, 228 , 2366-2390 (2009).


/////////////////////////////////////////// INCLUDES //////////////////////////////
#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif


/////////////////////////////////////////// NAMESPACE /////////////////////////////

namespace IBAMR
{

/////////////////////////////////////////// PUBLIC ///////////////////////////////

inline void
ConstraintIBMethod::setINSCycleNumberAndCounter(
    const int cycle_num)
{

    d_INS_current_cycle_num = cycle_num;
    if(cycle_num = (d_INS_num_cycles -1))
        ++d_timestep_counter;      
 
    return;

} //setINSCycleNumberAndCounter


inline void
ConstraintIBMethod::setFurmorpTime(
    const double new_time)
{

    d_FuRMoRP_time = new_time;
    return;
}










} //IBAMR
