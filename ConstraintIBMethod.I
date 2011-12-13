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
    if(cycle_num == (d_INS_num_cycles -1))
        ++d_timestep_counter;      
 
    return;

} //setINSCycleNumberAndCounter


inline void
ConstraintIBMethod::setFurmorpTime(
    const double current_time,
    const double new_time)
{
    d_FuRMoRP_current_time = current_time;
    d_FuRMoRP_new_time     = new_time;

    return;
} //setFurmorpTime

inline void
ConstraintIBMethod::createNewLagrangianWorkspace()
{
  
    const int coarsest_ln = 0;
    const int finest_ln   = d_hierarchy->getFinestLevelNumber();
    d_l_data_U_interp      .resize(finest_ln+1);
    d_l_data_U_correction  .resize(finest_ln+1);
    d_l_data_U_new         .resize(finest_ln+1);
    d_l_data_U_half        .resize(finest_ln+1);
    d_l_data_X_new_Euler   .resize(finest_ln+1);
    d_l_data_X_new_MidPoint.resize(finest_ln+1);
    if (d_INS_current_cycle_num == 0)  d_l_data_U_current.resize(finest_ln+1);
    
    for(int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if(!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
	d_l_data_U_interp[ln]       = d_l_data_manager->createLData(d_object_name + "interp_lag_vel", ln, NDIM, false);  
        d_l_data_U_correction[ln]   = d_l_data_manager->createLData(d_object_name + "correct_lag_vel", ln, NDIM, false);  
        d_l_data_U_new[ln]          = d_l_data_manager->createLData(d_object_name + "new_lag_vel", ln, NDIM, false); 
        d_l_data_U_half[ln]         = d_l_data_manager->createLData(d_object_name + "half_lag_vel", ln, NDIM, false); 
	d_l_data_X_new_Euler[ln]    = d_l_data_manager->createLData(d_object_name + "X_Euler", ln, NDIM, false); 
	d_l_data_X_new_MidPoint[ln] = d_l_data_manager->createLData(d_object_name + "X_MidPoint", ln, NDIM, false); 
        if (d_INS_current_cycle_num == 0) 
	   d_l_data_U_current[ln]   = d_l_data_manager->createLData(d_object_name + "current_lag_vel", ln, NDIM, false); 

    }
    return;
} //createLagrangianWorkspace


inline void
ConstraintIBMethod::destroyPreviousLagrangianWorkspace()
{
    d_l_data_U_interp      .clear();
    d_l_data_U_correction  .clear();
    d_l_data_U_new         .clear();
    d_l_data_U_half        .clear();
    d_l_data_X_new_Euler   .clear();
    d_l_data_X_new_MidPoint.clear();
    
    return;
} //destroyLagrangianWorkspace

inline void
ConstraintIBMethod::interpolateFluidSolveVelocity()
{
  
    const int coarsest_ln = 0;
    const int finest_ln   = d_hierarchy->getFinestLevelNumber();
    
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > F_data(finest_ln+1, SAMRAI::tbox::Pointer<IBTK::LData>(NULL) );
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > X_data(finest_ln+1, SAMRAI::tbox::Pointer<IBTK::LData>(NULL) );
  
    for(int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if(!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
	F_data[ln] = d_l_data_U_interp[ln];  
	X_data[ln] = d_l_data_manager->getLData("X",ln);
    }
  
    d_l_data_manager->interp(d_u_fluidSolve_idx,F_data,X_data,
        d_ib_hier_integrator->getCoarsenSchedules(d_object_name+"SYNC::u_fluidSolve"),
	d_ib_hier_integrator->getGhostfillRefineSchedules(d_object_name+"FILL_GHOSTCELL::u_fluidSolve"),
	d_FuRMoRP_new_time); 
    
    return; 
  
} //interpolateFluidSolveVelocity

inline void
ConstraintIBMethod::spreadCorrectedLagrangianVelocity()
{
    const int coarsest_ln = 0;
    const int finest_ln   = d_hierarchy->getFinestLevelNumber();
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > F_data(finest_ln+1,  SAMRAI::tbox::Pointer<IBTK::LData>(NULL) );
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > X_data(finest_ln+1,  SAMRAI::tbox::Pointer<IBTK::LData>(NULL) );
    
    for(int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if(!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
	F_data[ln] = d_l_data_manager->getLData(d_object_name + "correct_lag_vel", ln);  
	X_data[ln] = d_l_data_manager->getLData("X",ln);
    }
  
    d_l_data_manager->spread(d_u_fluidSolve_idx,F_data,X_data,
        d_ib_hier_integrator->getProlongRefineSchedules(d_object_name+"PROLONG::u_fluidSolve"),
	true,true); 
     
    return;
  
} //spreadCorrectedLagrangianVelocity

inline void
ConstraintIBMethod::synchronizeLevels()
{
    const int coarsest_ln = 0;
    const int finest_ln   = d_hierarchy->getFinestLevelNumber();
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& coarsen_schedules = d_ib_hier_integrator->getCoarsenSchedules(d_object_name+"SYNC::u_fluidSolve");
    
    // Do the Coarsening Operation.
    for(int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        coarsen_schedules[ln]->coarsenData();
    }
    
    return;
  
} //synchronizeLevels


inline void
ConstraintIBMethod::calculateMidPointVelocity()
{
  
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    int ierr;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecAXPBYPCZ(d_l_data_U_half[ln]->getVec(), 0.5, 0.5, 0.0, d_l_data_U_current[ln]->getVec(), d_l_data_U_new[ln]->getVec());  IBTK_CHKERRQ(ierr);
    }
    return;
  
}// calculateMidPointVelocity


} //IBAMR
