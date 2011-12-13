// Filename: ConstraintIBKinematics.I
// Written by amneet bhalla on meghnad@mech.northwestern.edu
// Created on 12/12/2011.

// This is an abstract class which encapsulates structure information and provides kinematics and 
// updated shape of the structure to ConstraintIBMethod class.
     

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

/////////////////////////////////// Class StructureParameters ////////////////////
int
ConstraintIBKinematics::StructureParameters::getTaggedPtIdx() const
{
    return d_tagged_pt_idx;

} //getTaggedPtIdx


bool
ConstraintIBKinematics::StructureParameters::getStructureIsSelfTranslating() const
{
    return d_struct_is_self_translating;

} //getStructureIsSelfTranslating


bool
ConstraintIBKinematics::StructureParameters::getStructureIsSelfRotating() const
{
    return d_struct_is_self_rotating;

} //getStructureIsSelfRotating


int
ConstraintIBKinematics::StructureParameters::getCoarsestLevelNumber() const
{
    return d_coarsest_ln;

} //getCoarsestLevelNumber


int
ConstraintIBKinematics::StructureParameters::getFinestLevelNumber() const
{
    return d_finest_ln;

} //getFinestLevelNumber
    

const std::vector<std::pair<int,int> >& 
ConstraintIBKinematics::StructureParameters::getLagIdxRange() const
{
    return d_idx_range;
 
} //getLagIdxRange


int 
ConstraintIBKinematics::StructureParameters::getTotalNodes() const
{
    return d_total_nodes;

} //getTotalNodes

std::string
ConstraintIBKinematics::StructureParameters::getPositionUpdateMethod() const
{

    return d_lag_position_update_method; 

} //getPositionUpdateMethod

/////////////////////////////////// Class ConstraintIBKinematics /////////////////

const ConstraintIBKinematics::StructureParameters&
ConstraintIBKinematics::getStructureParameters() const
{
    return d_struct_param;

} //getStructureParameters



} //IBAMR