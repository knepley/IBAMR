// Filename: MovingPlateKinematics.h
// Written by amneet bhalla on meghnad@mech.northwestern.edu
// Created on 12/20/2011.

// This is a concrete class which provides kinematics of Stokes' plate to ConstraintIBMethod class.


#include "MovingPlateKinematics.h"

/////////////////////////////////// INCLUDES /////////////////////////////////////



//SAMRAI INCLUDES
#include <tbox/Utilities.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/PIO.h>


//IBAMR INCLUDES
#include <ibamr/namespaces.h>

//IBTK INCLUDES


//C++ INCLUDES
#include <string>

namespace IBAMR
{
  
namespace
{
  
  
  
  
} //namespace anonymous

MovingPlateKinematics::MovingPlateKinematics(
    const std::string& object_name,
    Pointer<Database> input_db,
    LDataManager* l_data_manager,
    bool register_for_restart )
    : ConstraintIBKinematics(object_name,input_db,l_data_manager,register_for_restart),
      d_new_kinematics_vel(NDIM),
      d_new_shape()
{
    
    // NOTE: Parent class constructor registers class with the restart manager, sets object name. 
    
    const double vel_plate = input_db->getDoubleWithDefault("vel_plate",1.0);
    
    //Set the size of vectors.
    const StructureParameters& struct_param = getStructureParameters();
    const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    const int finest_ln   = struct_param.getFinestLevelNumber();
    const int total_nodes = struct_param.getTotalNodes();
    TBOX_ASSERT(coarsest_ln == finest_ln);
    
    for(int d = 0; d < NDIM; ++d)
    {
        d_new_kinematics_vel[d].resize(total_nodes);	
    }
    
    //fill in new and current velocity.
    const int rank       = SAMRAI_MPI::getRank();
    const int processors = SAMRAI_MPI::getNodes();
    
    //Fill the data arrays in a round-robin(interleaved/cyclic) fashion.
    for(int i = rank; i < total_nodes ; i += processors)
    {
         d_new_kinematics_vel[0][i]  = vel_plate;
    }
    SAMRAI_MPI::sumReduction(&d_new_kinematics_vel[0][0],total_nodes);
    
    return;
  
} //MovingPlateKinematics


MovingPlateKinematics::~MovingPlateKinematics()
{
    //intentionally left blank
    return;
  
}//~MovingPlateKinematics


void
MovingPlateKinematics::setNewKinematicsVelocity(
     const double /*Time*/,
     const std::vector<double>& /*incremented_angle_from_reference_axis*/,
     const std::vector<double>& /*center_of_mass*/,
     const std::vector<double>& /*tagged_pt_position*/)
{
  
    //intentionally  left blank
    return;
  
}//setNewKinematicsVelocity


const std::vector<std::vector<double> >&
MovingPlateKinematics::getNewKinematicsVelocity(
    const int level) const
{
  
    TBOX_ASSERT(level == getStructureParameters().getFinestLevelNumber());
        
    return d_new_kinematics_vel;

    
} //getNewKinematicsVelocity


const std::vector<std::vector<double> >&
MovingPlateKinematics::getCurrentKinematicsVelocity(const int level) const
{
    TBOX_ASSERT(level == getStructureParameters().getFinestLevelNumber());
    return d_new_kinematics_vel;
      
} //getCurrentKinematicsVelocity

void
MovingPlateKinematics::setNewShape()
{
    //intentionally left blank
    return;
  
} //setNewShape

const std::vector<std::vector<double> >&
MovingPlateKinematics::getNewShape(const int /*level*/) const
{
  
    return d_new_shape;
  
} //getNewShape




} //namespace IBAMR
