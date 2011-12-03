/******************************************************************************************************************************
   Filename : IBRigidStructureKinematics.C
   Written by amneet bhalla on Taylor@Courant Institute of Mathematical Sciences
   Created on 8/20/2011.


   This is a concrete class which sets the deformation velocity to zero and implements the shape update method of the abstract class 
   IBAMR::IBKinematics for a rigid structure needed by IBAMR::RigidityConstraint class.
 
*******************************************************************************************************************************/

////////////////////////////////////////// INCLUDES ////////////////////////////////////////////////

// APPLICATION INCLUDES
#include "IBRigidStructureKinematics.h"

// SAMRAI INCLUDES
#include <tbox/Utilities.h>


// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES


// C++ INCLUDES
#include <sstream>
#include <fstream>
#include <vector>

namespace IBAMR
{
  
    namespace 
    {
      static const bool DISCARD_COMMENTS = false;
      
     inline std::string
     discard_comments(const std::string& input_string)
     {
        // Create a copy of the input string, but without any text following a '!',
        // '#', or '%' character.
        std::string output_string = input_string;
        std::istringstream string_stream;

        // Discard any text following a '!' character.
        string_stream.str(output_string);
        std::getline(string_stream, output_string, '!');
        string_stream.clear();

        // Discard any text following a '#' character.
        string_stream.str(output_string);
        std::getline(string_stream, output_string, '#');
        string_stream.clear();

        // Discard any text following a '%' character.
        string_stream.str(output_string);
        std::getline(string_stream, output_string, '%');
        string_stream.clear();
        return output_string;
      }// discard_comments
      
    }//namespace anonymous
  
   IBRigidStructureKinematics::IBRigidStructureKinematics(Pointer< Database > input_db, 
				      Pointer< PatchHierarchy<NDIM> > patch_hierarchy)
   {

      // get the name of the object from the input database.
      d_object_name = input_db->getStringWithDefault("body_name","Rigid_Body");
      
      // get the name of file which has coordinates of the sphere.
      std::fstream coord_file_stream;
      std::string dir_name = input_db->getStringWithDefault("vertex_file_dir","./");
      std::string filename = input_db->getString( "vertex_filename");    
      if(filename.empty())
      {
	TBOX_ERROR(d_object_name << "ERROR :: base file containing coordinates of the rigid structure not found\n\n");
      }
      else
      {
	filename = dir_name + filename;
	coord_file_stream.open(filename.c_str(),std::fstream::in);
      }
     
     // total no. of material points.
     int num_material_pts;
     int global_lag_idx = -1;
     std::string line_string;
     std::vector<double> vec_coords(NDIM);
     std::vector<double> vec_def_vel(NDIM,0.0);
     //vec_def_vel[0] = 0.5;
     std::vector<double> vec_com(NDIM,0.0);
     
     // The first line in the file should contain the no. of material points.
     if( !std::getline(coord_file_stream,line_string))
     {
       TBOX_ERROR("Premature end to input file encountered before line 1 of file " << filename << "\n"); 
     }
     else
     {
       if(DISCARD_COMMENTS) line_string = discard_comments(line_string);
       std::istringstream line_stream(line_string);
       if(!(line_stream >> num_material_pts))
       {
	  TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << coord_file_stream 
	               << " .No. of material points expected." << std::endl);
       }
     }
     
     //Read-in the coordinates of material points from line number 2 onwards
     for (int k = 0; k < num_material_pts ; ++k)
     {
       
       if (!std::getline(coord_file_stream, line_string))
       {
          TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line  "<<  k+2 << " of file " <<  coord_file_stream << std::endl);
       }
       else
       {
           if(DISCARD_COMMENTS) line_string = discard_comments(line_string);
           std::istringstream line_stream(line_string);

           for(int dim = 0; dim < NDIM; ++dim)
	   {
	        if ( !(line_stream >> vec_coords[dim]) )
                {
                  TBOX_ERROR(d_object_name << ":\n Invalid entry in input file encountered on line " << k+2 << " of file " << coord_file_stream 
                             << " .Material point position component " << dim <<  "expected." << std::endl);
                }
	   
	   }
	    
	   // calculate a running sum of material points. 
	   for(int dim = 0; dim < NDIM; ++dim)
	   {
	       vec_com[dim] += vec_coords[dim];
	   }
	   
       }
       
       //insert the velocity in the map.
       d_map_new_shape[++global_lag_idx] = vec_coords;
       d_map_kinematics_vel[global_lag_idx] = vec_def_vel;
       
     }
     
     // find the com of the structure.
     for(int dim = 0; dim < NDIM; ++dim)
     {
	 vec_com[dim] /= num_material_pts;
     }
     
     // shift the structure to the origin.
     for(std::map<int, std::vector<double> >::iterator mitr = d_map_new_shape.begin(); mitr != d_map_new_shape.end();
	 ++mitr)
     {
	for(int dim = 0; dim < NDIM ; ++dim)
	{
	   (mitr->second)[dim] -= vec_com[dim];
	}
     }
     
     
     return;
     
   }//IBRigidStructureKinematics ctor

  void 
  IBRigidStructureKinematics::calculateGivenKinematicsVelocity(const double Time, const SAMRAI::tbox::Array<double>& incremented_angle_from_reference_axis, 
				        const std::vector<double>& X_COM, const std::vector<double>& X_tagged )
  {
      //intentionally left blank
      return;
    
  }//calculateGivenDeformationVelocity

  void 
  IBRigidStructureKinematics::updateLagrangianMarkerShapeAndOrientation(const SAMRAI::tbox::Array<double>& incremented_angle_from_reference_axis,
					                                const double new_time)
  {
      //intentionally left blank
      return;
  }//updateLagrangianMarkerShapeAndOrientation

  IBRigidStructureKinematics::~IBRigidStructureKinematics()
  {
    //intentionally left blank
    return;
  }//~IBRigidStructureKinematics

}//namespace IBAMR
