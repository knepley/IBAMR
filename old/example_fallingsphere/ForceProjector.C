/******************************************************************************************************************************
   Filename: ForceProjector.C
   Written by amneet bhalla on Taylor@mech.northwestern.edu
   Created on 09/30/2011.
 
*******************************************************************************************************************************/


// APPLICATION INCLUDES
#include "ForceProjector.h"

// SAMRAI INCLUDES
#include <tbox/Array.h>
#include <VariableDatabase.h>
#include <PatchLevel.h>
#include <Patch.h>
#include <Box.h>

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/LNodeIndex.h>
#include <ibtk/LNodeIndexData.h>

// C++ INCLUDES


namespace IBTK
{
  
   namespace 
   {
     
     
     
   }// namespace anonymous
   
  void 
  callForceProjectorCallBackFunction(const double current_time,
			             const double new_time,
                                     const int cycle_num, 
				     void* ctx)
  {
    
    ForceProjector* forceprojector_ptr = static_cast<ForceProjector*>(ctx);
    forceprojector_ptr->calculateLagrangianBodyForce(new_time,current_time);
    forceprojector_ptr->calculateEulerianBodyForce(new_time, current_time);
    
    return;
    
    
  } // callForceProjectorCallBackFunction
  
  
  
  ForceProjector::ForceProjector(const std::string& object_name,
		                 LDataManager* lag_data_manager,
				 Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
		                 Pointer<Database> input_db)
                                :d_object_name(object_name),
                                 d_lag_data_manager(lag_data_manager),
                                 d_patch_hierarchy(patch_hierarchy),
                                 d_grav_const(NDIM)
  {
    
    // put some default values.
    d_rho_fluid = 1.0;
    d_rho_body  = 1.0;
    
    // Initialize  variables & variable contexts associated with Eulerian forces.
    VariableDatabase<NDIM>* var_db         = VariableDatabase<NDIM>::getDatabase();
    d_body_force_context                   = var_db->getContext(d_object_name + "::BODYFORCE");
    d_body_force_var                       = new SideVariable<NDIM,double>(d_object_name + "::BodyForce_var");
    d_sc_body_force_idx                    = var_db->registerVariableAndContext(d_body_force_var, d_body_force_context, 0);
    
    // Initialize variable associated with Lagrangian forces.
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    d_lag_force         = d_lag_data_manager->createLNodeLevelData(d_object_name + "::lag_force_data",
							           finest_ln,NDIM,true);
    
    getFromInput(input_db);
    
    return; 
    
  } // ctor
  
  ForceProjector::~ForceProjector()
  {
    
    //intentionally left blank
    return;
  }//dtor
   
  void
  ForceProjector::getFromInput(Pointer<Database> input_db)
  {
    
    d_rho_fluid   = input_db->getDoubleWithDefault("rho_fluid",d_rho_fluid);
    d_rho_body    = input_db->getDoubleWithDefault("rho_body",d_rho_body);
    d_grav_const  = input_db->getDoubleArray("gravitational_constant");
    
    
    
    return;
    
  }// getFromInput
  
  void
  ForceProjector::registerLagrangianQuantityName(const std::string& lag_quantity_name)
  {
     
    registerLagrangianQuantitiesName(std::vector<std::string>(1,lag_quantity_name));
    
    return;
    
  }//
  
  
  void
  ForceProjector::registerLagrangianQuantitiesName(const std::vector<std::string>& lag_quantities_name)
  {
    
    const unsigned size = lag_quantities_name.size();
    for(unsigned i = 0; i < size; ++i)
    {
      d_lag_quantities_name.push_back(lag_quantities_name[i]); 
    }
    
    return;
    
  }// registerLagrangianQuantitiesName
  
  void
  ForceProjector::associateVolumeElement(const double vol_lag_pt)
  {
    d_vol_lag_pt = vol_lag_pt;
    
    return;
  } // associateVolumeElement

  
  void
  ForceProjector::calculateLagrangianBodyForce(const double new_time, const double current_time)
  {
    
    // Get the patch data descriptor index for the LNodeIndexData.
    const int finest_ln              = d_patch_hierarchy->getFinestLevelNumber();
    const int lag_node_index_idx     = d_lag_data_manager->getLNodeIndexPatchDescriptorIndex();
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(finest_ln);
    
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
      Pointer<Patch<NDIM> > patch = level->getPatch(p());
      const Box<NDIM>& patch_box  = patch->getBox();
      const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
      
      for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
      {
         const LNodeIndex& node_idx    = *it;
         const int local_petsc_index   = node_idx.getLocalPETScIndex();
	 for(int depth = 0; depth < NDIM ; depth++)
         {
	   (*d_lag_force)(local_petsc_index,depth) = d_rho_fluid*d_grav_const[depth]*d_vol_lag_pt;
	 }
      }
	 
     }
    
    return;
    
  }// calculateLagrangianBodyForce
  
  void
  ForceProjector::calculateEulerianBodyForce(const double new_time, const double current_time)
  {
    
    // allocate patch data for Eulerian forcing.
     const int finest_ln   = d_patch_hierarchy->getFinestLevelNumber();
     for(int ln = 0; ln <= finest_ln; ++ln)
      {
	Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(ln);
       // deallocate the prior allocated data.
	if (level->checkAllocated(d_sc_body_force_idx) )       
	   level->deallocatePatchData(d_sc_body_force_idx);

	// allocate the data.
	level->allocatePatchData(d_sc_body_force_idx,current_time);
	for(PatchLevel<NDIM>::Iterator p(level); p ; p++)
        {
	    Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM,double> > bodyforce_data = patch->getPatchData(d_sc_body_force_idx);
	 
	    bodyforce_data->fill(0.0);
	 
         }//iterate over patches
	
	
      }
      
    // spread the lagrangian force from finest level to the finest level.
     std::vector< Pointer< LNodeLevelData > > F_data(finest_ln+1,  Pointer<LNodeLevelData>(NULL) );
     std::vector< Pointer< LNodeLevelData > > X_data(finest_ln+1,  Pointer<LNodeLevelData>(NULL) );


   // Fill in the above vectors at the finest level.
    F_data[finest_ln] = d_lag_force;   
    X_data[finest_ln] = d_lag_data_manager->getLNodeLevelData("X", finest_ln);
 
  // Spread the deformation velocities.
   d_lag_data_manager->spread(d_sc_body_force_idx ,
                              F_data,
			      X_data,
                              std::vector< Pointer< RefineSchedule< NDIM > > > (),
			      true,
			      true,
			      finest_ln,
			      finest_ln);
   

   return;
    
    
    
    return;
  }// calculateEulerianBodyForce

   
} // namespace IBTK