/******************************************************************************************************************************
   Filename: ForceProjector.h
   Written by amneet bhalla on Taylor@mech.northwestern.edu
   Created on 09/30/2011.
 
*******************************************************************************************************************************/

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_ForceProjector
#define included_ForceProjector

///////////////////////////// INCLUDES ///////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Pointer.h>
#include <tbox/Database.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <PatchHierarchy.h>

// IBAMR INCLUDES


// IBTK INCLUDES
#include <ibtk/LDataManager.h>
#include <ibtk/LNodeLevelData.h>

// C++ INCLUDES
#include <string>

namespace IBTK
{
 
  
  /*!
   * Pre processing call back function to be hooked into IBAMR::INSStaggeredHierachyIntegrator class.
   * 
   * \param current_time is the time at t_n.
   * \param new_time is the time at t_n+1 = t_n + dt.
   * \param cycle_num is the cycle of predictor-corrector scheme.
   * \param ctx is the pointer to IBTK::ForceProjector class object.
   */
  
  

  void 
  callForceProjectorCallBackFunction(const double current_time,
			             const double new_time,
                                     const int cycle_num, 
				     void* ctx);
  
  
  class ForceProjector
  {
    
    /*!
     * \brief Class ForceProjector is a utility class which projects force from
     * Lagrangian points onto the background mesh.
     * 
     */
  public:
    
    /*!
     * The only constructor of this class.
     */
    ForceProjector(const std::string& object_name,
		   IBTK::LDataManager* lag_data_manager,
		   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
		   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
  
    /*!
     * Destructor for this class.
     */
     ~ForceProjector();
        
    /*!
     * Register the name of Lagrangian quantities to be used to calculate forces on Lagrangian
     * points.
     */
    void
    registerLagrangianQuantityName(const std::string& lag_quantity_name);
    
    void
    registerLagrangianQuantitiesName(const std::vector<std::string>& lag_quantities_name);
    
    /*!
     * Register volume associated with each material point.
     */
    void
    associateVolumeElement(const double vol_lag_pt);
    /*!
     * Calculate forcing on Lagrangian points.
     */
    void
    calculateLagrangianBodyForce(const double new_time, const double current_time);
    
    /*!
     * Spread the Lagrangian forcing on the background mesh.
     */
    void
    calculateEulerianBodyForce(const double new_time, const double current_time);
    
   /*!
    *  Get the patch index associated with Eulerian force.
    */
   int
   getEulerianForcePatchDataIndex() const
   {
     
     return d_sc_body_force_idx;
     
   }// getEulerianForcePatchDataIndex
    
    
    
   //////////////// PRIVATE /////////////////////////////
      
  private:
    
   /*!
    * Default ctor is not implemented and should not be used.
    */
    ForceProjector();
    
   /*!
    * Default assignment operator is not implemented and should not be used.
    */
    ForceProjector& 
    operator = (const ForceProjector& that);
    
   /*!
    * Default copy ctor is not implemented and should not be used. 
    */
    ForceProjector(const ForceProjector& from);
    
   /*!
    * Get the values from input_db.
    */
    void
    getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
   
  
   /*!
    * Name of this object.
    */
    std::string d_object_name;
    
   /*!
    * Pointer to LDataManager.
    */
    IBTK::LDataManager* d_lag_data_manager;
    
   /*!
    * Pointer to Patch Hierarchy.
    */ 
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_patch_hierarchy;
    
   /*!
    * Pointer to Lagrangian force data.
    */
   SAMRAI::tbox::Pointer< IBTK::LNodeLevelData > d_lag_force;
    
   /*!
    *  Variables and variable context associated with calculating Eulerian force.
    */
    SAMRAI::tbox::Pointer< SAMRAI::pdat::SideVariable< NDIM,double > > d_body_force_var;
    SAMRAI::tbox::Pointer< SAMRAI::hier::VariableContext >             d_body_force_context;
    int d_sc_body_force_idx;
    
   /*!
    * Name of Lagrangian quantities to be used in calculating forces.
    */
    std::vector<std::string> d_lag_quantities_name;
   
   /*!
    * Volume associated with each element.
    */
    double d_vol_lag_pt;
   
   /*!
    * Densities of fluid and body.
    */
    double d_rho_fluid, d_rho_body;
    
  /*!
   * Gravitational force constants.
   */
   SAMRAI::tbox::Array<double> d_grav_const;
   
  
    
  }; //class ForceProjector
  
  
  
  
}// namespace IBTK

#endif // #ifndef included_ForceProjector