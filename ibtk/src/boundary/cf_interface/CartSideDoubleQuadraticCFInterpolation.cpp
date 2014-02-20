// Filename: CartSideDoubleQuadraticCFInterpolation.cpp
// Created on 30 Apr 2008 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "SAMRAI/hier/BoundaryBox.h"
#include "CartSideDoubleQuadraticCFInterpolation.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/geom/CartesianSideDoubleConservativeLinearRefine.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "IBTK_config.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define SC_QUAD_TANGENTIAL_INTERPOLATION_FC IBTK_FC_FUNC(scquadtangentialinterpolation2d,SCQUADTANGENTIALINTERPOLATION2D)
#define SC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(scquadnormalinterpolation2d,SCQUADNORMALINTERPOLATION2D)
#endif
#if (NDIM == 3)
#define SC_QUAD_TANGENTIAL_INTERPOLATION_FC IBTK_FC_FUNC(scquadtangentialinterpolation3d,SCQUADTANGENTIALINTERPOLATION3D)
#define SC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(scquadnormalinterpolation3d,SCQUADNORMALINTERPOLATION3D)
#endif

// Function interfaces
extern "C"
{
    void
    SC_QUAD_TANGENTIAL_INTERPOLATION_FC(
        double* U_fine0, double* U_fine1,
#if (NDIM == 3)
        double* U_fine2,
#endif
        const int& U_fine_gcw,
        const double* U_crse0, const double* U_crse1,
#if (NDIM == 3)
        const double* U_crse2,
#endif
        const int& U_crse_gcw,
        const int* sc_indicator0, const int* sc_indicator1,
#if (NDIM == 3)
        const int* sc_indicator2,
#endif
        const int& sc_indicator_gcw,
        const int& ilowerf0, const int& iupperf0,
        const int& ilowerf1, const int& iupperf1,
#if (NDIM == 3)
        const int& ilowerf2, const int& iupperf2,
#endif
        const int& ilowerc0, const int& iupperc0,
        const int& ilowerc1, const int& iupperc1,
#if (NDIM == 3)
        const int& ilowerc2, const int& iupperc2,
#endif
        const int& loc_index, const int* ratio_to_coarser,
        const int* blower, const int* bupper);

    void
    SC_QUAD_NORMAL_INTERPOLATION_FC(
        double* U0, double* U1,
#if (NDIM == 3)
        double* U2,
#endif
        const int& U_gcw,
        const double* W0, const double* W1,
#if (NDIM == 3)
        const double* W2,
#endif
        const int& W_gcw,
        const int* sc_indicator0, const int* sc_indicator1,
#if (NDIM == 3)
        const int* sc_indicator2,
#endif
        const int& sc_indicator_gcw,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
#if (NDIM == 3)
        const int& ilower2, const int& iupper2,
#endif
        const int& loc_index, const int* ratio_to_coarser,
        const int* blower, const int* bupper);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 1;
static const int GHOST_WIDTH_TO_FILL     = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideDoubleQuadraticCFInterpolation::CartSideDoubleQuadraticCFInterpolation()
    : d_patch_data_indices(),
      d_refine_op(new CartesianSideDoubleConservativeLinearRefine()),
      d_hierarchy(),
      d_cf_boundary(),
      d_sc_indicator_var(boost::make_shared<SideVariable<int> >(Dimension(NDIM), "CartSideDoubleQuadraticCFInterpolation::sc_indicator_var"))
{
    // Setup scratch variables.
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<VariableContext> context = var_db->getContext("CartSideDoubleQuadraticCFInterpolation::CONTEXT");
    if (var_db->checkVariableExists(d_sc_indicator_var->getName()))
    {
        d_sc_indicator_var = boost::dynamic_pointer_cast<SideVariable<int> >(var_db->getVariable(d_sc_indicator_var->getName()));
        TBOX_ASSERT(d_sc_indicator_var);
        d_sc_indicator_idx = var_db->mapVariableAndContextToIndex(d_sc_indicator_var, context);
        TBOX_ASSERT(d_sc_indicator_idx >= 0);
    }
    else
    {
        d_sc_indicator_idx = var_db->registerVariableAndContext(d_sc_indicator_var, context, IntVector(Dimension(NDIM), GHOST_WIDTH_TO_FILL));
    }
    return;
}// CartSideDoubleQuadraticCFInterpolation

CartSideDoubleQuadraticCFInterpolation::~CartSideDoubleQuadraticCFInterpolation()
{
    clearPatchHierarchy();
    return;
}// ~CartSideDoubleQuadraticCFInterpolation

void
CartSideDoubleQuadraticCFInterpolation::setPhysicalBoundaryConditions(
    Patch& /*patch*/,
    const double /*fill_time*/,
    const IntVector& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
}// setPhysicalBoundaryConditions

IntVector
CartSideDoubleQuadraticCFInterpolation::getRefineOpStencilWidth(
    const Dimension& dim) const
{
    TBOX_ASSERT(d_refine_op->getStencilWidth(dim).max() <= REFINE_OP_STENCIL_WIDTH);
    return IntVector(dim, REFINE_OP_STENCIL_WIDTH);
}// getRefineOpStencilWidth

void
CartSideDoubleQuadraticCFInterpolation::preprocessRefine(
    Patch& /*fine*/,
    const Patch& /*coarse*/,
    const Box& /*fine_box*/,
    const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
}// preprocessRefine

void
CartSideDoubleQuadraticCFInterpolation::postprocessRefine(
    Patch& fine,
    const Patch& coarse,
    const Box& fine_box,
    const IntVector& ratio)
{
    TBOX_ASSERT(d_hierarchy);

    // Ensure that the fine patch is located on the expected destination level;
    // if not, we are not guaranteed to have appropriate coarse-fine interface
    // boundary box information.
    if (!fine.inHierarchy())
    {
        const Dimension& dim = fine.getDim();
        const int ndim = dim.getValue();
        std::vector<BoxContainer> side_boxes(ndim);
        for (int d = 0; d < ndim; ++d)
        {
            side_boxes[d] = BoxContainer(SideGeometry::toSideBox(fine_box,d));
        }
        SideOverlap overlap(side_boxes, Transformation(IntVector::getZero(dim)));
        for (const auto& patch_data_index : d_patch_data_indices)
        {
            d_refine_op->refine(fine, coarse, patch_data_index, patch_data_index, overlap, ratio);
        }
        return;
    }

    // Get the co-dimension 1 cf boundary boxes.
    const GlobalId& patch_id = fine.getGlobalId();
    const int fine_patch_level_num = fine.getPatchLevelNumber();
    const std::vector<BoundaryBox>& cf_bdry_codim1_boxes = d_cf_boundary[fine_patch_level_num]->getBoundaries(patch_id, 1);
    if (cf_bdry_codim1_boxes.size() == 0) return;

    // Get the patch data.
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        auto fdata = BOOST_CAST<SideData<double> >(fine  .getPatchData(patch_data_index));
        auto cdata = BOOST_CAST<SideData<double> >(coarse.getPatchData(patch_data_index));
        auto indicator_data = BOOST_CAST<SideData<int> >(fine.getPatchData(d_sc_indicator_idx));
        TBOX_ASSERT(fdata);
        TBOX_ASSERT(cdata);
        TBOX_ASSERT(indicator_data);
        TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
        const int U_fine_ghosts = (fdata->getGhostCellWidth()).max();
        const int U_crse_ghosts = (cdata->getGhostCellWidth()).max();
        const int indicator_ghosts = (indicator_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_fine_ghosts != (fdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideDoubleQuadraticCFInterpolation::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
        if (U_crse_ghosts != (cdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideDoubleQuadraticCFInterpolation::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
#endif
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).max() == GHOST_WIDTH_TO_FILL);
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).min() == GHOST_WIDTH_TO_FILL);
        const int data_depth = fdata->getDepth();
        const IntVector ghost_width_to_fill(fdata->getDim(), GHOST_WIDTH_TO_FILL);
        auto pgeom_fine = BOOST_CAST<CartesianPatchGeometry>(fine.getPatchGeometry());
        TBOX_ASSERT(pgeom_fine);
        const Box& patch_box_fine = fine.getBox();
        const Box& patch_box_crse = coarse.getBox();
        for (int k = 0; k < cf_bdry_codim1_boxes.size(); ++k)
        {
            const BoundaryBox& bdry_box = cf_bdry_codim1_boxes[k];
            const Box bc_fill_box = pgeom_fine->getBoundaryFillBox(bdry_box, patch_box_fine, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            const int* const indicator0 = indicator_data->getPointer(0);
            const int* const indicator1 = indicator_data->getPointer(1);
#if (NDIM == 3)
            const int* const indicator2 = indicator_data->getPointer(2);
#endif
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U_fine0 = fdata->getPointer(0,depth);
                double* const U_fine1 = fdata->getPointer(1,depth);
#if (NDIM == 3)
                double* const U_fine2 = fdata->getPointer(2,depth);
#endif
                const double* const U_crse0 = cdata->getPointer(0,depth);
                const double* const U_crse1 = cdata->getPointer(1,depth);
#if (NDIM == 3)
                const double* const U_crse2 = cdata->getPointer(2,depth);
#endif
                SC_QUAD_TANGENTIAL_INTERPOLATION_FC(
                    U_fine0, U_fine1,
#if (NDIM == 3)
                    U_fine2,
#endif
                    U_fine_ghosts,
                    U_crse0, U_crse1,
#if (NDIM == 3)
                    U_crse2,
#endif
                    U_crse_ghosts,
                    indicator0, indicator1,
#if (NDIM == 3)
                    indicator2,
#endif
                    indicator_ghosts,
                    patch_box_fine.lower(0), patch_box_fine.upper(0),
                    patch_box_fine.lower(1), patch_box_fine.upper(1),
#if (NDIM == 3)
                    patch_box_fine.lower(2), patch_box_fine.upper(2),
#endif
                    patch_box_crse.lower(0), patch_box_crse.upper(0),
                    patch_box_crse.lower(1), patch_box_crse.upper(1),
#if (NDIM == 3)
                    patch_box_crse.lower(2), patch_box_crse.upper(2),
#endif
                    location_index, &ratio(0),
                    &bc_fill_box.lower(0), &bc_fill_box.upper(0));
            }
        }
    }
    return;
}// postprocessRefine

void
CartSideDoubleQuadraticCFInterpolation::setPatchDataIndex(
    const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
}// setPatchDataIndex

void
CartSideDoubleQuadraticCFInterpolation::setPatchDataIndices(
    const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
}// setPatchDataIndices

void
CartSideDoubleQuadraticCFInterpolation::setPatchDataIndices(
    const ComponentSelector& patch_data_indices)
{
    std::set<int> patch_data_index_set;
    for (int l = 0; l < patch_data_indices.getSize(); ++l)
    {
        if (patch_data_indices.isSet(l))
        {
            const int patch_data_index = l;
            patch_data_index_set.insert(patch_data_index);
        }
    }
    setPatchDataIndices(patch_data_index_set);
    return;
}// setPatchDataIndices

void
CartSideDoubleQuadraticCFInterpolation::setPatchHierarchy(
    const boost::shared_ptr<PatchHierarchy>& hierarchy)
{
    TBOX_ASSERT(hierarchy);
    if (d_hierarchy) clearPatchHierarchy();
    d_hierarchy = hierarchy;
    const int finest_level_number = d_hierarchy->getFinestLevelNumber();
    d_cf_boundary.resize(finest_level_number+1);
    const IntVector max_ghost_width(hierarchy->getDim(), GHOST_WIDTH_TO_FILL);
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        d_cf_boundary[ln] = boost::make_shared<CoarseFineBoundary>(*d_hierarchy, ln, max_ghost_width);
    }

    boost::shared_ptr<RefineAlgorithm> refine_alg;
    boost::shared_ptr<RefineOperator > refine_op;
    refine_alg->registerRefine(d_sc_indicator_idx, d_sc_indicator_idx, d_sc_indicator_idx, refine_op);
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        boost::shared_ptr<PatchLevel> level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_sc_indicator_idx))
        {
            level->allocatePatchData(d_sc_indicator_idx,0.0);
        }
        else
        {
            level->setTime(0.0,d_sc_indicator_idx);
        }
        for (const auto& patch : *level)
        {
            auto sc_indicator_data = boost::dynamic_pointer_cast<SideData<int> >(patch->getPatchData(d_sc_indicator_idx));
            sc_indicator_data->fillAll(0,sc_indicator_data->getGhostBox());
            sc_indicator_data->fillAll(1,sc_indicator_data->getBox()     );
        }
        refine_alg->createSchedule(level)->fillData(0.0);
    }
    return;
}// setPatchHierarchy

void
CartSideDoubleQuadraticCFInterpolation::clearPatchHierarchy()
{
    d_hierarchy.reset();
    d_cf_boundary.clear();
    return;
}// clearPatchHierarchy

void
CartSideDoubleQuadraticCFInterpolation::computeNormalExtension(
    Patch& patch,
    const IntVector& ratio,
    const IntVector& /*ghost_width_to_fill*/)
{
    TBOX_ASSERT(d_hierarchy);

    // Ensure that the fine patch is located on the expected destination level;
    // if not, we are not guaranteed to have appropriate coarse-fine interface
    // boundary box information.
    if (!patch.inHierarchy())
    {
        return;
    }
    
    // Get the co-dimension 1 cf boundary boxes.
    const GlobalId& patch_id = patch.getGlobalId();
    const int patch_level_num = patch.getPatchLevelNumber();
    const std::vector<BoundaryBox >& cf_bdry_codim1_boxes = d_cf_boundary[patch_level_num]->getBoundaries(patch_id, 1);
    const int n_cf_bdry_codim1_boxes = cf_bdry_codim1_boxes.size();

    // Check to see if there are any co-dimension 1 coarse-fine boundary boxes
    // associated with the patch; if not, there is nothing to do.
    if (n_cf_bdry_codim1_boxes == 0) return;

    // Get the patch data.
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        auto data = BOOST_CAST<SideData<double> >(patch.getPatchData(patch_data_index));
        TBOX_ASSERT(data);
        SideData<double> data_copy(data->getBox(), data->getDepth(), data->getGhostCellWidth());
        data_copy.copyOnBox(*data,data->getGhostBox());
        auto indicator_data = BOOST_CAST<SideData<int> >(patch.getPatchData(d_sc_indicator_idx));
        TBOX_ASSERT(data);
        TBOX_ASSERT(indicator_data);
        const int U_ghosts = (data->getGhostCellWidth()).max();
        const int W_ghosts = (data_copy.getGhostCellWidth()).max();
        const int indicator_ghosts = (indicator_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_ghosts != (data->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideDoubleQuadraticCFInterpolation::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
        if (W_ghosts != (data_copy.getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideDoubleQuadraticCFInterpolation::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
#endif
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).max() == GHOST_WIDTH_TO_FILL);
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).min() == GHOST_WIDTH_TO_FILL);
        const int data_depth = data->getDepth();
        const IntVector ghost_width_to_fill(data->getDim(), GHOST_WIDTH_TO_FILL);
        auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());
        TBOX_ASSERT(pgeom);
        const Box& patch_box = patch.getBox();
        for (const auto& bdry_box : cf_bdry_codim1_boxes)
        {
            const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            const int* const indicator0 = indicator_data->getPointer(0);
            const int* const indicator1 = indicator_data->getPointer(1);
#if (NDIM == 3)
            const int* const indicator2 = indicator_data->getPointer(2);
#endif
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U0 = data->getPointer(0,depth);
                double* const U1 = data->getPointer(1,depth);
#if (NDIM == 3)
                double* const U2 = data->getPointer(2,depth);
#endif
                const double* const W0 = data_copy.getPointer(0,depth);
                const double* const W1 = data_copy.getPointer(1,depth);
#if (NDIM == 3)
                const double* const W2 = data_copy.getPointer(2,depth);
#endif
                SC_QUAD_NORMAL_INTERPOLATION_FC(
                    U0, U1,
#if (NDIM == 3)
                    U2,
#endif
                    U_ghosts,
                    W0, W1,
#if (NDIM == 3)
                    W2,
#endif
                    W_ghosts,
                    indicator0, indicator1,
#if (NDIM == 3)
                    indicator2,
#endif
                    indicator_ghosts,
                    patch_box.lower(0), patch_box.upper(0),
                    patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
                    patch_box.lower(2), patch_box.upper(2),
#endif
                    location_index, &ratio(0),
                    &bc_fill_box.lower(0), &bc_fill_box.upper(0));
            }
        }
    }
    return;
}// computeNormalExtension

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
