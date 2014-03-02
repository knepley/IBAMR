// Filename: CartCellDoubleQuadraticCFInterpolation.cpp
// Created on 29 Oct 2007 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
#include "CartCellDoubleQuadraticCFInterpolation.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/geom/GridGeometry.h"
#include "IBTK_config.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/SAMRAI_config.h"
#include "boost/array.hpp"
#include "boost/array.hpp"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_QUAD_TANGENTIAL_INTERPOLATION_FC IBTK_FC_FUNC(ccquadtangentialinterpolation2d,CCQUADTANGENTIALINTERPOLATION2D)
#define CC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(ccquadnormalinterpolation2d,CCQUADNORMALINTERPOLATION2D)
#endif
#if (NDIM == 3)
#define CC_QUAD_TANGENTIAL_INTERPOLATION_FC IBTK_FC_FUNC(ccquadtangentialinterpolation3d,CCQUADTANGENTIALINTERPOLATION3D)
#define CC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(ccquadnormalinterpolation3d,CCQUADNORMALINTERPOLATION3D)
#endif

// Function interfaces
extern "C"
{
    void
    CC_QUAD_TANGENTIAL_INTERPOLATION_FC(
        double* U_fine, const int& U_fine_gcw,
        const double* U_coarse, const int& U_crse_gcw,
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
    CC_QUAD_NORMAL_INTERPOLATION_FC(
        double* U, const int& U_gcw,
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

CartCellDoubleQuadraticCFInterpolation::CartCellDoubleQuadraticCFInterpolation()
    : d_patch_data_indices(),
      d_refine_op(boost::make_shared<CartesianCellDoubleLinearRefine>()),
      d_hierarchy(),
      d_cf_boundary()
{
    // intentionally blank
    return;
}// CartCellDoubleQuadraticCFInterpolation

CartCellDoubleQuadraticCFInterpolation::~CartCellDoubleQuadraticCFInterpolation()
{
    clearPatchHierarchy();
    return;
}// ~CartCellDoubleQuadraticCFInterpolation

void
CartCellDoubleQuadraticCFInterpolation::setPhysicalBoundaryConditions(
    Patch& /*patch*/,
    const double /*fill_time*/,
    const IntVector& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
}// setPhysicalBoundaryConditions

IntVector
CartCellDoubleQuadraticCFInterpolation::getRefineOpStencilWidth(
    const Dimension& dim) const
{
    TBOX_ASSERT(d_refine_op->getStencilWidth(dim).max() <= REFINE_OP_STENCIL_WIDTH);
    return IntVector(dim, REFINE_OP_STENCIL_WIDTH);
}// getRefineOpStencilWidth

void
CartCellDoubleQuadraticCFInterpolation::preprocessRefine(
    Patch& /*fine*/,
    const Patch& /*coarse*/,
    const Box& /*fine_box*/,
    const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
}// preprocessRefine

void
CartCellDoubleQuadraticCFInterpolation::postprocessRefine(
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
        CellOverlap overlap(BoxContainer(fine_box), Transformation(IntVector::getZero(fine_box.getDim())));
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
        auto fdata = BOOST_CAST<CellData<double> >(fine  .getPatchData(patch_data_index));
        auto cdata = BOOST_CAST<CellData<double> >(coarse.getPatchData(patch_data_index));
        TBOX_ASSERT(fdata);
        TBOX_ASSERT(cdata);
        TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
        const int U_fine_ghosts = (fdata->getGhostCellWidth()).max();
        const int U_crse_ghosts = (cdata->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_fine_ghosts != (fdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleQuadraticCFInterpolation::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
        if (U_crse_ghosts != (cdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleQuadraticCFInterpolation::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
#endif
        const int data_depth = fdata->getDepth();
        const IntVector ghost_width_to_fill(fdata->getDim(), GHOST_WIDTH_TO_FILL);
        auto pgeom_fine = BOOST_CAST<CartesianPatchGeometry>(fine.getPatchGeometry());
        TBOX_ASSERT(pgeom_fine);
        const Box& patch_box_fine = fine.getBox();
        const Box& patch_box_crse = coarse.getBox();
        for (const auto& bdry_box : cf_bdry_codim1_boxes)
        {
            const Box bc_fill_box = pgeom_fine->getBoundaryFillBox(bdry_box, patch_box_fine, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U_fine = fdata->getPointer(depth);
                const double* const U_crse = cdata->getPointer(depth);
                CC_QUAD_TANGENTIAL_INTERPOLATION_FC(
                    U_fine, U_fine_ghosts,
                    U_crse, U_crse_ghosts,
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
CartCellDoubleQuadraticCFInterpolation::setPatchDataIndex(
    const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
}// setPatchDataIndex

void
CartCellDoubleQuadraticCFInterpolation::setPatchDataIndices(
    const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
}// setPatchDataIndices

void
CartCellDoubleQuadraticCFInterpolation::setPatchDataIndices(
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
CartCellDoubleQuadraticCFInterpolation::setPatchHierarchy(
    const boost::shared_ptr<PatchHierarchy>& hierarchy)
{
    TBOX_ASSERT(hierarchy);
    if (d_hierarchy) clearPatchHierarchy();
    d_hierarchy = hierarchy;
    const int finest_level_number = d_hierarchy->getFinestLevelNumber();
    d_cf_boundary.resize(finest_level_number+1);
    const IntVector& max_ghost_width = getRefineOpStencilWidth(hierarchy->getDim());
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        d_cf_boundary[ln] = boost::make_shared<CoarseFineBoundary>(*d_hierarchy, ln, max_ghost_width);
    }
    return;
}// setPatchHierarchy

void
CartCellDoubleQuadraticCFInterpolation::clearPatchHierarchy()
{
    d_hierarchy.reset();
    d_cf_boundary.clear();
    return;
}// clearPatchHierarchy

void
CartCellDoubleQuadraticCFInterpolation::computeNormalExtension(
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
    const std::vector<BoundaryBox>& cf_bdry_codim1_boxes = d_cf_boundary[patch_level_num]->getBoundaries(patch_id, 1);
    const int n_cf_bdry_codim1_boxes = cf_bdry_codim1_boxes.size();

    // Check to see if there are any co-dimension 1 coarse-fine boundary boxes
    // associated with the patch; if not, there is nothing to do.
    if (n_cf_bdry_codim1_boxes == 0) return;

    // Get the patch data.
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        auto data = BOOST_CAST<CellData<double> >(patch.getPatchData(patch_data_index));
        TBOX_ASSERT(data);
        const int U_ghosts = (data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_ghosts != (data->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleQuadraticCFInterpolation::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
#endif
        const int data_depth = data->getDepth();
        const IntVector ghost_width_to_fill(data->getDim(), GHOST_WIDTH_TO_FILL);
        auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());
        TBOX_ASSERT(pgeom);
        const Box& patch_box = patch.getBox();
        for (const auto& bdry_box : cf_bdry_codim1_boxes)
        {
            const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U = data->getPointer(depth);
                CC_QUAD_NORMAL_INTERPOLATION_FC(
                    U, U_ghosts,
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
