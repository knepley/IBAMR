// Filename: CartSideRobinPhysBdryOp.cpp
// Created on 21 May 2008 by Boyce Griffith
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

#include <map>
#include <ostream>
#include <string>

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "CartSideRobinPhysBdryOp.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "IBTK_config.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "boost/array.hpp"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "boost/array.hpp"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Utilities.h"

// FORTRAN ROUTINES
#define CC_ROBIN_PHYS_BDRY_OP_1_X_FC2D IBTK_FC_FUNC(ccrobinphysbdryop1x2d, CCROBINPHYSBDRYOP1X2D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Y_FC2D IBTK_FC_FUNC(ccrobinphysbdryop1y2d, CCROBINPHYSBDRYOP1Y2D)
#define SC_ROBIN_PHYS_BDRY_OP_1_X_FC2D IBTK_FC_FUNC(scrobinphysbdryop1x2d, SCROBINPHYSBDRYOP1X2D)
#define SC_ROBIN_PHYS_BDRY_OP_1_Y_FC2D IBTK_FC_FUNC(scrobinphysbdryop1y2d, SCROBINPHYSBDRYOP1Y2D)
#define SC_ROBIN_PHYS_BDRY_OP_2_FC2D IBTK_FC_FUNC(scrobinphysbdryop22d, SCROBINPHYSBDRYOP22D)

#define CC_ROBIN_PHYS_BDRY_OP_1_X_FC3D IBTK_FC_FUNC(ccrobinphysbdryop1x3d, CCROBINPHYSBDRYOP1X3D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Y_FC3D IBTK_FC_FUNC(ccrobinphysbdryop1y3d, CCROBINPHYSBDRYOP1Y3D)
#define CC_ROBIN_PHYS_BDRY_OP_1_Z_FC3D IBTK_FC_FUNC(ccrobinphysbdryop1z3d, CCROBINPHYSBDRYOP1Z3D)
#define SC_ROBIN_PHYS_BDRY_OP_1_X_FC3D IBTK_FC_FUNC(scrobinphysbdryop1x3d, SCROBINPHYSBDRYOP1X3D)
#define SC_ROBIN_PHYS_BDRY_OP_1_Y_FC3D IBTK_FC_FUNC(scrobinphysbdryop1y3d, SCROBINPHYSBDRYOP1Y3D)
#define SC_ROBIN_PHYS_BDRY_OP_1_Z_FC3D IBTK_FC_FUNC(scrobinphysbdryop1z3d, SCROBINPHYSBDRYOP1Z3D)
#define SC_ROBIN_PHYS_BDRY_OP_2_FC3D IBTK_FC_FUNC(scrobinphysbdryop23d, SCROBINPHYSBDRYOP22D)
#define CC_ROBIN_PHYS_BDRY_OP_2_FC3D IBTK_FC_FUNC(ccrobinphysbdryop23d, CCROBINPHYSBDRYOP23D)
#define SC_ROBIN_PHYS_BDRY_OP_3_FC3D IBTK_FC_FUNC(scrobinphysbdryop33d, SCROBINPHYSBDRYOP32D)

extern "C"
{
    void
    CC_ROBIN_PHYS_BDRY_OP_1_X_FC2D(
        double* u, const int& u_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& blower1, const int& bupper1,
        const double* dx, const int& adjoint_op);

    void
    CC_ROBIN_PHYS_BDRY_OP_1_Y_FC2D(
        double* u, const int& u_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& blower0, const int& bupper0,
        const double* dx, const int& adjoint_op);

    void
    SC_ROBIN_PHYS_BDRY_OP_1_X_FC2D(
        double* u0, const int& u_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& blower1, const int& bupper1,
        const double* dx, const int& adjoint_op);

    void
    SC_ROBIN_PHYS_BDRY_OP_1_Y_FC2D(
        double* u1, const int& u_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& blower0, const int& bupper0,
        const double* dx, const int& adjoint_op);

    void
    SC_ROBIN_PHYS_BDRY_OP_2_FC2D(
        double* u0, double* u1,
        const int& u_gcw,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
        const int& adjoint_op);

    void
    CC_ROBIN_PHYS_BDRY_OP_1_X_FC3D(
        double* u, const int& u_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower1, const int& bupper1,
        const int& blower2, const int& bupper2,
        const double* dx, const int& adjoint_op);

    void
    CC_ROBIN_PHYS_BDRY_OP_1_Y_FC3D(
        double* u, const int& u_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower0, const int& bupper0,
        const int& blower2, const int& bupper2,
        const double* dx, const int& adjoint_op);

    void
    CC_ROBIN_PHYS_BDRY_OP_1_Z_FC3D(
        double* U, const int& U_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
        const double* dx, const int& adjoint_op);

    void
    SC_ROBIN_PHYS_BDRY_OP_1_X_FC3D(
        double* u0, const int& u_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower1, const int& bupper1,
        const int& blower2, const int& bupper2,
        const double* dx, const int& adjoint_op);

    void
    SC_ROBIN_PHYS_BDRY_OP_1_Y_FC3D(
        double* u1, const int& u_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower0, const int& bupper0,
        const int& blower2, const int& bupper2,
        const double* dx, const int& adjoint_op);

    void
    SC_ROBIN_PHYS_BDRY_OP_1_Z_FC3D(
        double* u2, const int& u_gcw,
        const double* acoef, const double* bcoef, const double* gcoef,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
        const double* dx, const int& adjoint_op);

    void
    SC_ROBIN_PHYS_BDRY_OP_2_FC3D(
        double* u0, double* u1,
        double* u2,
        const int& u_gcw,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
        const int& blower2, const int& bupper2,
        const int& adjoint_op);

    void
    CC_ROBIN_PHYS_BDRY_OP_2_FC3D(
        double* U, const int& U_gcw,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
        const int& blower2, const int& bupper2,
        const int& adjoint_op);

    void
    SC_ROBIN_PHYS_BDRY_OP_3_FC3D(
        double* u0, double* u1, double* u2,
        const int& u_gcw,
        const int& location_index,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
        const int& ilower2, const int& iupper2,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
        const int& blower2, const int& bupper2,
        const int& adjoint_op);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 1;

inline Box
compute_tangential_extension(
    const Box& box,
    const int data_axis)
{
    Box extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
}// compute_tangential_extension
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp()
    : RobinPhysBdryPatchStrategy()
{
    // intentionally blank
    return;
}// CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp(
    const int patch_data_index,
    const std::vector<boost::shared_ptr<RobinBcCoefStrategy> >& bc_coefs,
    const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndex(patch_data_index);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp(
    const std::set<int>& patch_data_indices,
    const std::vector<boost::shared_ptr<RobinBcCoefStrategy> >& bc_coefs,
    const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::CartSideRobinPhysBdryOp(
    const ComponentSelector& patch_data_indices,
    const std::vector<boost::shared_ptr<RobinBcCoefStrategy> >& bc_coefs,
    const bool homogeneous_bc)
    : RobinPhysBdryPatchStrategy()
{
    setPatchDataIndices(patch_data_indices);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// CartSideRobinPhysBdryOp

CartSideRobinPhysBdryOp::~CartSideRobinPhysBdryOp()
{
    // intentionally blank
    return;
}// ~CartSideRobinPhysBdryOp

void
CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions(
    Patch& patch,
    const double fill_time,
    const IntVector& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector::getZero(patch.getDim())) return;

    const Dimension& dim = patch.getDim();
    const int ndim = dim.getValue();
#if !defined(NDEBUG)
    // Ensure the target patch data corresponds to a side centered variable and
    // that the proper number of boundary condition objects have been provided.
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        auto patch_data = boost::dynamic_pointer_cast<SideData<double> >(patch.getPatchData(patch_data_idx));
        if (!patch_data)
        {
            TBOX_ERROR("CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  patch data index " << patch_data_idx << " does not correspond to a side-centered double precision variable." << std::endl);
        }
        if (ndim*patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartSideRobinPhysBdryOp::setPhysicalBoundaryConditions():\n"
                       << "  data depth for patch data index " << patch_data_idx << " is " << patch_data->getDepth() << "\n"
                       << "  but " << d_bc_coefs.size() << " boundary condition coefficient objects were provided to the class constructor." << std::endl);
        }
    }
#endif

    // Set the boundary conditions along the co-dimension one boundary boxes,
    // then extrapolate those values to the co-dimension two and three boundary
    // boxes.
    static const bool adjoint_op = false;
    const std::vector<BoundaryBox> physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim1Normal(patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim1Transverse(patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    const std::vector<BoundaryBox> physical_codim2_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim2(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
    const std::vector<BoundaryBox> physical_codim3_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim3(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
    return;
}// setPhysicalBoundaryConditions

IntVector
CartSideRobinPhysBdryOp::getRefineOpStencilWidth(
    const Dimension& dim) const
{
    return IntVector(dim, REFINE_OP_STENCIL_WIDTH);
}// getRefineOpStencilWidth

void
CartSideRobinPhysBdryOp::accumulateFromPhysicalBoundaryData(
    Patch& patch,
    const double fill_time,
    const IntVector& ghost_width_to_fill)
{
    if (ghost_width_to_fill == IntVector::getZero(patch.getDim())) return;

    const Dimension& dim = patch.getDim();
    const int ndim = dim.getValue();
#if !defined(NDEBUG)
    // Ensure the target patch data corresponds to a side centered variable and
    // that the proper number of boundary condition objects have been provided.
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        auto patch_data = boost::dynamic_pointer_cast<SideData<double> >(patch.getPatchData(patch_data_idx));
        if (!patch_data)
        {
            TBOX_ERROR("CartSideRobinPhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  patch data index " << patch_data_idx << " does not correspond to a side-centered double precision variable." << std::endl);
        }
        if (ndim*patch_data->getDepth() != static_cast<int>(d_bc_coefs.size()))
        {
            TBOX_ERROR("CartSideRobinPhysBdryOp::accumulateFromPhysicalBoundaryData():\n"
                       << "  data depth for patch data index " << patch_data_idx << " is " << patch_data->getDepth() << "\n"
                       << "  but " << d_bc_coefs.size() << " boundary condition coefficient objects were provided to the class constructor." << std::endl);
        }
    }
#endif

    // Set the boundary conditions along the co-dimension one boundary boxes,
    // then extrapolate those values to the co-dimension two and three boundary
    // boxes.
    static const bool adjoint_op = true;
    const std::vector<BoundaryBox> physical_codim3_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim3(patch_data_idx, physical_codim3_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
    const std::vector<BoundaryBox> physical_codim2_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim2(patch_data_idx, physical_codim2_boxes, ghost_width_to_fill, patch, adjoint_op);
    }
    const std::vector<BoundaryBox> physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim1Transverse(patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    for (const auto& patch_data_idx : d_patch_data_indices)
    {
        fillGhostCellValuesCodim1Normal(patch_data_idx, physical_codim1_boxes, fill_time, ghost_width_to_fill, patch, adjoint_op);
    }
    return;
}// accumulateFromPhysicalBoundaryData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CartSideRobinPhysBdryOp::fillGhostCellValuesCodim1Normal(
    const int patch_data_idx,
    const std::vector<BoundaryBox>& physical_codim1_boxes,
    const double fill_time,
    const IntVector& ghost_width_to_fill,
    Patch& patch,
    const bool adjoint_op)
{
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    if (n_physical_codim1_boxes == 0) return;

    const Dimension& dim = patch.getDim();
    const int ndim = dim.getValue();
    const Box& patch_box = patch.getBox();
    auto pgeom = BOOST_CAST<CartesianPatchGeometry >(patch.getPatchGeometry());
    TBOX_ASSERT(pgeom);
    const double* const dx = pgeom->getDx();
    auto patch_data = BOOST_CAST<SideData<double> >(patch.getPatchData(patch_data_idx));
    TBOX_ASSERT(patch_data);
    const int patch_data_depth = patch_data->getDepth();
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable > var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartSideRobinPhysBdryOp::fillGhostCellValuesCodim1Normal():\n"
                   "  patch data for patch data index " << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector gcw_to_fill = IntVector::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    // Set the boundary condition coefficients and then set the ghost cell
    // values.
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox& bdry_box   = physical_codim1_boxes[n];
        const unsigned int location_index   = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), location_index);
        const Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
        auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
        auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
        auto gcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            auto bc_coef = d_bc_coefs[ndim*d+bdry_normal_axis];
            auto extended_bc_coef = boost::dynamic_pointer_cast<ExtendedRobinBcCoefStrategy>(bc_coef);
            if (extended_bc_coef)
            {
                extended_bc_coef->setTargetPatchDataIndex(patch_data_idx);
                extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
            }
            bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, patch, trimmed_bdry_box, fill_time);
            if (d_homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
            if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();
            if (location_index == 0 || location_index == 1)  // lower x, upper x
            {              
                if (ndim == 2)
                {
                    SC_ROBIN_PHYS_BDRY_OP_1_X_FC2D(
                        patch_data->getPointer(bdry_normal_axis,d), patch_data_gcw,
                        acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                        location_index,
                        patch_box.lower(0), patch_box.upper(0),
                        patch_box.lower(1), patch_box.upper(1),
                        bc_coef_box.lower(1), bc_coef_box.upper(1),
                        dx, adjoint_op ? 1 : 0);
                }
                else if (ndim == 3)
                {
                    SC_ROBIN_PHYS_BDRY_OP_1_X_FC3D(
                        patch_data->getPointer(bdry_normal_axis,d), patch_data_gcw,
                        acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                        location_index,
                        patch_box.lower(0), patch_box.upper(0),
                        patch_box.lower(1), patch_box.upper(1),
                        patch_box.lower(2), patch_box.upper(2),
                        bc_coef_box.lower(1), bc_coef_box.upper(1),
                        bc_coef_box.lower(2), bc_coef_box.upper(2),
                        dx, adjoint_op ? 1 : 0);
                }
                else
                {
                    TBOX_ERROR("unsupported dim\n");
                }
            }
            else if (location_index == 2 || location_index == 3)  // lower y, upper y
            {              
                if (ndim == 2)
                {
                    SC_ROBIN_PHYS_BDRY_OP_1_Y_FC2D(
                        patch_data->getPointer(bdry_normal_axis,d), patch_data_gcw,
                        acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                        location_index,
                        patch_box.lower(0), patch_box.upper(0),
                        patch_box.lower(1), patch_box.upper(1),
                        bc_coef_box.lower(0), bc_coef_box.upper(0),
                        dx, adjoint_op ? 1 : 0);
                }
                else if (ndim == 3)
                {
                    SC_ROBIN_PHYS_BDRY_OP_1_Y_FC3D(
                        patch_data->getPointer(bdry_normal_axis,d), patch_data_gcw,
                        acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                        location_index,
                        patch_box.lower(0), patch_box.upper(0),
                        patch_box.lower(1), patch_box.upper(1),
                        patch_box.lower(2), patch_box.upper(2),
                        bc_coef_box.lower(0), bc_coef_box.upper(0),
                        bc_coef_box.lower(2), bc_coef_box.upper(2),
                        dx, adjoint_op ? 1 : 0);
                }
                else
                {
                    TBOX_ERROR("unsupported dim\n");
                }
            }
            else if (location_index == 4 || location_index == 5)  // lower z, upper z
            {              
                if (ndim == 3)
                {
                    SC_ROBIN_PHYS_BDRY_OP_1_Z_FC3D(
                        patch_data->getPointer(bdry_normal_axis,d), patch_data_gcw,
                        acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                        location_index,
                        patch_box.lower(0), patch_box.upper(0),
                        patch_box.lower(1), patch_box.upper(1),
                        patch_box.lower(2), patch_box.upper(2),
                        bc_coef_box.lower(0), bc_coef_box.upper(0),
                        bc_coef_box.lower(1), bc_coef_box.upper(1),
                        dx, adjoint_op ? 1 : 0);
                }
                else
                {
                    TBOX_ERROR("unsupported dim\n");
                }
            }
            else
            {
                TBOX_ERROR("unsupported location index\n");
            }
        }
    }    
    return;
}// fillGhostCellValuesCodim1Normal
        
void
CartSideRobinPhysBdryOp::fillGhostCellValuesCodim1Transverse(
    const int patch_data_idx,
    const std::vector<BoundaryBox>& physical_codim1_boxes,
    const double fill_time,
    const IntVector& ghost_width_to_fill,
    Patch& patch,
    const bool adjoint_op)
{
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    if (n_physical_codim1_boxes == 0) return;

    const Dimension& dim = patch.getDim();
    const int ndim = dim.getValue();
    const Box& patch_box = patch.getBox();
    auto pgeom = BOOST_CAST<CartesianPatchGeometry >(patch.getPatchGeometry());
    TBOX_ASSERT(pgeom);
    const double* const dx = pgeom->getDx();
    const double* const patch_x_lower = pgeom->getXLower();
    const double* const patch_x_upper = pgeom->getXUpper();
    auto patch_data = BOOST_CAST<SideData<double> >(patch.getPatchData(patch_data_idx));
    TBOX_ASSERT(patch_data);
    const int patch_data_depth = patch_data->getDepth();
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable > var;
    var_db->mapIndexToVariable(patch_data_idx, var);
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartSideRobinPhysBdryOp::fillGhostCellValuesCodim1Transverse():\n"
                   "  patch data for patch data index " << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector gcw_to_fill = IntVector::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    std::vector<Box> side_box(ndim, Box(dim));
    for (unsigned int axis = 0; axis < ndim; ++axis)
    {
        side_box[axis] = SideGeometry::toSideBox(patch_box,axis);
    }

    const IntVector& ratio_to_level_zero = pgeom->getRatio();
    PatchGeometry::TwoDimBool touches_regular_bdry(dim);
    for (unsigned int axis = 0; axis < ndim; ++axis)
    {
        for (int upperlower = 0; upperlower < 2; ++upperlower)
        {
            touches_regular_bdry (axis,upperlower) = pgeom->getTouchesRegularBoundary(axis,upperlower);
        }
    }

    // Set the boundary condition coefficients and then set the ghost cell
    // values.
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox& bdry_box   = physical_codim1_boxes[n];
        const unsigned int location_index   = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), location_index);
        for (unsigned int axis = 0; axis < ndim; ++axis)
        {
            if (axis != bdry_normal_axis)
            {
                const Box bc_coef_box = compute_tangential_extension(PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);
                auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
                auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
                auto gcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);

                // Temporarily reset the patch geometry object associated with
                // the patch so that boundary conditions are set at the correct
                // spatial locations.
                boost::array<double,SAMRAI::MAX_DIM_VAL> shifted_patch_x_lower, shifted_patch_x_upper;
                for (unsigned int d = 0; d < ndim; ++d)
                {
                    shifted_patch_x_lower[d] = patch_x_lower[d];
                    shifted_patch_x_upper[d] = patch_x_upper[d];
                }
                shifted_patch_x_lower[axis] -= 0.5*dx[axis];
                shifted_patch_x_upper[axis] -= 0.5*dx[axis];
                patch.setPatchGeometry(boost::make_shared<CartesianPatchGeometry>(ratio_to_level_zero, touches_regular_bdry, dx, shifted_patch_x_lower.data(), shifted_patch_x_upper.data()));

                // Set the boundary condition coefficients.
                for (int d = 0; d < patch_data_depth; ++d)
                {
                    auto bc_coef = d_bc_coefs[ndim*d+axis];
                    auto extended_bc_coef = boost::dynamic_pointer_cast<ExtendedRobinBcCoefStrategy>(bc_coef);
                    if (extended_bc_coef)
                    {
                        extended_bc_coef->setTargetPatchDataIndex(patch_data_idx);
                        extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
                    }
                    bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, patch, trimmed_bdry_box, fill_time);
                    if (d_homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
                    if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();

                    // Restore the original patch geometry object.
                    patch.setPatchGeometry(pgeom);

                    // Set the boundary values.
                    if (location_index == 0 || location_index == 1)  // lower x, upper x
                    {
                        if (ndim == 2)
                        {
                            CC_ROBIN_PHYS_BDRY_OP_1_X_FC2D(
                                patch_data->getPointer(axis,d), patch_data_gcw,
                                acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                                location_index,
                                side_box[axis].lower(0), side_box[axis].upper(0),
                                side_box[axis].lower(1), side_box[axis].upper(1),
                                bc_coef_box.lower(1), bc_coef_box.upper(1),
                                dx, adjoint_op ? 1 : 0);
                        }
                        else if (ndim == 3)
                        {
                            CC_ROBIN_PHYS_BDRY_OP_1_X_FC3D(
                                patch_data->getPointer(axis,d), patch_data_gcw,
                                acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                                location_index,
                                side_box[axis].lower(0), side_box[axis].upper(0),
                                side_box[axis].lower(1), side_box[axis].upper(1),
                                side_box[axis].lower(2), side_box[axis].upper(2),
                                bc_coef_box.lower(1), bc_coef_box.upper(1),
                                bc_coef_box.lower(2), bc_coef_box.upper(2),
                                dx, adjoint_op ? 1 : 0);
                        }
                        else
                        {
                            TBOX_ERROR("unsupported dim\n");
                        }
                    }
                    else if (location_index == 2 || location_index == 3)  // lower y, upper y
                    {
                        if (ndim == 2)
                        {
                            CC_ROBIN_PHYS_BDRY_OP_1_Y_FC2D(
                                patch_data->getPointer(axis,d), patch_data_gcw,
                                acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                                location_index,
                                side_box[axis].lower(0), side_box[axis].upper(0),
                                side_box[axis].lower(1), side_box[axis].upper(1),
                                bc_coef_box.lower(0), bc_coef_box.upper(0),
                                dx, adjoint_op ? 1 : 0);
                        }
                        else if (ndim == 3)
                        {
                            CC_ROBIN_PHYS_BDRY_OP_1_Y_FC3D(
                                patch_data->getPointer(axis,d), patch_data_gcw,
                                acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                                location_index,
                                side_box[axis].lower(0), side_box[axis].upper(0),
                                side_box[axis].lower(1), side_box[axis].upper(1),
                                side_box[axis].lower(2), side_box[axis].upper(2),
                                bc_coef_box.lower(0), bc_coef_box.upper(0),
                                bc_coef_box.lower(2), bc_coef_box.upper(2),
                                dx, adjoint_op ? 1 : 0);
                        }
                        else
                        {
                            TBOX_ERROR("unsupported dim\n");
                        }
                    }
                    else if (location_index == 4 || location_index == 5)  // lower z, upper z
                    {
                        if (ndim == 3)
                        {
                            CC_ROBIN_PHYS_BDRY_OP_1_Z_FC3D(
                                patch_data->getPointer(axis,d), patch_data_gcw,
                                acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
                                location_index,
                                side_box[axis].lower(0), side_box[axis].upper(0),
                                side_box[axis].lower(1), side_box[axis].upper(1),
                                side_box[axis].lower(2), side_box[axis].upper(2),
                                bc_coef_box.lower(0), bc_coef_box.upper(0),
                                bc_coef_box.lower(1), bc_coef_box.upper(1),
                                dx, adjoint_op ? 1 : 0);
                        }
                        else
                        {
                            TBOX_ERROR("unsupported dim\n");
                        }
                    }                
                    else
                    {
                        TBOX_ERROR("unsupported location index\n");
                    }
                }
            }
        }
    }
    return;
}// fillGhostCellValuesCodim1Transverse

void
CartSideRobinPhysBdryOp::fillGhostCellValuesCodim2(
    const int patch_data_idx,
    const std::vector<BoundaryBox>& physical_codim2_boxes,
    const IntVector& ghost_width_to_fill,
    const Patch& patch,
    const bool adjoint_op)
{
    const int n_physical_codim2_boxes = physical_codim2_boxes.size();
    if (n_physical_codim2_boxes == 0) return;

    const Dimension& dim = patch.getDim();
    const int ndim = dim.getValue();
    const Box& patch_box = patch.getBox();
    auto pgeom = BOOST_CAST<CartesianPatchGeometry >(patch.getPatchGeometry());
    TBOX_ASSERT(pgeom);
    auto patch_data = BOOST_CAST<SideData<double> >(patch.getPatchData(patch_data_idx));
    TBOX_ASSERT(patch_data);
    const int patch_data_depth = patch_data->getDepth();
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartSideRobinPhysBdryOp::fillGhostCellValuesCodim2():\n"
                   "  patch data for patch data index " << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector gcw_to_fill = IntVector::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    std::vector<Box> side_box(ndim, Box(dim));
    for (unsigned int axis = 0; axis < ndim; ++axis)
    {
        side_box[axis] = SideGeometry::toSideBox(patch_box,axis);
    }

    for (int n = 0; n < n_physical_codim2_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim2_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            if (ndim == 2)
            {
                SC_ROBIN_PHYS_BDRY_OP_2_FC2D(
                    patch_data->getPointer(0,d),
                    patch_data->getPointer(1,d),
                    patch_data_gcw,
                    location_index,
                    patch_box.lower(0), patch_box.upper(0),
                    patch_box.lower(1), patch_box.upper(1),
                    bc_fill_box.lower(0), bc_fill_box.upper(0),
                    bc_fill_box.lower(1), bc_fill_box.upper(1),
                    adjoint_op ? 1 : 0);
            }
            else if (ndim == 3)
            {
                SC_ROBIN_PHYS_BDRY_OP_2_FC3D(
                    patch_data->getPointer(0,d),
                    patch_data->getPointer(1,d),
                    patch_data->getPointer(2,d),
                    patch_data_gcw,
                    location_index,
                    patch_box.lower(0), patch_box.upper(0),
                    patch_box.lower(1), patch_box.upper(1),
                    patch_box.lower(2), patch_box.upper(2),
                    bc_fill_box.lower(0), bc_fill_box.upper(0),
                    bc_fill_box.lower(1), bc_fill_box.upper(1),
                    bc_fill_box.lower(2), bc_fill_box.upper(2),
                    adjoint_op ? 1 : 0);
                
                if (location_index < 4)
                {
                    static const unsigned int axis = 0;
                    CC_ROBIN_PHYS_BDRY_OP_2_FC3D(
                        patch_data->getPointer(axis,d), patch_data_gcw,
                        location_index,
                        side_box[axis].lower(0), side_box[axis].upper(0),
                        side_box[axis].lower(1), side_box[axis].upper(1),
                        side_box[axis].lower(2), side_box[axis].upper(2),
                        bc_fill_box.lower(0), bc_fill_box.upper(0),
                        bc_fill_box.lower(1), bc_fill_box.upper(1),
                        bc_fill_box.lower(2), bc_fill_box.upper(2),
                        adjoint_op ? 1 : 0);
                }
                else if (location_index >= 4 && location_index < 8)
                {
                    static const unsigned int axis = 1;
                    CC_ROBIN_PHYS_BDRY_OP_2_FC3D(
                        patch_data->getPointer(axis,d), patch_data_gcw,
                        location_index,
                        side_box[axis].lower(0), side_box[axis].upper(0),
                        side_box[axis].lower(1), side_box[axis].upper(1),
                        side_box[axis].lower(2), side_box[axis].upper(2),
                        bc_fill_box.lower(0), bc_fill_box.upper(0),
                        bc_fill_box.lower(1), bc_fill_box.upper(1),
                        bc_fill_box.lower(2), bc_fill_box.upper(2),
                        adjoint_op ? 1 : 0);
                }
                else if (location_index >= 8 && location_index < 12)
                {
                    static const unsigned int axis = 2;
                    CC_ROBIN_PHYS_BDRY_OP_2_FC3D(
                        patch_data->getPointer(axis,d), patch_data_gcw,
                        location_index,
                        side_box[axis].lower(0), side_box[axis].upper(0),
                        side_box[axis].lower(1), side_box[axis].upper(1),
                        side_box[axis].lower(2), side_box[axis].upper(2),
                        bc_fill_box.lower(0), bc_fill_box.upper(0),
                        bc_fill_box.lower(1), bc_fill_box.upper(1),
                        bc_fill_box.lower(2), bc_fill_box.upper(2),
                        adjoint_op ? 1 : 0);
                }
                else
                {
                    TBOX_ERROR("unsupported location index\n");
                }
            }
            else
            {
                TBOX_ERROR("unsupported dim\n");
            }            
        }
    }
    return;
}// fillGhostCellValuesCodim2

void
CartSideRobinPhysBdryOp::fillGhostCellValuesCodim3(
    const int patch_data_idx,
    const std::vector<BoundaryBox>& physical_codim3_boxes,
    const IntVector& ghost_width_to_fill,
    const Patch& patch,
    const bool adjoint_op)
{
    const int n_physical_codim3_boxes = physical_codim3_boxes.size();
    if (n_physical_codim3_boxes == 0) return;

    const Dimension& dim = patch.getDim();
    const int ndim = dim.getValue();
    const Box& patch_box = patch.getBox();
    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());
    TBOX_ASSERT(pgeom);
    auto patch_data = BOOST_CAST<SideData<double> >(patch.getPatchData(patch_data_idx));
    TBOX_ASSERT(patch_data);
    const int patch_data_depth = patch_data->getDepth();
    const int patch_data_gcw = (patch_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
    if (patch_data_gcw != (patch_data->getGhostCellWidth()).min())
    {
        TBOX_ERROR("CartSideRobinPhysBdryOp::fillGhostCellValuesCodim3():\n"
                   "  patch data for patch data index " << patch_data_idx << " does not have uniform ghost cell widths." << std::endl);
    }
#endif
    const IntVector gcw_to_fill = IntVector::min(patch_data->getGhostCellWidth(), ghost_width_to_fill);

    for (int n = 0; n < n_physical_codim3_boxes; ++n)
    {
        const BoundaryBox& bdry_box = physical_codim3_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        for (int d = 0; d < patch_data_depth; ++d)
        {
            if (ndim == 3)
            {
                SC_ROBIN_PHYS_BDRY_OP_3_FC3D(
                    patch_data->getPointer(0,d),
                    patch_data->getPointer(1,d),
                    patch_data->getPointer(2,d),
                    patch_data_gcw,
                    location_index,
                    patch_box.lower(0), patch_box.upper(0),
                    patch_box.lower(1), patch_box.upper(1),
                    patch_box.lower(2), patch_box.upper(2),
                    bc_fill_box.lower(0), bc_fill_box.upper(0),
                    bc_fill_box.lower(1), bc_fill_box.upper(1),
                    bc_fill_box.lower(2), bc_fill_box.upper(2),
                    adjoint_op ? 1 : 0);
            }
            else
            {
                TBOX_ERROR("unsupported dim\n");
            }
        }
    }
    return;
}// fillGhostCellValuesCodim3

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
