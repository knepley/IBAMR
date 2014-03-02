// Filename: LMarkerCoarsen.cpp
// Created on 30 Sep 2006 by Boyce Griffith
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

#include <vector>

#include "SAMRAI/hier/Index.h"
#include "LMarkerCoarsen.h"
#include "SAMRAI/hier/Patch.h"
#include "ibtk/LMarkerSet.h"
#include "ibtk/LMarkerSetData.h"
#include "ibtk/LMarkerSetVariable.h"
#include "ibtk/LMarker-inl.h"
#include "ibtk/LSetData.h"
#include "ibtk/LSet-inl.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace SAMRAI {
namespace hier {
class Variable;
}  // namespace hier
}  // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string LMarkerCoarsen::s_op_name = "LMARKER_COARSEN";

namespace
{
static const int COARSEN_OP_PRIORITY = 0;
static const int COARSEN_OP_STENCIL_WIDTH = 0;

inline int
coarsen(
    const int index,
    const int ratio)
{
    return (index < 0 ? (index+1)/ratio-1 : index/ratio);
}// coarsen

inline Index
coarsen_index(
    const Index& i,
    const IntVector& ratio)
{
    const Dimension& dim = i.getDim();
    Index coarse_i(dim);
    for (unsigned int d = 0; d < dim.getValue(); ++d)
    {
        coarse_i(d) = coarsen(i(d), ratio(d));
    }
    return coarse_i;
}// coarsen_index
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

LMarkerCoarsen::LMarkerCoarsen()
    : CoarsenOperator(s_op_name)
{
    // intentionally blank
    return;
}// LMarkerCoarsen

LMarkerCoarsen::~LMarkerCoarsen()
{
    // intentionally blank
    return;
}// ~LMarkerCoarsen

int
LMarkerCoarsen::getOperatorPriority() const
{
    return COARSEN_OP_PRIORITY;
}// getOperatorPriority

IntVector
LMarkerCoarsen::getStencilWidth(
    const Dimension& dim) const
{
    return IntVector(dim,COARSEN_OP_STENCIL_WIDTH);
}// getStencilWidth

void
LMarkerCoarsen::coarsen(
    Patch& coarse,
    const Patch& fine,
    const int dst_component,
    const int src_component,
    const Box& coarse_box,
    const IntVector& ratio) const
{
    auto dst_mark_data = BOOST_CAST<LMarkerSetData>(coarse.getPatchData(dst_component));
    auto src_mark_data = BOOST_CAST<LMarkerSetData>(fine  .getPatchData(src_component));

    const Box fine_box = Box::refine(coarse_box,ratio);
    for (auto    it = LMarkerSetData::SetIterator(*src_mark_data,true),
             it_end = LMarkerSetData::SetIterator(*src_mark_data,false);
         it != it_end; ++it)
    {
        const Index& fine_i = it.getIndex();
        const Index coarse_i = coarsen_index(fine_i,ratio);
        if (fine_box.contains(fine_i) && coarse_box.contains(coarse_i))
        {
            const auto& fine_mark_set = it->getDataSet();
            if (!dst_mark_data->isElement(coarse_i))
            {
                dst_mark_data->appendItemPointer(coarse_i, new LMarkerSet());
            }
            auto& coarse_mark_set = dst_mark_data->getItem(coarse_i)->getDataSet();
            coarse_mark_set.insert(coarse_mark_set.end(), fine_mark_set.begin(), fine_mark_set.end());
        }
    }
    return;
}// coarsen

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
