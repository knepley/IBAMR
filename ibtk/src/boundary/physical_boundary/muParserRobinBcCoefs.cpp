// Filename: muParserRobinBcCoefs.cpp
// Created on 25 Aug 2007 by Boyce Griffith
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

#include <iosfwd>
#include <ostream>
#include <sstream>
#include <utility>

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/SAMRAI_config.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "muParserRobinBcCoefs.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {
class Variable;
}  // namespace hier
}  // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int EXTENSIONS_FILLABLE = 128;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

muParserRobinBcCoefs::muParserRobinBcCoefs(
    const std::string& object_name,
    boost::shared_ptr<Database> input_db,
    boost::shared_ptr<CartesianGridGeometry > grid_geom)
    : d_grid_geom(grid_geom),
      d_constants(),
      d_acoef_function_strings(),
      d_bcoef_function_strings(),
      d_gcoef_function_strings(),
      d_acoef_parsers(2*NDIM),
      d_bcoef_parsers(2*NDIM),
      d_gcoef_parsers(2*NDIM),
      d_parser_time(new double),
      d_parser_posn(new Point)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
#else
    NULL_USE(object_name);
#endif
    // Read in user-provided constants.
    std::vector<std::string> db_key_names = input_db->getAllKeys();
    for (const auto& name : db_key_names)
    {
        if (input_db->isDouble(name))
        {
            d_constants[name] = input_db->getDouble(name);
        }
        else if (input_db->isFloat(name))
        {
            d_constants[name] = input_db->getFloat(name);
        }
        else if (input_db->isInteger(name))
        {
            d_constants[name] = input_db->getInteger(name);
        }
    }

    // Initialize the parsers with data read in from the input database.
    for (int d = 0; d < 2*NDIM; ++d)
    {
        std::string key_name;
        std::ostringstream stream;
        stream << "_function_" << d;
        const std::string postfix = stream.str();

        key_name = "acoef" + postfix;
        if (input_db->isString(key_name))
        {
            d_acoef_function_strings.push_back(input_db->getString(key_name));
        }
        else
        {
            d_acoef_function_strings.push_back("0.0");
            TBOX_WARNING("muParserRobinBcCoefs::muParserRobinBcCoefs():\n"
                         << "  no function corresponding to key ``" << key_name <<"'' found for side = " << d << "; using acoef = 0.0." << std::endl);
        }
        try
        {
            d_acoef_parsers[d].SetExpr(d_acoef_function_strings.back());
        }
        catch (mu::ParserError& e)
        {
            TBOX_ERROR("muParserRobinBcCoefs::setDataOnPatch():\n" <<
                       "  error: " << e.GetMsg()  << "\n" <<
                       "  in:    " << e.GetExpr() << "\n");
        }
        catch (...)
        {
            TBOX_ERROR("muParserRobinBcCoefs::setDataOnPatch():\n" <<
                       "  unrecognized exception generated by muParser library.\n");
        }

        key_name = "bcoef" + postfix;
        if (input_db->isString(key_name))
        {
            d_bcoef_function_strings.push_back(input_db->getString(key_name));
        }
        else
        {
            d_bcoef_function_strings.push_back("0.0");
            TBOX_WARNING("muParserRobinBcCoefs::muParserRobinBcCoefs():\n"
                         << "  no function corresponding to key ``" << key_name <<"'' found for side = " << d << "; using bcoef = 0.0." << std::endl);
        }
        try
        {
            d_bcoef_parsers[d].SetExpr(d_bcoef_function_strings.back());
        }
        catch (mu::ParserError& e)
        {
            TBOX_ERROR("muParserRobinBcCoefs::setDataOnPatch():\n" <<
                       "  error: " << e.GetMsg()  << "\n" <<
                       "  in:    " << e.GetExpr() << "\n");
        }
        catch (...)
        {
            TBOX_ERROR("muParserRobinBcCoefs::setDataOnPatch():\n" <<
                       "  unrecognized exception generated by muParser library.\n");
        }

        key_name = "gcoef" + postfix;
        if (input_db->isString(key_name))
        {
            d_gcoef_function_strings.push_back(input_db->getString(key_name));
        }
        else
        {
            d_gcoef_function_strings.push_back("0.0");
            TBOX_WARNING("muParserRobinBcCoefs::muParserRobinBcCoefs():\n"
                         << "  no function corresponding to key ``" << key_name <<"'' found for side = " << d << "; using gcoef = 0.0." << std::endl);
        }
        try
        {
            d_gcoef_parsers[d].SetExpr(d_gcoef_function_strings.back());
        }
        catch (mu::ParserError& e)
        {
            TBOX_ERROR("muParserRobinBcCoefs::setDataOnPatch():\n" <<
                       "  error: " << e.GetMsg()  << "\n" <<
                       "  in:    " << e.GetExpr() << "\n");
        }
        catch (...)
        {
            TBOX_ERROR("muParserRobinBcCoefs::setDataOnPatch():\n" <<
                       "  unrecognized exception generated by muParser library.\n");
        }
    }

    // Define the default and user-provided constants.
    std::vector<mu::Parser*> all_parsers(3*2*NDIM);
    for (int d = 0; d < 2*NDIM; ++d)
    {
        all_parsers[3*d  ] = &d_acoef_parsers[d];
        all_parsers[3*d+1] = &d_bcoef_parsers[d];
        all_parsers[3*d+2] = &d_gcoef_parsers[d];
    }
    const double pi = 3.1415926535897932384626433832795;
    const double* const xLower = grid_geom->getXLower();
    const double* const xUpper = grid_geom->getXUpper();
    for (auto it = begin(all_parsers); it != end(all_parsers); ++it)
    {
        // Various names for pi.
        (*it)->DefineConst("pi", pi);
        (*it)->DefineConst("Pi", pi);
        (*it)->DefineConst("PI", pi);

        // The extents of the domain.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            const std::string postfix = stream.str();

            (*it)->DefineConst("X_LOWER" + postfix, xLower[d]);
            (*it)->DefineConst("X_lower" + postfix, xLower[d]);
            (*it)->DefineConst("x_lower" + postfix, xLower[d]);
            (*it)->DefineConst("x_LOWER" + postfix, xLower[d]);
            (*it)->DefineConst("X_Lower" + postfix, xLower[d]);
            (*it)->DefineConst("X_lower" + postfix, xLower[d]);
            (*it)->DefineConst("XLower"  + postfix, xLower[d]);
            (*it)->DefineConst("Xlower"  + postfix, xLower[d]);
            (*it)->DefineConst("x_Lower" + postfix, xLower[d]);
            (*it)->DefineConst("x_lower" + postfix, xLower[d]);
            (*it)->DefineConst("xLower"  + postfix, xLower[d]);
            (*it)->DefineConst("xlower"  + postfix, xLower[d]);

            (*it)->DefineConst("X_LOWER_" + postfix, xLower[d]);
            (*it)->DefineConst("X_lower_" + postfix, xLower[d]);
            (*it)->DefineConst("x_lower_" + postfix, xLower[d]);
            (*it)->DefineConst("x_LOWER_" + postfix, xLower[d]);
            (*it)->DefineConst("X_Lower_" + postfix, xLower[d]);
            (*it)->DefineConst("X_lower_" + postfix, xLower[d]);
            (*it)->DefineConst("XLower_"  + postfix, xLower[d]);
            (*it)->DefineConst("Xlower_"  + postfix, xLower[d]);
            (*it)->DefineConst("x_Lower_" + postfix, xLower[d]);
            (*it)->DefineConst("x_lower_" + postfix, xLower[d]);
            (*it)->DefineConst("xLower_"  + postfix, xLower[d]);
            (*it)->DefineConst("xlower_"  + postfix, xLower[d]);

            (*it)->DefineConst("X_UPPER" + postfix, xUpper[d]);
            (*it)->DefineConst("X_upper" + postfix, xUpper[d]);
            (*it)->DefineConst("x_upper" + postfix, xUpper[d]);
            (*it)->DefineConst("x_UPPER" + postfix, xUpper[d]);
            (*it)->DefineConst("X_Upper" + postfix, xUpper[d]);
            (*it)->DefineConst("X_upper" + postfix, xUpper[d]);
            (*it)->DefineConst("XUpper"  + postfix, xUpper[d]);
            (*it)->DefineConst("Xupper"  + postfix, xUpper[d]);
            (*it)->DefineConst("x_Upper" + postfix, xUpper[d]);
            (*it)->DefineConst("x_upper" + postfix, xUpper[d]);
            (*it)->DefineConst("xUpper"  + postfix, xUpper[d]);
            (*it)->DefineConst("xupper"  + postfix, xUpper[d]);

            (*it)->DefineConst("X_UPPER_" + postfix, xUpper[d]);
            (*it)->DefineConst("X_upper_" + postfix, xUpper[d]);
            (*it)->DefineConst("x_upper_" + postfix, xUpper[d]);
            (*it)->DefineConst("x_UPPER_" + postfix, xUpper[d]);
            (*it)->DefineConst("X_Upper_" + postfix, xUpper[d]);
            (*it)->DefineConst("X_upper_" + postfix, xUpper[d]);
            (*it)->DefineConst("XUpper_"  + postfix, xUpper[d]);
            (*it)->DefineConst("Xupper_"  + postfix, xUpper[d]);
            (*it)->DefineConst("x_Upper_" + postfix, xUpper[d]);
            (*it)->DefineConst("x_upper_" + postfix, xUpper[d]);
            (*it)->DefineConst("xUpper_"  + postfix, xUpper[d]);
            (*it)->DefineConst("xupper_"  + postfix, xUpper[d]);
        }

        // User-provided constants.
        for (auto map_it = begin(d_constants); map_it != end(d_constants); ++map_it)
        {
            (*it)->DefineConst(map_it->first, map_it->second);
        }

        // Variables.
        (*it)->DefineVar("T", d_parser_time);
        (*it)->DefineVar("t", d_parser_time);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            const std::string postfix = stream.str();
            (*it)->DefineVar("X" + postfix, &d_parser_posn->data()[d]);
            (*it)->DefineVar("x" + postfix, &d_parser_posn->data()[d]);
            (*it)->DefineVar("X_" + postfix, &d_parser_posn->data()[d]);
            (*it)->DefineVar("x_" + postfix, &d_parser_posn->data()[d]);
        }
    }
    return;
}// muParserRobinBcCoefs

muParserRobinBcCoefs::~muParserRobinBcCoefs()
{
    delete d_parser_time;
    delete d_parser_posn;
    return;
}// ~muParserRobinBcCoefs

void
muParserRobinBcCoefs::setBcCoefs(
    const boost::shared_ptr<ArrayData<double> >& acoef_data,
    const boost::shared_ptr<ArrayData<double> >& bcoef_data,
    const boost::shared_ptr<ArrayData<double> >& gcoef_data,
    const boost::shared_ptr<Variable >& /*variable*/,
    const Patch& patch,
    const BoundaryBox& bdry_box,
    double fill_time) const
{
    const Dimension& dim = patch.getDim();
    const Box& patch_box = patch.getBox();
    const Index& patch_lower = patch_box.lower();
    auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch.getPatchGeometry());

    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    // Loop over the boundary box and set the coefficients.
    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis =  location_index / 2;
    const Box& bc_coef_box = (acoef_data ? acoef_data->getBox() :
                              bcoef_data ? bcoef_data->getBox() :
                              gcoef_data ? gcoef_data->getBox() :
                              Box(dim));

    const mu::Parser& acoef_parser = d_acoef_parsers[location_index];
    const mu::Parser& bcoef_parser = d_bcoef_parsers[location_index];
    const mu::Parser& gcoef_parser = d_gcoef_parsers[location_index];
    *d_parser_time = fill_time;
    for (const auto& i : bc_coef_box)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d != bdry_normal_axis)
            {
                (*d_parser_posn)[d] = x_lower[d] + dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
            }
            else
            {
                (*d_parser_posn)[d] = x_lower[d] + dx[d]*(static_cast<double>(i(d)-patch_lower(d)));
            }
        }
        try
        {
            if (acoef_data) (*acoef_data)(i,0) = acoef_parser.Eval();
            if (bcoef_data) (*bcoef_data)(i,0) = bcoef_parser.Eval();
            if (gcoef_data) (*gcoef_data)(i,0) = gcoef_parser.Eval();
        }
        catch (mu::ParserError& e)
        {
            TBOX_ERROR("muParserRobinBcCoefs::setDataOnPatch():\n" <<
                       "  error: " << e.GetMsg()  << "\n" <<
                       "  in:    " << e.GetExpr() << "\n");
        }
        catch (...)
        {
            TBOX_ERROR("muParserRobinBcCoefs::setDataOnPatch():\n" <<
                       "  unrecognized exception generated by muParser library.\n");
        }
    }
    return;
}// setBcCoefs

IntVector
muParserRobinBcCoefs::numberOfExtensionsFillable() const
{
    return IntVector(Dimension(NDIM), EXTENSIONS_FILLABLE);
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
