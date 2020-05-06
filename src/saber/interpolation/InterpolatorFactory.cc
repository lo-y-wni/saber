/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <memory>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/LocalConfiguration.h"

#include "oops/generic/Interpolator.h"
#include "oops/generic/InterpolatorFactory.h"
#include "saber/interpolation/InterpolatorBump.h"
#include "saber/interpolation/InterpolatorFactory.h"

namespace saber {

//-------------------------------------------------------------------------------------

std::unique_ptr<oops::Interpolator> InterpolatorFactory::create(
                                   eckit::LocalConfiguration & config,
                                   const atlas::FunctionSpace & grid1,
                                   const atlas::FunctionSpace & grid2,
                                   const atlas::field::FieldSetImpl * masks)
                                   const {
if (config.getString("interpolator") == "bump") {
    return std::unique_ptr<InterpolatorBump>(new InterpolatorBump(config, grid1, grid2, masks));
}

//-------------------------------------------------------------------------------------

}  // namespace saber
