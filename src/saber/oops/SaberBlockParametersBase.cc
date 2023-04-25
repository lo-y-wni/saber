/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/SaberBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

bool SaberBlockParametersBase::doCalibration() const {
  const auto & readConf = readParams.value();
  const auto & calibrationConf = calibrationParams.value();
  if (readConf != boost::none && calibrationConf != boost::none) {
    ABORT("read and calibration configuration cannot be present at the same time in block " +
    saberBlockName.value());
  }
  const auto & blockName = saberBlockName.value();
  return (calibrationConf != boost::none || blockName == "Ensemble");
}

// -----------------------------------------------------------------------------

bool SaberBlockParametersBase::doRead() const {
  const auto & readConf = readParams.value();
  const auto & calibrationConf = calibrationParams.value();
  if (readConf != boost::none && calibrationConf != boost::none) {
    ABORT("read and calibration configuration cannot be present at the same time in block " +
    saberBlockName.value());
  }
  return (readConf != boost::none);
}

// -----------------------------------------------------------------------------

}  // namespace saber
