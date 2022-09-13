/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_SABERBLOCKTEST_H_
#define SABER_OOPS_SABERBLOCKTEST_H_

#include <memory>
#include <string>
#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>

#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> class SaberBlockTestParameters
  : public oops::ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(SaberBlockTestParameters, oops::ApplicationParameters)

 public:
  typedef typename oops::Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef typename oops::State<MODEL>::Parameters_    StateParameters_;

  /// Geometry parameters
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Tested variables
  oops::RequiredParameter<oops::Variables> variables{"variables", this};

  /// Background parameters
  oops::RequiredParameter<StateParameters_> background{"background", this};

  /// SABER blocks
  oops::RequiredParameter<std::vector<SaberBlockParametersWrapper>>
    saberBlocks{"saber blocks", this};

  /// Adjoint test tolerance
  oops::Parameter<double> adjointTolerance{"adjoint test tolerance", 1.0e-12, this};

  /// Inverse test tolerance
  oops::Parameter<double> inverseTolerance{"inverse test tolerance", 1.0e-12, this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class SaberBlockTest : public oops::Application {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef oops::State<MODEL>                              State_;

 public:
  static const std::string classname() {return "saber::SaberBlockTest";}
  explicit SaberBlockTest(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}

  virtual ~SaberBlockTest() {}

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    util::Timer timer(classname(), "execute");

    // Deserialize parameters
    SaberBlockTestParameters<MODEL> params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup geometry
    const Geometry_ geom(params.geometry, this->getComm());

    // Setup variables
    const oops::Variables vars(params.variables);

    // Setup background state
    const State_ xx(geom, params.background);
/*
    // Vector of FieldSets
    std::vector<atlas::FieldSet> fsvec;

    // Initialize first FieldSet
    Increment_ dx(geom, vars, xx.validTime());
    fsvec.push_back(dx.fieldset());

    // Loop over SABER blocks
    std::unique_ptr<SaberBlockBase_> saberBlock_;
    for (SaberBlockParametersWrapper_::reverse_iterator it = params.saberBlocks.value().rbegin();
      it != params.saberBlocks.value().rend(); ++it) {
      const SaberBlockParametersBase & saberBlockParams = it.saberBlockParameters;

      // Create block
      saberBlock_.reset(SaberBlockFactory<MODEL>::create(geom, saberBlockParams, xx, xx));

      // Adjoint test

      // Apply adjoint block
      saberBlock_->multiplyAD();

      // 
      Increment_ dx2TLAD(dx1TLAD);
      dx2TLAD.random();
      Increment_ dx2TLADsave(dx2TLAD, true);
  
      // Apply forward block
      saberBlock_->multiply(dx2TLAD.fieldSet());

      // ATLAS fieldset to Increment_
      dx1TLAD.synchronizeFieldsAD();
      dx2TLAD.synchronizeFields();

      // Compute adjoint test
      const double dp1 = dx1TLAD.dot_product_with(dx2TLADsave);
      const double dp2 = dx2TLAD.dot_product_with(dx1TLADsave);

      oops::Log::test() << "Adjoint test for block " << saberCentralBlock_->name() <<
                           ": y^t (Ax) = " << dp1 <<
                           ": x^t (Ay) = " << dp2 << std::endl;
      ASSERT(0.5*abs(dp1-dp2)/(dp1+dp2) < params.adjointTolerance.value());

      // Left inverse test

      // Apply block
      it->multiply(dx1Inv.fieldSet());
      it->inverseMultiply(dx1Inv.fieldSet());
      it->multiplyAD(dx2Inv.fieldSet());
      it->inverseMultiplyAD(dx2Inv.fieldSet());

      // ATLAS fieldset to Increment_
      dx1Inv.synchronizeFields();
      dx2Inv.synchronizeFields();

      // Compute inverse test
      dx1Inv -= dx1Invsave;
      dx2Inv -= dx2Invsave;
      const double dp1 = dx1Inv.norm()/dx1Invsave.norm();
      const double dp2 = dx2Inv.norm()/dx2Invsave.norm();
      oops::Log::test() << "Inverse test for block " << it->name() << std::endl;
      ASSERT(dp1 < params.inverseTolerance.value());
      ASSERT(dp2 < params.inverseTolerance.value());
    }
*/
    return 0;
  }

 private:
  std::string appname() const override {
    return "saber::SaberBlockTest<" + MODEL::name() + ">";
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_SABERBLOCKTEST_H_
