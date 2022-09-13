/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/Run.h"
#include "quench/Traits.h"
#include "saber/oops/instantiateSaberBlockFactory.h"
#include "saber/oops/SaberBlockTest.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::instantiateSaberBlockFactory();
  saber::SaberBlockTest<quench::Traits> dir;
  return run.execute(dir);
}
