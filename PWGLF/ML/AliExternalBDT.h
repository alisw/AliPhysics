// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file AliExternalBDT.h
/// \brief Implementation of a C++ interface in AliPhysics for XGboost,
////       LightGBM and SciKit BDT models.
/// \author maximiliano.puccio@cern.ch, pietro.fecchio@cern.ch

#ifndef ALIEXTERNALBDT_H
#define ALIEXTERNALBDT_H

#include "treelite/c_api.h"
#include "treelite/c_api_runtime.h"
#include <string>
#include <vector>

class AliExternalBDT {
public:
  AliExternalBDT(std::string name = "");
  virtual ~AliExternalBDT(){};

  bool LoadLightGBMModel(std::string path);
  bool LoadModelLibrary(std::string path);
  bool LoadXGBoostModel(std::string path);

  double Predict(double *features, int size, bool useRaw = false);

private:
  bool CompileAndLoadModelLibrary();
  bool CreateModelCode();
  std::string GetUniquePath();
  bool LoadModel(const std::string &path, int type);

  std::string fBDTname;       /// Unique name of this external BDT handler
  ModelHandle fModel;
  std::string fModelPath;
  std::string fModelName;
  CompilerHandle fCompiler;
  PredictorHandle fPredictor;
};

#endif
