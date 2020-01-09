// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file AliMLModelHandler.cxx
/// \author pietro.fecchio@cern.ch, maximiliano.puccio@cern.ch, fabio.catalano@cern.ch

#include "AliMLModelHandler.h"

#include <map>
#include "assert.h"
#include "yaml-cpp/yaml.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TGrid.h>
#include <TSystem.h>
#include "AliExternalBDT.h"

namespace {

std::map<std::string, int> kLibraryMap = {{"kXGBoost", AliMLModelHandler::kXGBoost}, {"kLightGBM", AliMLModelHandler::kLightGBM}, 
                                          {"kModelLibrary", AliMLModelHandler::kModelLibrary}};

std::string ImportFile(std::string path) {
  std::string modelname = path.substr(path.find_last_of("/") + 1);

  if (path.find("alien:") != std::string::npos) {
    if (gGrid == nullptr) {
      TGrid::Connect("alien://");
      assert(gGrid != nullptr && "Connection to GRID not established! Exit");
    }
  }

  std::string newpath = gSystem->pwd() + std::string("/") + modelname.data();
  std::string oldpath = gDirectory->GetPath();

  bool cpStatus = TFile::Cp(path.data(), newpath.data());
  assert(cpStatus && "Error in coping file in the working directory! Exit");

  gDirectory->Cd(oldpath.data());

  return newpath;
}
}    // namespace

//_______________________________________________________________________________
AliMLModelHandler::AliMLModelHandler() :  fModel{nullptr}, fPath{}, fLibrary{}, fScoreCut{} {
  //
  // Default constructor
  //
}

//_______________________________________________________________________________
AliMLModelHandler::AliMLModelHandler(const YAML::Node &node)
    : fModel{nullptr}, fPath{node["path"].as<std::string>()}, fLibrary{node["library"].as<std::string>()},
      fScoreCut{node["cut"].as<double>()} {
  //
  // Standard constructor
  //
  fModel = new AliExternalBDT();
}

AliMLModelHandler::~AliMLModelHandler() {
  //
  // Destructor
  //
  if(fModel)
    delete fModel;
}

//_______________________________________________________________________________
AliMLModelHandler::AliMLModelHandler(const AliMLModelHandler &source)
    : fModel{nullptr}, fPath{source.fPath}, fLibrary{source.fLibrary}, fScoreCut{source.fScoreCut} {
  //
  // Copy constructor
  //
  fModel = new AliExternalBDT(*source.fModel);
}

AliMLModelHandler &AliMLModelHandler::operator=(const AliMLModelHandler &source) {
  //
  // Assignment operator
  //
  if (&source == this) return *this;

  if(fModel)
    delete fModel;
  fModel = new AliExternalBDT(*source.fModel);

  fPath      = source.fPath;
  fLibrary   = source.fLibrary;
  fScoreCut  = source.fScoreCut;

  return *this;
}

//_______________________________________________________________________________
bool AliMLModelHandler::CompileModel() {
  std::string localpath = ImportFile(this->fPath);

  switch (kLibraryMap[GetLibrary()]) {
    case kXGBoost: {
      return this->fModel->LoadXGBoostModel(localpath.data());
      break;
    }
    case kLightGBM: {
      return this->fModel->LoadLightGBMModel(localpath.data());
      break;
    }
    case kModelLibrary: {
      return this->fModel->LoadModelLibrary(localpath.data());
      break;
    }
    default: {
      return this->fModel->LoadXGBoostModel(localpath.data());
      break;
    }
  }
}