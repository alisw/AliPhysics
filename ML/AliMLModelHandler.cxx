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
#include "yaml-cpp/yaml.h"

#include <TFile.h>
#include <TGrid.h>
#include <TSystem.h>
#include "AliLog.h"
#include "AliExternalBDT.h"

/// \cond CLASSIMP
ClassImp(AliMLModelHandler);
/// \endcond

//_______________________________________________________________________________
AliMLModelHandler::AliMLModelHandler() : TNamed(), fModel{nullptr}, fPath{}, fLibrary{}, fScoreCut{} {
  //
  // Default constructor
  //
}

//_______________________________________________________________________________
AliMLModelHandler::AliMLModelHandler(const YAML::Node &node)
    : TNamed(), fModel{nullptr}, fPath{node["path"].as<std::string>()},
      fLibrary{node["library"].as<std::string>()}, fScoreCut{node["cut"].as<double>()} {
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
    : TNamed(source.GetName(), source.GetTitle()), fModel{nullptr}, fPath{source.fPath},
      fLibrary{source.fLibrary}, fScoreCut{source.fScoreCut} {
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

  TNamed::operator=(source);

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

  std::map<std::string, int> libraryMap = {{"kXGBoost", AliMLModelHandler::kXGBoost}, 
                                           {"kLightGBM", AliMLModelHandler::kLightGBM},
                                           {"kModelLibrary", AliMLModelHandler::kModelLibrary}};

  std::string localpath = ImportFile(fPath);

  switch (libraryMap[GetLibrary()]) {
    case kXGBoost: {
      return fModel->LoadXGBoostModel(localpath.data());
      break;
    }
    case kLightGBM: {
      return fModel->LoadLightGBMModel(localpath.data());
      break;
    }
    case kModelLibrary: {
      return fModel->LoadModelLibrary(localpath.data());
      break;
    }
    default: {
      return fModel->LoadXGBoostModel(localpath.data());
      break;
    }
  }
}

//_______________________________________________________________________________
std::string AliMLModelHandler::ImportFile(std::string path) {
  std::string modelname = path.substr(path.find_last_of("/") + 1);

  // check if file is in current directory
  if (path.find("/") == std::string::npos) {
    bool checkFile = gSystem->AccessPathName(gSystem->ExpandPathName(path.c_str()));
    if (checkFile) {
      AliFatalClass(Form("Error file %s not found! Exit", path.data()));
    }
    return path;
  }
    
  // check if file is on alien
  if (path.find("alien:") != std::string::npos) {
    if (gGrid == nullptr) {
      TGrid::Connect("alien://");
      if (gGrid == nullptr) {
        AliFatalClass("Connection to GRID not established! Exit");
      }
    }
  }

  std::string newpath = gSystem->pwd() + std::string("/") + modelname.data();
  std::string oldpath = gDirectory->GetPath();

  bool cpStatus = TFile::Cp(path.data(), newpath.data());
  if (!cpStatus) {
    AliFatalClass(Form("Error in coping file %s in the working directory! Exit", path.data()));
  }

  gDirectory->Cd(oldpath.data());

  return newpath;
}