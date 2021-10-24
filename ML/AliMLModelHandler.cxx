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
#include <AliDataFile.h>
#include "AliLog.h"
#include "AliExternalBDT.h"

/// \cond CLASSIMP
ClassImp(AliMLModelHandler);
/// \endcond

//_______________________________________________________________________________
AliMLModelHandler::AliMLModelHandler() : TNamed(), fModel{nullptr}, fPath{}, fLibrary{}, fScoreCut{}, fScoreCutOpt{} {
  //
  // Default constructor
  //
}

//_______________________________________________________________________________
AliMLModelHandler::AliMLModelHandler(const YAML::Node &node)
    : TNamed(), fModel{nullptr}, fPath{node["path"].as<std::string>()},
      fLibrary{node["library"].as<std::string>()}, fScoreCut{}, fScoreCutOpt{} {
  //
  // Standard constructor
  //
  if(node["cut"].IsSequence())
    fScoreCut = node["cut"].as<std::vector<double> >();
  else
    fScoreCut.push_back(node["cut"].as<double>());

  if(node["cut_opt"]) {
    if(node["cut_opt"].IsSequence()) {
      if(fScoreCut.size() != node["cut_opt"].as<std::vector<std::string> >().size())
        AliFatal("Number of score cuts different from number of score cut options! Exit");
      for (const auto &cutOpt: node["cut_opt"].as<std::vector<std::string> >()) {
        if(cutOpt == "lower")
          fScoreCutOpt.push_back(kLowerCut);
        else if(cutOpt == "upper")
          fScoreCutOpt.push_back(kUpperCut);
        else
          AliFatal("Only upper or lower cuts on output scores possible, check your yaml config file! Exit");        
      }
    }
    else {
      if(node["cut_opt"].as<std::string>() == "lower")
        fScoreCutOpt.push_back(kLowerCut);
      else if(node["cut_opt"].as<std::string>() == "upper")
        fScoreCutOpt.push_back(kUpperCut);
      else
        AliFatal("Only upper or lower cuts on output scores possible, check your yaml config file! Exit");        
    }
  }
  else { // if not specified it assumes that is always lower cut
    for (unsigned int iCut=0; iCut<fScoreCut.size(); iCut++)
      fScoreCutOpt.push_back(kLowerCut);
  }

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
      fLibrary{source.fLibrary}, fScoreCut{source.fScoreCut}, fScoreCutOpt{source.fScoreCutOpt} {
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

  fPath        = source.fPath;
  fLibrary     = source.fLibrary;
  fScoreCut    = source.fScoreCut;
  fScoreCutOpt = source.fScoreCutOpt;

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

  // check if file is on cvmfs
  if (path.find("/cvmfs") != std::string::npos) {
    bool checkFile = gSystem->AccessPathName(path.c_str());
    if (checkFile) {
      AliFatalClass(Form("Error file %s not found on CVMFS! Exit", path.data()));
    }
    // only skip copying for certain file types
    std::string fileform = path.substr(path.find_last_of(".")+1,path.size());
    if(fileform == "so" || fileform == "yml" || fileform ==  "yaml"){
      return path;
    }
  }

  // check if file is on cvmfs, but still requires to set the full path (so starts with PWG or AODB)
  if (path.rfind("PWG", 0) != std::string::npos || path.rfind("AODB", 0) != std::string::npos) {
    std::string path_cvmfs = AliDataFile::GetFileName(path);
    if (path_cvmfs == "") {
      AliFatalClass(Form("Error file %s not found on CVMFS! \n  When running locally, "
                         "please export ALICE_DATA=root://eospublic.cern.ch//eos/experiment/alice/analysis-data"
                         "\n  or use a different folder structure in case a path on CVMFS was not intended. Exit", path.data()));
    }
    // only skip copying for certain file types
    std::string fileform = path.substr(path.find_last_of(".")+1,path.size());
    if(path.find("root:") == std::string::npos && (fileform == "so" || fileform == "yml" || fileform ==  "yaml")){
      return path_cvmfs;
    }
    path = path_cvmfs;
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
