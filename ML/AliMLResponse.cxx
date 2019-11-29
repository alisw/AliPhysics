// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliMLResponse
// \brief helper class to handle application of ML models trained with python libraries
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// P. Fecchio, pietro.fecchio@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TDirectory.h>
#include <TFile.h>
#include <TGrid.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliMLResponse.h"

#include "assert.h"

namespace {

enum kLibrary { kXGBoost, kLightGBM, kModelLibrary };

map<string, int> kLibraryMap = {{"kXGBoost", kXGBoost}, {"kLightGBM", kLightGBM}, {"kModelLibrary", kModelLibrary}};

string ImportFile(string path) {
  if (path.find("alien:") != string::npos) {
    string modelname = path.substr(path.find_last_of("/") + 1);

    if (gGrid == nullptr) {
      TGrid::Connect("alien://");
      assert(gGrid != nullptr && "Connection to GRID not established! Exit");
    }

    string newpath = gSystem->pwd() + string("/") + modelname.data();
    string oldpath = gDirectory->GetPath();

    bool cpStatus = TFile::Cp(path.data(), newpath.data());
    assert(cpStatus && "Error in coping file from Alien! Exit");

    gDirectory->Cd(oldpath.data());

    return newpath;

  } else {
    return path;
  }
}
}    // namespace

bool ModelHandler::CompileModel() {
  // TODO: try without local copy of the model
  string localpath = ImportFile(path);

  switch (kLibraryMap[GetLibrary()]) {
  case kXGBoost: {
    return model.LoadXGBoostModel(localpath);
    break;
  }
  case kLightGBM: {
    return model.LoadLightGBMModel(localpath);
    break;
  }
  case kModelLibrary: {
    return model.LoadModelLibrary(localpath);
    break;
  }
  default: {
    return model.LoadXGBoostModel(localpath);
    break;
  }
  }
}

/// \cond CLASSIMP
ClassImp(AliMLResponse);
/// \endcond

//________________________________________________________________
AliMLResponse::AliMLResponse()
    : TObject(), fConfigFilePath{}, fModels{}, fCentClasses{}, fBins{}, fVariableNames{}, fNBins{}, fNVariables{},
      fBinsBegin{}, fRaw{} {
  //
  // Default constructor
  //
}

//________________________________________________________________
AliMLResponse::AliMLResponse(string configfilename)
    : TObject(), fConfigFilePath{configfilename}, fModels{}, fCentClasses{}, fBins{}, fVariableNames{}, fNBins{},
      fNVariables{}, fBinsBegin{}, fRaw{} {
  //
  // Standard constructor
  //
}

//________________________________________________________________
AliMLResponse::~AliMLResponse() {
  //
  // Destructor
  //
}

//--------------------------------------------------------------------------
AliMLResponse::AliMLResponse(const AliMLResponse &source)
    : TObject(source), fConfigFilePath{source.fConfigFilePath}, fModels{source.fModels},
      fCentClasses{source.fCentClasses}, fBins{source.fBins}, fVariableNames{source.fVariableNames},
      fNBins{source.fNBins}, fNVariables{source.fNVariables}, fBinsBegin{source.fBinsBegin}, fRaw{source.fRaw} {
  //
  // Copy constructor
  //
}

AliMLResponse &AliMLResponse::operator=(const AliMLResponse &source) {
  //
  // assignment operator
  //
  if (&source == this) return *this;

  TObject::operator=(source);

  fConfigFilePath = source.fConfigFilePath;
  fModels         = source.fModels;
  fCentClasses    = source.fCentClasses;
  fBins           = source.fBins;
  fVariableNames  = source.fVariableNames;
  fNBins          = source.fNBins;
  fNVariables     = source.fNVariables;
  fBinsBegin      = source.fBinsBegin;
  fRaw            = source.fRaw;

  return *this;
}

//_________________________________________________________________________
void AliMLResponse::CheckConfigFile(YAML::Node nodelist) {
  /// error for empty config file
  if (nodelist.IsNull()) {
    AliFatal("Empty .yaml config file, please check it! Exit");
  }
  /// error for bin/model number inconsistencies
  if ((nodelist["BINS"].as<vector<float>>().size() - 1) != nodelist["N_MODELS"].as<float>() ||
      (nodelist["N_MODELS"].as<float>() != nodelist["MODELS"].size())) {
    AliFatal("Inconsistency found in the number of bins/models, please check it! Exit");
  }
  /// error for variables/numberofvariable inconsistency
  if (nodelist["NUM_VAR"].as<float>() != nodelist["VAR_NAMES"].size()) {
    AliFatal("Inconsistency found in the number of varibles, please check it! Exit");
  }
  return;
}

//________________________________________________________________
void AliMLResponse::MLResponseInit() {
  YAML::Node nodeList;
  /// manage wrong config file path
  try {
    nodeList = YAML::LoadFile(fConfigFilePath);
  } catch (std::exception &e) {
    AliFatal(Form("Yaml-ccp error: %s! Exit", e.what()));
  }
  /// manage inconsistencies in config file
  CheckConfigFile(nodeList);

  fVariableNames = nodeList["VAR_NAMES"].as<vector<string>>();
  fBins          = nodeList["BINS"].as<vector<float>>();
  fNBins         = nodeList["N_MODELS"].as<int>();
  fNVariables    = nodeList["NUM_VARS"].as<int>();
  fRaw           = nodeList["RAW_SCORE"].as<bool>();
  fBinsBegin     = fBins.begin();

  for (const auto &model : nodeList["MODELS"]) {
    fModels.push_back(ModelHandler{model});
  }

  for (auto &model : fModels) {
    bool comp = model.CompileModel();
    if (!comp) {
      AliFatal("Error in model compilation! Exit");
    }
  }
}

//________________________________________________________________
int AliMLResponse::FindBin(double binvar) {
  vector<float>::iterator low;
  low = std::lower_bound(fBins.begin(), fBins.end(), binvar);
  return low - fBinsBegin;
}

//________________________________________________________________
double AliMLResponse::Predict(double binvar, map<string, double> varmap, bool useraw) {
  if ((int)varmap.size() >= fNVariables) {
    AliFatal("The variables map you provided to the predictor have a size different from the variable list size! Exit");
  }

  vector<double> features;
  for (const auto &varname : fVariableNames) {
    if (varmap.find(varname) == varmap.end()) {
      AliFatal(Form("Variable |%s| not found in variable list provided in config! Exit", varname.data()));
    }
    features.push_back(varmap[varname]);
  }

  int bin = FindBin(binvar);
  if (bin == 0 || bin == fNBins) {
    AliWarning("Binned variable outside range, no model available! Exit");
    return -999.;
  }

  return fModels[bin - 1].GetModel().Predict(&features[0], fNVariables, useraw);
}

// TODO: far fare la compilazione e testare il tutto con i file presi dal mio alien