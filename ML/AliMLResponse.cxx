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

#include "yaml-cpp/yaml.h"

#include "AliLog.h"
#include "AliMLResponse.h"

/// \cond CLASSIMP
ClassImp(AliMLResponse);
/// \endcond

//________________________________________________________________
AliMLResponse::AliMLResponse() : TObject(), fConfigFilePath{}, fModel{}, fCentClasses{}, fPtBins{}, fCtBins{} {
  //
  // Default constructor
  //
}

//________________________________________________________________
AliMLResponse::AliMLResponse(string configfilename)
    : TObject(), fConfigFilePath{configfilename}, fModel{}, fCentClasses{}, fPtBins{}, fCtBins{} {
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
    : TObject(source), fConfigFilePath{source.fConfigFilePath}, fModel{source.fModel},
      fCentClasses{source.fCentClasses}, fPtBins{source.fPtBins}, fCtBins{source.fCtBins} {
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
  fModel          = source.fModel;
  fCentClasses    = source.fCentClasses;
  fPtBins         = source.fPtBins;
  fCtBins         = source.fCtBins;
  // fModelLibraries  = source.fModelLibraries;
  // fModelVarNames   = source.fModelVarNames;
  // fModelPaths      = source.fModelPaths;
  // fModelOutputCuts = source.fModelOutputCuts;
  // fPtBinsModel     = source.fPtBinsModel;
  // fPtBinCand       = source.fPtBinCand;
  // fVars            = source.fVars;

  return *this;
}

//_________________________________________________________________________
void AliMLResponse::ImportConfigFile() {
  // if (path.find("alien:") != string::npos) {
  string modelName = fConfigFilePath.substr(fConfigFilePath.find_last_of("/") + 1);

  if (gGrid == nullptr) {
    TGrid::Connect("alien://");
    if (gGrid == nullptr) {
      AliFatal("Connection to GRID not established! Exit");
    }
  }

  string newPath    = gSystem->pwd() + string("/") + modelName.data();
  string oldRootDir = gDirectory->GetPath();

  bool cpStatus = TFile::Cp(fConfigFilePath.data(), newPath.data());
  gDirectory->Cd(oldRootDir.data());

  if (!cpStatus) {
    AliFatal("Error in coping file from Alien! Exit");
  }
  // } else {
  //   return path;
  // }
}

//________________________________________________________________
void AliMLResponse::Config() {
  YAML::Node nodeList;
  /// manage wrong config file path and empty config file
  try {
    nodeList = YAML::LoadFile(fConfigFilePath);
  } catch (std::exception &e) {
    AliFatal(Form("Yaml-ccp error: %s! Exit", e.what());
  }
  /// manage empty config file
  if (nodeList.IsNull()) AliFatal("Empty .yaml config file! Exit");

  // ora devo gestire i 2 casi di configurazioni
  // creo una funzione che genera le configurazioni base (no cent)
  // e poi la chiamo in un loop sulle centralità

  // fModelVarNames   =nodeList["VarNames"].as<vector<string>>();
  // fModelPaths      =nodeList["ModelNames"].as<vector<string>>();
  // fModelOutputCuts =nodeList["ModelOutputCuts"].as<vector<double>>();
  // fPtBinsModel     =nodeList["PtBins"].as<vector<double>>();
  // fModelLibraries  =nodeList["ModelLibraries"].as<vector<string>>();

  // // for consistency check
  // unsigned int numModels =nodeList["NumModels"].as<unsigned int>();
  // unsigned int numVars   =nodeList["NumVars"].as<unsigned int>();

  // if (numModels != fModelPaths.size() || numModels != fModelLibraries.size() || numModels != fModelOutputCuts.size()
  // ||
  //     numModels != fPtBinsModel.size() - 1)
  //   AliFatal("Inconsistency found in the number of models loaded from your yaml config file, please check it! Exit");

  // if (numVars != fModelVarNames.size())
  //   AliFatal("Inconsistency found in the number of variables (features) loaded from your yaml config file, please "
  //            "check it! Exit");
}

//_________________________________________________________________________
// void AliMLResponse::InitModels() {
//   for (auto iMod = 0; iMod < fModelPaths.size(); iMod++) {
//     string localpath     = GetFile(fModelPaths[iMod]);
//     AliExternalBDT model = AliExternalBDT();

//     bool loadmodel = false;

//     switch (kLibMap[fModelLibraries[iMod]]) {
//     case kXGBoost: {
//       loadmodel = model.LoadXGBoostModel(localpath);
//       break;
//     }
//     case kLightGBM: {
//       loadmodel = model.LoadLightGBMModel(localpath);
//       break;
//     }
//     case kModelLibrary: {
//       loadmodel = model.LoadModelLibrary(localpath);
//       break;
//     }
//     default: {
//       loadmodel = model.LoadXGBoostModel(localpath);
//       break;
//     }
//     }
//     if (!loadmodel) AliFatal("Problem in loading model");
//     fModels.push_back(model);
//   }
// }

//_________________________________________________________________________
// int AliMLResponse::FindPtBin(double pt) {
//   int bin = TMath::BinarySearch(fPtBinsModel.size() - 1, &fPtBinsModel[0], pt);
//   if (bin < 0) // underflow --> equal to min value
//     bin = 0;

//   return bin;
// }

//________________________________________________________________
// double AliMLResponse::PredictProbaML(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int masshypo) {
//   fPtBinCand = FindPtBin(cand->Pt());
//   SetMapOfVariables(cand, bfield, pidHF, masshypo);
//   if (fVars.empty()) {
//     AliWarning("Map of features empty!");
//     return -999.;
//   }

//   vector<double> vecOfFeatures;
//   for (auto varname : fModelVarNames) {
//     if (fVars.find(varname) == fVars.end()) { // variable name not found in variables map!
//       AliWarning("Variable (feature) used for ML training not implemented for model application!");
//       return -999.;
//     }

//     vecOfFeatures.push_back(fVars[varname]);
//   }

//   if (fPtBinCand > fModels.size() - 1) {
//     AliWarning(Form("Model for pT bin %d not loaded!", fPtBinCand));
//     return -999.;
//   }

//   double modelPred = fModels[fPtBinCand].Predict(&vecOfFeatures[0], vecOfFeatures.size());

//   return modelPred;
// }

//________________________________________________________________
// double AliMLResponse::PredictProbaML(double pt, vector<double> variables) {
//   fPtBinCand = FindPtBin(pt);
//   if (fPtBinCand > fModels.size() - 1) {
//     AliWarning(Form("Model for pT bin %d not loaded!", fPtBinCand));
//     return -999.;
//   }

//   double modelPred = fModels[fPtBinCand].Predict(&variables[0], variables.size());

//   return modelPred;
// }

//________________________________________________________________
// bool AliMLResponse::IsSelectedML(double &prob, AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF,
//                                  int masshypo) {
//   prob = PredictProbaML(cand, bfield, pidHF, masshypo);
//   if (prob > fModelOutputCuts[fPtBinCand]) return true;

//   return false;
// }

//________________________________________________________________
// bool AliMLResponse::IsSelectedML(double &prob, double pt, vector<double> variables) {
//   prob = PredictProbaML(pt, variables);
//   if (prob > fModelOutputCuts[fPtBinCand]) return true;

//   return false;
// }
