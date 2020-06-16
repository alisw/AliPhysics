// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file AliMLResponse.cxx
/// \author pietro.fecchio@cern.ch, maximiliano.puccio@cern.ch

#include "AliMLResponse.h"

#include "yaml-cpp/yaml.h"

#include "AliExternalBDT.h"
#include "AliLog.h"

using std::map;
using std::pair;
using std::string;
using std::vector;

/// \cond CLASSIMP
ClassImp(AliMLResponse);
/// \endcond

//_______________________________________________________________________________
AliMLResponse::AliMLResponse()
    : TNamed(), fConfigFilePath{}, fModels{}, fCentClasses{}, fBins{}, fVariableNames{}, fNBins{}, fNVariables{},
      fBinsBegin{}, fRaw{} {
  //
  // Default constructor
  //
}

//_______________________________________________________________________________
AliMLResponse::AliMLResponse(const Char_t *name, const Char_t *title)
    : TNamed(name, title), fConfigFilePath{""}, fModels{}, fCentClasses{}, fBins{}, fVariableNames{}, fNBins{},
      fNVariables{}, fBinsBegin{}, fRaw{} {
  //
  // Standard constructor
  //
}

//_______________________________________________________________________________
AliMLResponse::~AliMLResponse() {
  //
  // Destructor
  //
}

//_______________________________________________________________________________
AliMLResponse::AliMLResponse(const AliMLResponse &source)
    : TNamed(source.GetName(), source.GetTitle()), fConfigFilePath{source.fConfigFilePath}, fModels{source.fModels},
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

  TNamed::operator=(source);

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

//_______________________________________________________________________________
void AliMLResponse::CheckConfigFile(YAML::Node nodelist) {
  /// error for empty config file
  if (nodelist.IsNull()) {
    AliFatal("Empty .yaml config file, please check it! Exit");
  }
  /// error for bin/model number inconsistencies
  if ((nodelist["BINS"].as<vector<float>>().size() - 1) != nodelist["N_MODELS"].as<unsigned int>() ||
      (nodelist["N_MODELS"].as<unsigned int>() != nodelist["MODELS"].size())) {
    AliFatal("Inconsistency found in the number of bins/models, please check it! Exit");
  }
  /// error for variables/numberofvariable inconsistency
  if (nodelist["NUM_VAR"].as<unsigned int>() != nodelist["VAR_NAMES"].size()) {
    AliFatal("Inconsistency found in the number of variables, please check it! Exit");
  }
  return;
}

//_______________________________________________________________________________
void AliMLResponse::CompileModels(std::string configLocalPath) {
  YAML::Node nodeList;
  /// manage wrong config file path
  try {
    nodeList = YAML::LoadFile(configLocalPath);
  } catch (std::exception &e) {
    AliFatal(Form("Yaml-ccp error: %s! Exit", e.what()));
  }
  /// manage inconsistencies in config file
  CheckConfigFile(nodeList);

  fVariableNames = nodeList["VAR_NAMES"].as<vector<string>>();
  fBins          = nodeList["BINS"].as<vector<float>>();
  fNBins         = nodeList["N_MODELS"].as<int>();
  fNVariables    = nodeList["NUM_VAR"].as<int>();
  fRaw           = nodeList["RAW_SCORE"].as<bool>();

  fBinsBegin = fBins.begin();

  for (const auto &model : nodeList["MODELS"]) {
    fModels.push_back(AliMLModelHandler{model});
  }

  for (auto &model : fModels) {
    bool comp = model.CompileModel();
    if (!comp) {
      AliFatal("Error in model compilation! Exit");
    }
  }
}

//_______________________________________________________________________________
void AliMLResponse::MLResponseInit() {
  /// import config file from alien path
  string configLocalPath = ImportConfigFile();
  CompileModels(configLocalPath);
}

//_______________________________________________________________________________
int AliMLResponse::FindBin(double binvar) {
  vector<float>::iterator low = std::lower_bound(fBins.begin(), fBins.end(), binvar);
  int bin = low - fBinsBegin;
  if (bin == 0 || bin >= fNBins) {
    AliWarning("Binned variable outside range, no model available!");
    return -1;
  }
  return bin;
}

//_______________________________________________________________________________
double AliMLResponse::Predict(double binvar, map<string, double> varmap) {
  if ((int)varmap.size() < fNVariables) {
    AliFatal("The variable map you provided to the predictor has a size smaller than the variable list size! Exit");
  }

  vector<double> features;
  for (const auto &varname : fVariableNames) {
    if (varmap.find(varname) == varmap.end()) {
      AliFatal(Form("Variable |%s| not found in variable list provided in config! Exit", varname.data()));
    }
    features.push_back(varmap[varname]);
  }

  int bin = FindBin(binvar);
  if (bin < 0)
    return -999.;

  return fModels.at(bin - 1).GetModel()->Predict(&features[0], fNVariables, fRaw);
}

//_______________________________________________________________________________
double AliMLResponse::Predict(double binvar, vector<double> variables) {
  if ((int)variables.size() != fNVariables) {
    AliFatal(Form("Number of variables passed (%d) different from the one used in the model (%d)! Exit",
                  (int)variables.size(), fNVariables));
  }

  int bin = FindBin(binvar);
  if (bin < 0)
    return -999.;

  return fModels.at(bin - 1).GetModel()->Predict(&variables[0], fNVariables, fRaw);
}

//_______________________________________________________________________________
bool AliMLResponse::IsSelected(double binvar, std::map<std::string, double> varmap) {
  double score{0.};
  return IsSelected(binvar, varmap, score);
}

//_______________________________________________________________________________
bool AliMLResponse::IsSelected(double binvar, std::vector<double> variables) {
  double score{0.};
  return IsSelected(binvar, variables, score);
}