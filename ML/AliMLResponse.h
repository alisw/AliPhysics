#ifndef ALIMLRESPONSE_H
#define ALIMLRESPONSE_H

// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file AliMLResponse.h
/// \brief Implementation of a C++ interface in AliPhysics for the
///        configuration of ML application on grid
/// \author pietro.fecchio@cern.ch, maximiliano.puccio@cern.ch

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "yaml-cpp/yaml.h"

#include "AliExternalBDT.h"

#include "TNamed.h"

using std::map;
using std::pair;
using std::string;
using std::vector;

/////////////////////////////////////////////////////////////////////////////////////////

class ModelHandler {
public:
  ModelHandler() : model(), path(), library(), scorecut() {}
  ModelHandler(const YAML::Node &node)
      : model(), path(node["path"].as<string>()), library(node["library"].as<string>()),
        scorecut(node["cut"].as<double>()) {}

  string const &GetPath() const { return path; }
  string const &GetLibrary() const { return library; }
  double const &GetScoreCut() const { return scorecut; }

  AliExternalBDT &GetModel() { return model; }

  bool CompileModel();

private:
  AliExternalBDT model;

  string path;
  string library;

  double scorecut;
};

/////////////////////////////////////////////////////////////////////////////////////////

class AliMLResponse : public TNamed {
public:
  AliMLResponse();
  AliMLResponse(const Char_t* name, const Char_t* title);
  virtual ~AliMLResponse();

  AliMLResponse(const AliMLResponse &source);
  AliMLResponse &operator=(const AliMLResponse &source);

  /// method to set the object name
  void SetName(const string name) { fMLResponseName = name; }
  /// method to set yaml config file
  void SetConfigFilePath(const string configfilepath) { fConfigFilePath = configfilepath; }

  /// method to check whether the config file is formally correct
  void CheckConfigFile(YAML::Node nodelist);
  /// method to configure the AliMLResponse object from the config file and compile the models usign treelite
  void MLResponseInit();    /// (it has to be done run time)

  /// return the bin index
  int FindBin(double binvar);
  /// return the MLModel predicted score (raw or proba, depending on useraw)
  double Predict(double binvar, map<string, double> varmap);
  /// return true if predicted score for map is above the threshold given in the config
  bool IsSelected(double binvar, map<string, double> varmap);

protected:
  string fMLResponseName;    /// unique name of this ML response object
  string fConfigFilePath;    /// path of the config file

  vector<ModelHandler> fModels;
  vector<int> fCentClasses;         /// centrality classes ([cent_min, cent_max])
  vector<float> fBins;              /// bin edges for the binned variable (pt/ct)
  vector<string> fVariableNames;    /// bin edges for the binned variable (pt/ct)

  int fNBins;         /// number of bins stored for consistency checks
  int fNVariables;    /// number of variables (features) stored for checks

  vector<float>::iterator fBinsBegin;    /// evaluate just once is better

  bool fRaw;

  /// \cond CLASSIMP
  ClassDef(AliMLResponse, 1);    ///
  /// \endcond
};
#endif
