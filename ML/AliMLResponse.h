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

//**************************************************************************************
// \class AliMLResponse
// \brief helper class to handle application of ML models trained with python libraries
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// P. Fecchio, pietro.fecchio@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AliExternalBDT.h"

using std::map;
using std::pair;
using std::string;
using std::vector;

class AliMLResponse : public TObject {
public:
  enum libraries { kXGBoost, kLightGBM, kModelLibrary };

  AliMLResponse();
  AliMLResponse(string configfilename);
  virtual ~AliMLResponse();

  AliMLResponse(const AliMLResponse &source);
  AliMLResponse &operator=(const AliMLResponse &source);

  /// method to initialise and compile models (it has to be done run time)
  // void InitModels();

  /// method to set yaml config file
  string SetConfigFilePath(const string configfilename);
  /// method to set the list of varibles used for predictions
  // void SetVariableList(const vector<string> &variablelist);

  /// methods to get ML response
  //   bool IsSelectedML(double &prob, AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF = nullptr,
  //                     int masshypo = 0);
  //   bool IsSelectedML(double &prob, double pt, vector<double> variables);
  //   double PredictProbaML(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF = nullptr, int masshypo = 0);
  //   double PredictProbaML(double pt, vector<double> variables);

  /// method to get variable (feature) from map
  //   double GetVariable(string name = "") { return fVars[name]; }

protected:
  // string GetFilePath(const string path);
  // int FindModelForCandidate(double pt);
  // int FindModelForCandidate(double pt, pair<int, int>);
  // int FindModelForCandidate(double pt, pair<int, int>, double ct);

  /// method used to define map of name <-> variables (features) --> to be implemented for each derived class
  // virtual void SetMapOfVariables(AliAODRecoDecayHF * /*cand*/, double /*bfield*/, AliAODPidHF * /*pidHF*/,
  //                                int /*masshypo*/) {
  //   return;
  // }

  // map<string, int> kLibMap = {{"kXGBoost", kXGBoost}, {"kLightGBM", kLightGBM}, {"kModelLibrary", kModelLibrary}};
  string fConfigFilePath;         /// path of the config file
  pair <float [2], vector<AliExternalBDT>> fModels; /// vector of ML models
  // vector<string> fModelLibraries;  /// python libraries used to train ML model (kXGBoost, kLightGBM, or
  // kModelLibrary) vector<string> fModelVarNames;   /// vector with names of variables (features) used in ML model (in
  // the correct order) vector<string> fModelPaths;      /// vector of paths of the files containing the ML model
  // vector<double> fModelOutputCuts; /// vector of model output cuts
  // vector<double> fPtBinsModel;     /// vector with pT bins defined for the ML model application
  // vector<pair<int, int>> fCentClassModel; /// vector with cent classes for the ML model application
  // vector<double> fCtBinsModel;            /// vector with pT bins defined for the ML model application
  // map<string, double> fVars;              /// map of variables (features) that can be used for the ML model
  // application

  /// \cond CLASSIMP
  ClassDef(AliMLResponse, 1); ///
                              /// \endcond
};
#endif
