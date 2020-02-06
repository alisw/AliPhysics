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

#include <map>
#include <string>
#include <vector>

#include "TNamed.h"

#include "AliMLModelHandler.h"

namespace YAML {
  class Node;
}

class AliMLResponse : public TNamed {
public:
  AliMLResponse();
  AliMLResponse(const Char_t *name, const Char_t *title);
  virtual ~AliMLResponse();

  AliMLResponse(const AliMLResponse &source);
  AliMLResponse &operator=(const AliMLResponse &source);

  /// method to set yaml config file
  void SetConfigFilePath(const std::string configfilepath) { fConfigFilePath = configfilepath; }

  /// method to check whether the config file is formally correct
  void CheckConfigFile(YAML::Node nodelist);
  /// method to configure the AliMLResponse object from the config file and compile the models usign treelite
  void MLResponseInit();    /// (it has to be done run time)

  /// return the bin index
  int FindBin(double binvar);
  /// return the ML model predicted score (raw or proba, depending on useraw)
  double Predict(double binvar, std::map<std::string, double> varmap);
  /// overload to pass directly a vector of variables
  double Predict(double binvar, std::vector<double> variables);
  /// return true if predicted score for map is above the threshold given in the config
  bool IsSelected(double binvar, std::map<std::string, double> varmap);
  /// overload for getting the model score too
  template<typename F>
  bool IsSelected(double binvar, std::map<std::string, double> varmap, F &score);
  /// overload to pass directly a vector of variables
  bool IsSelected(double binvar, std::vector<double> variables);
  /// overload for getting the model score too
  template<typename F>
  bool IsSelected(double binvar, std::vector<double> variables, F &score);

protected:
  std::string fConfigFilePath;              /// path of the config file

  std::vector<AliMLModelHandler> fModels;   //!<! vector of models
  std::vector<int> fCentClasses;            /// centrality classes ([cent_min, cent_max])
  std::vector<float> fBins;                 /// bin edges for the binned variable (pt/ct)
  std::vector<std::string> fVariableNames;  /// bin edges for the binned variable (pt/ct)

  int fNBins;                               /// number of bins stored for consistency checks
  int fNVariables;                          /// number of variables (features) stored for checks

  std::vector<float>::iterator fBinsBegin;  //!<!  evaluate just once is better

  bool fRaw;                                /// set to true to use raw score instead of probability

  /// \cond CLASSIMP
  ClassDef(AliMLResponse, 2);    ///
  /// \endcond
};

template<typename F>
bool AliMLResponse::IsSelected(double binvar, std::map<std::string, double> varmap, F &score) {
  int bin = FindBin(binvar);
  score = Predict(binvar, varmap);
  return score >= fModels[bin - 1].GetScoreCut();
}

template<typename F>
bool AliMLResponse::IsSelected(double binvar, std::vector<double> variables, F &score) {
  int bin = FindBin(binvar);
  score = Predict(binvar, variables);
  return score >= fModels[bin - 1].GetScoreCut();
}

#endif
