#ifndef ALIMLMODELHANDLER_H
#define ALIMLMODELHANDLER_H

// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file AliMLModelHandler.h
/// \brief Utility class to store the compiled model and it's information
/// \author pietro.fecchio@cern.ch, maximiliano.puccio@cern.ch, fabio.catalano@cern.ch

#include <string>

#include "TNamed.h"

namespace YAML {
  class Node;
}
class AliExternalBDT;

class AliMLModelHandler : public TNamed {
public:
  enum {kXGBoost, kLightGBM, kModelLibrary};

  AliMLModelHandler();
  AliMLModelHandler(const YAML::Node &node);
  virtual ~AliMLModelHandler();

  AliMLModelHandler(const AliMLModelHandler &source);
  AliMLModelHandler &operator=(const AliMLModelHandler &source);

  AliExternalBDT *GetModel() { return fModel; }
  std::string const &GetPath() const { return fPath; }
  std::string const &GetLibrary() const { return fLibrary; }
  double const &GetScoreCut() const { return fScoreCut; }

  bool CompileModel();
  static std::string ImportFile(std::string path);

private:
  AliExternalBDT *fModel;  //!<!

  std::string fPath;       ///
  std::string fLibrary;    ///

  double fScoreCut;        ///

/// \cond CLASSIMP
ClassDef(AliMLModelHandler, 1);    ///
/// \endcond
};

#endif
