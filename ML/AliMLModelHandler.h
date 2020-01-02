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
/// \brief Implementation of a C++ interface in AliPhysics for the
///        configuration of ML application on grid
/// \author pietro.fecchio@cern.ch, maximiliano.puccio@cern.ch, fabio.catalano@cern.ch

#include <string>

#include "AliExternalBDT.h"

namespace YAML {
  class Node;
}

class AliMLModelHandler {
public:
  enum {kXGBoost, kLightGBM, kModelLibrary};

  AliMLModelHandler();
  AliMLModelHandler(const YAML::Node &node);

  std::string const &GetPath() const { return path; }
  std::string const &GetLibrary() const { return library; }
  double const &GetScoreCut() const { return scorecut; }
  AliExternalBDT &GetModel() { return model; }

  bool CompileModel();

private:
  AliExternalBDT model;

  std::string path;
  std::string library;

  double scorecut;
};

#endif
