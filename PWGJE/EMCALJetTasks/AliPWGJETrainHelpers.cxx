/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TString.h>
#include <TError.h>

#include "AliPWGJETrainHelpers.h"

std::vector<std::string> AliPWGJETrainHelpers::ExtractAliEnProductionValuesForLEGOTrain()
{
  // NOTE: `Getenv(...)` returns `nullptr` if the environment variable isn't available.
  //       std::string _cannot_ be set to `nullptr`, so we retrieve the result as a `const char*`
  //       and then convert it to a `std::string` if it is not null. We use std::string for
  //       convenience and to simplify string comparison.

  // Automatically set shared common variables
  // The run period (ex. "LHC15o")
  const char * tempPeriod = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTAG");
  std::string period = tempPeriod ? tempPeriod : "";
  // Either "pp", "pPb" or "PbPb"
  const char * tempCollType = gSystem->Getenv("ALIEN_JDL_LPMINTERACTIONTYPE");
  std::string collType = tempCollType ? tempCollType : "";
  // Used to get MC.
  const char * tempProdType = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTYPE");
  std::string prodType = tempProdType ? tempProdType : "";

  // In the case of the metadataset, the standard environment variables above will not be set.
  // Instead, we attempt to retrieve the variable from the first enabled metadataset. This requires looping
  // over the available metadatasets until we find a non-null value.
  // We assume that all children in the metadataset should have the same collision system and production type.
  // Although this isn't the case for the period, since the variables are supposed to be the same
  // in all metadatasets (according to Markus Z.), the track cuts presumably must all be the same for
  // all periods in a metadataset.
  const char * tempChildDatasets = gSystem->Getenv("CHILD_DATASETS");
  // Used to get MC.
  const std::string childDatasets = tempChildDatasets ? tempChildDatasets : "";
  const int nChildDatasets = std::atoi(childDatasets.c_str());
  for (int i = 1; i <= nChildDatasets; i++) {
    if (period == "") {
      tempPeriod = gSystem->Getenv(TString::Format("ALIEN_JDL_child_%i_LPMPRODUCTIONTAG", i));
      period = tempPeriod ? tempPeriod : "";
    }
    if (collType == "") {
      tempCollType = gSystem->Getenv(TString::Format("ALIEN_JDL_child_%i_LPMINTERACTIONTYPE", i));
      collType = tempCollType ? tempCollType : "";
    }
    if (prodType == "") {
      tempProdType = gSystem->Getenv(TString::Format("ALIEN_JDL_child_%i_LPMPRODUCTIONTYPE", i));
      prodType = tempProdType ? tempProdType : "";
    }
  }
  // Validate the extract variables. Each one must be set at this point.
  if (period == "" || collType == "" || prodType == "") {
    // Somehow failed to extract the variables which should always be available.
    ::Fatal("AliPWGJETrainHelpers", "Somehow failed to extract the period, collision type, or production type.\n");
  }

  // Determine if it's an MC production.
  const bool mc = (prodType == "MC");

  // Determine if it's Run 2
  std::string yearString = period.substr(3, 2);
  int productionYear = std::atoi(yearString.c_str());
  const bool isRun2 = productionYear > 14;

  // We need to return multiple values, so we use a vector and convert the bools to strings temporarily.
  // They should be converted back after returning.
  // (Not that the bool text doesn't really matter as long as it's non-zero, but we select "true" as it
  // corresponds to the variables being true).
  return std::vector<std::string>{period, collType, mc ? "true" : "", isRun2 ? "true" : ""};
}

