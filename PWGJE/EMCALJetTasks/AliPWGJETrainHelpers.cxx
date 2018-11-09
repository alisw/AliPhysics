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

#include <iostream>
#include <string>
#include <cstdlib>

#include <TSystem.h>
#include <TString.h>
#include <TError.h>

#include "AliPWGJETrainHelpers.h"

void AliPWGJETrainHelpers::ExtractAliEnProductionValuesForLEGOTrain(std::string& period, std::string& collType,
                                  bool& mc, bool& isRun2)
{
  // Automatically set shared common variables
  // The run period (ex. "LHC15o")
  period = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTAG");
  // either "pp", "pPb" or "PbPb"
  collType = gSystem->Getenv("ALIEN_JDL_LPMINTERACTIONTYPE");
  // Used to get mc
  std::string prodType = gSystem->Getenv("ALIEN_JDL_LPMPRODUCTIONTYPE");

  // In the case of the metadataset, the standard environment variables above will not be set.
  // Instead, we attempt to retrieve the variable from the first enabled metadataset. This requires looping
  // over the available metadatasets untils we find a non-null value.
  // We assume that all children in the metadataset should have the same collision system and production type.
  // Although this isn't the case for the period, since the variables are supposed to be the same
  // in all metadatasets (according to Markus Z.), the track cuts presumably must all be the same for
  // all periods in a metadataset.
  const std::string childDatasets = gSystem->Getenv("CHILD_DATASETS"); // Used to get mc
  const int nChildDatasets = std::atoi(childDatasets.c_str());
  for (int i = 1; i <= nChildDatasets; i++) {
    if (period == "") {
      period = gSystem->Getenv(TString::Format("ALIEN_JDL_child_%i_LPMPRODUCTIONTAG", i));
    }
    if (collType == "") {
      collType = gSystem->Getenv(TString::Format("ALIEN_JDL_child_%i_LPMINTERACTIONTYPE", i));
    }
    if (prodType == "") {
      prodType = gSystem->Getenv(TString::Format("ALIEN_JDL_child_%i_LPMPRODUCTIONTYPE", i));
    }
  }
  // Validate the extract variables. Each one must be set at this point.
  if (period == "" || collType == "" || prodType == "") {
    // Somehow failed to extract the vaariables which should always be available.
    ::Fatal("AliPWGJETrainHelpers", "Somehow failed to extract the period, collision type, or production type.\n");
  }

  // Determine if it's an MC production.
  mc = (prodType == "MC");

  // Determine if it's Run2
  std::string yearString = period.substr(3, 2);
  int productionYear = std::atoi(yearString.c_str());
  isRun2 = productionYear > 14;

  // Values are returned via the function arguments, so there is nothing else to do here.
}
