#ifndef ALIPWGJETRAINHELPERS_H
#define ALIPWGJETRAINHELPERS_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
 * @class AliPWGJETrainHelpers
 * @brief Helpers for configure the PWGJE trains.
 *
 * Set of helper functions to aid in configuration of PWGJE trains.
 *
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University 
 * @date Nov 09, 2018
 */

class AliPWGJETrainHelpers
{
 public:
  /**
   * Automatically extract AliEn production parameters from environment variables for use with the LEGO train.
   * These parameters will be used to configure the shared global variables on each train.
   *
   * Usage should look like:
   *
   * ~~~{.cxx}
   * std::string period = "";
   * std::string collType = "";
   * bool kMC = false;
   * bool kIsRun2 = false;
   * AliPWGJETrainHelpers::ExtractAliEnProductionValuesForLEGOTrain(period, collType, kMC, kIsRun2);
   * // Set the final string variables, which are expected to be c strings.
   * const char* kPeriod = period.c_str();
   * const char* kColType = collType.c_str();
   * ~~~
   *
   * @param[out] period The run period. For example, "LHC15o".
   * @param[out] collType The collision type. Can be "pp", "pPb", "PbPb", etc.
   * @param[out] mc True if this is an MC production.
   * @param[out] isRun2 True if the run period is in Run 2.
   */
  static void ExtractAliEnProductionValuesForLEGOTrain(std::string& period, std::string& collType,
                             bool& mc, bool& isRun2);
};

#endif /* AliPWGJETrainHelpers.h */
