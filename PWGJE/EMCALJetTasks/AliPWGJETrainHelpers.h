#ifndef ALIPWGJETRAINHELPERS_H
#define ALIPWGJETRAINHELPERS_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**
 * @class AliPWGJETrainHelpers
 * @brief Helpers for configure the PWGJE trains.
 * @ingroup PWGJEBASE
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
   * // Determine the variables.
   * std::vector<std::string> tempVariables = AliPWGJETrainHelpers::ExtractAliEnProductionValuesForLEGOTrain();
   * // Assign them to their final global variables.
   * const char* kPeriod = tempVariables.at(0).c_str();
   * const char* kColType = tempVariables.at(1).c_str();
   * const bool kMC = (tempVariables.at(2) == "true" ? true : false);
   * const bool kIsRun2 = (tempVariables.at(3) == "true" ? true : false);
   * ~~~
   *
   * Note that we cannot return the values by taking them by reference in the arguments because that would involve
   * reassigning global variables, which is not allowed.
   *
   * The returned parameters correspond to:
   *
   * - period (std::string): The run period. For example, "LHC15o".
   * - collType (std::string): The collision type. Can be "pp", "pPb", "PbPb", etc.
   * - mc (bool): True if this is an MC production.
   * - isRun2 (bool): True if the run period is in Run 2.
   * @return std::vector<std::string> containing {run period, collision type, mc, isRun2}.
   */
  static std::vector<std::string> ExtractAliEnProductionValuesForLEGOTrain();
};

#endif /* AliPWGJETrainHelpers.h */
