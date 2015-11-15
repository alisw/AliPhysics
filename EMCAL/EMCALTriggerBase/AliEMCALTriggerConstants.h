/**
 * @file AliEMCALTriggerConstants.h
 * @date Nov. 15, 2015
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 */
#ifndef ALIEMCALTRIGGERCONSTANTS_H
#define ALIEMCALTRIGGERCONSTANTS_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TString.h>

namespace EMCALTrigger {
/**
 * \enum EMCalTriggerType_t
 * \brief Definition of different trigger patch types
 *
 * This enumeration defines the different trigger patch types
 * processed by the trigger maker. Each trigger patch type has
 * a certain patch size and therefore a certain length and
 * geometric center
 */
enum EMCalTriggerType_t {
  kEMCalL0            = 0,    ///< EMCAL Level0 trigger patches
  kEMCalL1GammaHigh   = 1,    ///< EMCAL Gamma trigger patches, high threshold
  kEMCalL1GammaLow    = 2,    ///< EMCAL Gamma trigger patches, low threshold
  kEMCalL1JetHigh     = 3,    ///< EMCAL Jet trigger patches, high threshold
  kEMCalL1JetLow      = 4,    ///< EMCAL Jet trigger patches, low threshold
  // leave some space for future implementation
  kEMCalRecalcL0      = 15,   ///< EMCAL Level0 patches, recalculated
  kEMCalRecalcL1Gamma = 16,   ///< EMCAL Gamma patches, recalculated
  kEMCalRecalcL1Jet   = 17,   ///< EMCAL Jet patches, recalculated
  kEMCalRecalcL1Bkg   = 18,   ///< EMCAL Background patches, recalculated
};

const TString kEMCalTriggerNames[32] = {
    "EMCL0",
    "EMCGAH",
    "EMCGAL",
    "EMCJEH",
    "EMCJEL",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "EMCREL0",
    "EMCREGA",
    "EMCREJE",
    "EMCREBKG",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    ""
};

}

#endif
