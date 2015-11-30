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

const Double_t kEMCL1ADCtoGeV = 0.07874;             ///< Conversion from EMCAL Level1 ADC to energy


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
  kTMEMCalLevel0 = 0,     ///< EMCAL Level0 patches
  kTMEMCalGammaL = 1,     ///< EMCAL Gamma trigger
  kTMEMCalGammaH = 2,     ///< EMCAL Gamma trigger
  kTMEMCalJetL   = 3,     ///< EMCAL Jet trigger
  kTMEMCalJetH   = 4,     ///< EMCAL Jet trigger
  kTMEMCalBkg    = 5      ///< EMCAL background
};

const TString kEMCalTriggerNames[6] = {
    "EMCL0",
    "EMCGAL",
    "EMCGAH",
    "EMCJEL",
    "EMCJEH",
    "EMCBKG"
};

}

#endif
