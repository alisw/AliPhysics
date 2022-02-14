#ifndef ALIEMCALTRIGGERANAHELPER_H
#define ALIEMCALTRIGGERANAHELPER_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

namespace PWGJE{ 

namespace EMCALJetTasks {

/**
 * \enum ETATriggerType
 * \brief Trigger types defined for this analysis
 *
 * This enumeration declared switches for then trigger types defined for
 * and used in this analysis.
 */
enum ETATriggerType{
  kTAEMCJHigh       = 0,            ///< Jet trigger, high threshold
  kTAEMCJLow        = 1,            ///< Jet trigger, low threshold
  kTAEMCGHigh       = 2,            ///< Gamma trigger, high threshold
  kTAEMCGLow        = 3,            ///< Gamma trigger, low threshold
  kTAUndef          = -1            ///< Default, undefined
};

/**
 * \enum EPatchEnergyType
 * Definition of the method to obtain the energy for an additional energy cut in the
 * selection of trigger patches.
 */
enum EPatchEnergyType_t{
  kAmplitudeOnline,                 ///< Online amplitude of the patch, from L0 time sums
  kAmplitudeOffline,                ///< Offline amplitude, estimated from EMCAL cells
  kEnergyOnline,                    ///< Online energy, estimated from L0 time sums
  kEnergyOffline                    ///< Offline energy, from cells, calibrated, exluding hot towers
};

/**
 * \enum ETriggerMethod_t
 * \brief Methods available to select event as triggered events
 *
 * This enumeration defines the possible methods to select events
 * as triggered events using the information stored in the trigger
 * decision object.
 */
enum ETriggerMethod_t{
  kTriggerString,                   ///< kTriggerString
  kTriggerPatches,                  ///< kPatches
  kTriggerMixed                     ///< kMixed
};

}

}

#endif /* ALIEMCALTRIGGERANAHELPER_H */
