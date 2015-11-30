/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliEMCALTriggerBitConfig.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerBitConfig)
ClassImp(AliEMCALTriggerBitConfigOld)
ClassImp(AliEMCALTriggerBitConfigNew)
/// \endcond

/**
 * Dummy constructor for the configuraiton base classes, not to be callled
 */
AliEMCALTriggerBitConfig::AliEMCALTriggerBitConfig():
    TNamed(),
    fL0Bit(-1),
    fJHighBit(-1),
    fJLowBit(-1),
    fGHighBit(-1),
    fGLowBit(-1),
    fBkgBit(-1),
    fTriggerTypesEnd(-1)
{
}

/**
 * Constructor initialising the configurations. Used by the inheriting classes
 * @param l0bit Trigger bit for Level0 triggers
 * @param jhighbit Trigger bit for jet trigger, high threshold
 * @param jlowbit Trigger bit for jet trigger, low threshold
 * @param ghighbit Trigger bit for gamma trigger, high threshold
 * @param glowbit Trigger bit for gamma trigger, low threshold
 * @param mcoffset Offset for MC
 */
AliEMCALTriggerBitConfig::AliEMCALTriggerBitConfig(
    Int_t l0bit,
    Int_t jhighbit,
    Int_t jlowbit,
    Int_t ghighbit,
    Int_t glowbit,
    Int_t bkgbit,
    Int_t mcoffset):
    TNamed("EmcalTriggerBitConfigUninit", ""),
    fL0Bit(l0bit),
    fJHighBit(jhighbit),
    fJLowBit(jlowbit),
    fGHighBit(ghighbit),
    fGLowBit(glowbit),
    fBkgBit(bkgbit),
    fTriggerTypesEnd(mcoffset)
{
}

/**
 * Initialise from other object
 * @param ref Reference to initialize this trigger bit configuaration from
 */
void AliEMCALTriggerBitConfig::Initialise(const AliEMCALTriggerBitConfig& ref) {
  SetName(ref.GetName());
  fL0Bit = ref.GetLevel0Bit();
  fJHighBit = ref.GetJetHighBit();
  fJLowBit = ref.GetJetLowBit();
  fGHighBit = ref.GetGammaHighBit();
  fGLowBit = ref.GetGammaLowBit();
  fBkgBit = ref.GetBkgBit();
  fTriggerTypesEnd = ref.GetTriggerTypesEnd();
}

/**
 * Constructor, initializing the configuration
 */
AliEMCALTriggerBitConfigOld::AliEMCALTriggerBitConfigOld():
    AliEMCALTriggerBitConfig(0,2,2,1,1,2,3)       // To be checked
{
  SetName("EmcalTriggerBitConfigOld");
}

/**
 * Constructor, initializing the configuration
 */
AliEMCALTriggerBitConfigNew::AliEMCALTriggerBitConfigNew():
    AliEMCALTriggerBitConfig(0,3,4,1,2,5,5)       // To be checked
{
  SetName("EmcalTriggerBitConfigNew");
}

