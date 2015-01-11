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

/*
 * Trigger bit configuration used in the trigger patch maker and by the trigger patches
 * themselves in order to identify of which type the trigger patch is. Can be adapted to different
 * trigger bit configurations use in different reconstructions
 *
 *   Author: Markus Fasel
 */

#include "AliEmcalTriggerBitConfig.h"

ClassImp(AliEmcalTriggerBitConfig)
ClassImp(AliEmcalTriggerBitConfigOld)
ClassImp(AliEmcalTriggerBitConfigNew)

//________________________________________________________________________
AliEmcalTriggerBitConfig::AliEmcalTriggerBitConfig():
    TNamed(),
    fL0Bit(-1),
    fJHighBit(-1),
    fJLowBit(-1),
    fGHighBit(-1),
    fGLowBit(-1),
    fTriggerTypesEnd(-1)
{
  /*
   * Dummy constructor for the configuraiton base classes, not to be callled
   */
}

//________________________________________________________________________
AliEmcalTriggerBitConfig::AliEmcalTriggerBitConfig(
    Int_t l0bit,
    Int_t jhighbit,
    Int_t jlowbit,
    Int_t ghighbit,
    Int_t glowbit,
    Int_t mcoffset):
    TNamed("EmcalTriggerBitConfigUninit", ""),
    fL0Bit(l0bit),
    fJHighBit(jhighbit),
    fJLowBit(jlowbit),
    fGHighBit(ghighbit),
    fGLowBit(glowbit),
    fTriggerTypesEnd(mcoffset)
{
  /*
   * Constructor initialising the configurations. Used by the inheriting classes
   */
}

//________________________________________________________________________
void AliEmcalTriggerBitConfig::Initialise(const AliEmcalTriggerBitConfig& ref) {
  /*
   * Initialise from other object
   */
  SetName(ref.GetName());
  fL0Bit = ref.GetLevel0Bit();
  fJHighBit = ref.GetJetHighBit();
  fJLowBit = ref.GetJetLowBit();
  fGHighBit = ref.GetGammaHighBit();
  fGLowBit = ref.GetGammaLowBit();
  fTriggerTypesEnd = ref.GetTriggerTypesEnd();
}

//________________________________________________________________________
AliEmcalTriggerBitConfigOld::AliEmcalTriggerBitConfigOld():
    AliEmcalTriggerBitConfig(0,2,2,1,1,3)       // To be checked
{
  /*
   * Settings for the 2-bit configuration
   */
  SetName("EmcalTriggerBitConfigOld");
}

//________________________________________________________________________
AliEmcalTriggerBitConfigNew::AliEmcalTriggerBitConfigNew():
    AliEmcalTriggerBitConfig(0,3,4,1,2,5)       // To be checked
{
  /*
   * Settings for the 4-bit configuration
   */
  SetName("EmcalTriggerBitConfigNew");
}

