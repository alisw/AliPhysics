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
#include <AliEMCALTriggerBitConfig.h>

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerBitConfig)
ClassImp(AliEMCALTriggerBitConfigOld)
ClassImp(AliEMCALTriggerBitConfigNew)
/// \endcond

AliEMCALTriggerBitConfig::AliEMCALTriggerBitConfig():
    TNamed(),
    fL0Bit(-1),
    fJHighBit(-1),
    fJLowBit(-1),
    fGHighBit(-1),
    fGLowBit(-1),
    fTriggerTypesEnd(-1)
{
}

AliEMCALTriggerBitConfig::AliEMCALTriggerBitConfig(
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
}

void AliEMCALTriggerBitConfig::Initialise(const AliEMCALTriggerBitConfig& ref) {
  SetName(ref.GetName());
  fL0Bit = ref.GetLevel0Bit();
  fJHighBit = ref.GetJetHighBit();
  fJLowBit = ref.GetJetLowBit();
  fGHighBit = ref.GetGammaHighBit();
  fGLowBit = ref.GetGammaLowBit();
  fTriggerTypesEnd = ref.GetTriggerTypesEnd();
}

AliEMCALTriggerBitConfigOld::AliEMCALTriggerBitConfigOld():
    AliEMCALTriggerBitConfig(0,2,2,1,1,3)       // To be checked
{
  SetName("EmcalTriggerBitConfigOld");
}

AliEMCALTriggerBitConfigNew::AliEMCALTriggerBitConfigNew():
    AliEMCALTriggerBitConfig(0,3,4,1,2,5)       // To be checked
{
  SetName("EmcalTriggerBitConfigNew");
}

