/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliEMCALTriggerChannelContainer.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerChannelContainer)
ClassImp(AliEMCALTriggerChannelContainer::AliEMCALTriggerChannelPosition)
/// \endcond

void AliEMCALTriggerChannelContainer::AddChannel(int col, int row){
  if(HasChannel(col, row)) return;
  fChannels.Add(new AliEMCALTriggerChannelPosition(col, row));
}

Bool_t AliEMCALTriggerChannelContainer::HasChannel(int col, int row){
  AliEMCALTriggerChannelPosition refChannel;
  if(fChannels.FindObject(&refChannel)) return true;
  return false;
}


Bool_t AliEMCALTriggerChannelContainer::AliEMCALTriggerChannelPosition::IsEqual(const TObject *ref) const {
  const AliEMCALTriggerChannelPosition *refpos = dynamic_cast<const AliEMCALTriggerChannelPosition *>(ref);
  if(!refpos) return false;
  return fCol == refpos->fCol && fRow == refpos->fRow;
}

Int_t AliEMCALTriggerChannelContainer::AliEMCALTriggerChannelPosition::Compare(const TObject *ref) const {
  const AliEMCALTriggerChannelPosition *refpos = dynamic_cast<const AliEMCALTriggerChannelPosition *>(ref);
  if(!refpos) return 1;
  if(fCol == refpos->fCol){
    if(fRow == refpos->fRow) return 0;
    else if(fRow < refpos->fRow) return -1;
    else return 1;
  }
  else if(fCol < refpos->fCol) return -1;
  else return 1;
}
