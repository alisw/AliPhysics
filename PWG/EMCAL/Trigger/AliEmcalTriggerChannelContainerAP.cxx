/**
 * @file AliEmcalTriggerChannelContainerAP.cxx
 * @since
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <AliEmcalTriggerChannelContainerAP.h>

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerChannelContainerAP)
ClassImp(AliEmcalTriggerChannelContainerAP::AliEmcalTriggerChannelPositionAP)
/// \endcond


void AliEmcalTriggerChannelContainerAP::AddChannel(int col, int row){
  if(HasChannel(col, row)) return;
  fChannels.Add(new AliEmcalTriggerChannelPositionAP(col, row));
}


Bool_t AliEmcalTriggerChannelContainerAP::HasChannel(int col, int row){
  AliEmcalTriggerChannelPositionAP refChannel;
  if(fChannels.FindObject(&refChannel)) return true;
  return false;
}

Bool_t AliEmcalTriggerChannelContainerAP::AliEmcalTriggerChannelPositionAP::IsEqual(const TObject *ref){
  const AliEmcalTriggerChannelPositionAP *refpos = dynamic_cast<const AliEmcalTriggerChannelPositionAP *>(ref);
  if(!refpos) return false;
  return fCol == refpos->fCol && fRow == refpos->fRow;
}

Int_t AliEmcalTriggerChannelContainerAP::AliEmcalTriggerChannelPositionAP::Compare(const TObject *ref){
  const AliEmcalTriggerChannelPositionAP *refpos = dynamic_cast<const AliEmcalTriggerChannelPositionAP *>(ref);
  if(!refpos) return 1;
  if(fCol == refpos->fCol){
    if(fRow == refpos->fRow) return 0;
    else if(fRow < refpos->fRow) return -1;
    else return 1;
  }
  else if(fCol < refpos->fCol) return -1;
  else return 1;
}
