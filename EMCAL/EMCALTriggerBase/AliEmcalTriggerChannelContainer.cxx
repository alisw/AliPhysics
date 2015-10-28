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
#include "AliEmcalTriggerChannelContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerChannelContainer)
ClassImp(AliEmcalTriggerChannelContainer::AliEmcalTriggerChannelPosition)
/// \endcond

/**
 * Add a new channel with the postion in column and row to the container, In case the channel
 * is already listed in the trigger channel container we don't add it again.
 * \param col Column of the channel
 * \param row Row of the channel
 */
void AliEmcalTriggerChannelContainer::AddChannel(int col, int row){
  if(HasChannel(col, row)) return;
  fChannels.Add(new AliEmcalTriggerChannelPosition(col, row));
}

/**
 * Check whether channel with the position (col, row) is listed in the trigger channel container
 * \param col Column of the channel
 * \param row Row of the channel
 * \return True if the channel is listed, false otherwise
 */
Bool_t AliEmcalTriggerChannelContainer::HasChannel(int col, int row){
  AliEmcalTriggerChannelPosition refChannel;
  if(fChannels.FindObject(&refChannel)) return true;
  return false;
}

/**
 * Check if the object is equal to object ref. Object can only be equal if ref is of
 * the same type (AliEmcalTrigger channel position). If this is the case, col and row
 * of the two objects have to match.
 * \param ref The object to check
 * \return True if objects are equal, false otherwise
 */
Bool_t AliEmcalTriggerChannelContainer::AliEmcalTriggerChannelPosition::IsEqual(const TObject *ref){
  const AliEmcalTriggerChannelPosition *refpos = dynamic_cast<const AliEmcalTriggerChannelPosition *>(ref);
  if(!refpos) return false;
  return fCol == refpos->fCol && fRow == refpos->fRow;
}

/**
 * Compare objects. If objects differ, return always greater (+1). Otherwise compare col and
 * row of the object. Col has priority with respect to row.
 * \param ref The object ot comparte to
 * \return 0 if objects are equal, -1 if this object is smaller, +1 if this object is larger.
 */
Int_t AliEmcalTriggerChannelContainer::AliEmcalTriggerChannelPosition::Compare(const TObject *ref){
  const AliEmcalTriggerChannelPosition *refpos = dynamic_cast<const AliEmcalTriggerChannelPosition *>(ref);
  if(!refpos) return 1;
  if(fCol == refpos->fCol){
    if(fRow == refpos->fRow) return 0;
    else if(fRow < refpos->fRow) return -1;
    else return 1;
  }
  else if(fCol < refpos->fCol) return -1;
  else return 1;
}
