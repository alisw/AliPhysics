/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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
 * Basic event selection component: Selects events according to the pA cut and a vertex-z cut
 * For more sophisticated event selection the method IsEventSelected has to be overwritten
 *
 *   Author: Markus Fasel
 */
#include "AliEMCalTriggerEventSelection.h"
#include "AliEMCalTriggerEventData.h"
#include <TString.h>
#include "AliAnalysisUtils.h"
#include "AliVEvent.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerEventSelection)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerEventSelection::AliEMCalTriggerEventSelection():
  TObject(),
  fVertexCut(-10., 10.)
{
  /*
   * Main Constructor
   */

}

//______________________________________________________________________________
bool AliEMCalTriggerEventSelection::IsEventSelected(const AliEMCalTriggerEventData* const ev) const {
  /*
   * Apply basic event selection
   *
   * Can be overwritten by inheriting classes
   *
   * @param ev: Combined event container
   * @return: event selection decision (true if event is selected)
   */
  AliAnalysisUtils evutils;
  AliVEvent *recEvent = ev->GetRecEvent();
  if(!evutils.IsVertexSelected2013pA(recEvent)) return kFALSE;
  if(evutils.IsPileUpEvent(recEvent)) return kFALSE;
  if(!fVertexCut.IsInRange(recEvent->GetPrimaryVertex()->GetZ())) return kFALSE;
  return true;
}

} /* namespace EMCalTriggerPtAnalysis */
