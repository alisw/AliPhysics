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
#include "AliEMCalTriggerEventSelection.h"
#include "AliEMCalTriggerEventData.h"
#include <TString.h>
#include "AliAnalysisUtils.h"
#include "AliVEvent.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerEventSelection)

using namespace PWGJE::EMCALJetTasks;

/**
 * Main Constructor
 */
AliEMCalTriggerEventSelection::AliEMCalTriggerEventSelection():
  TObject(),
  fVertexCut(-10., 10.),
  fOldPileupSelection(kFALSE),
  fOldVertexSelection(kFALSE)
{

}

/**
 * Apply basic event selection
 *
 * Can be overwritten by inheriting classes
 *
 * @param ev Combined event container
 * @return event selection decision (true if event is selected)
 */
bool AliEMCalTriggerEventSelection::IsEventSelected(const AliEMCalTriggerEventData* const ev) const {
  AliAnalysisUtils evutils;
  AliVEvent *recEvent = ev->GetRecEvent();
  if(!fOldVertexSelection && !evutils.IsVertexSelected2013pA(recEvent)) return kFALSE;
  if(fOldVertexSelection && !FalseVertexSelectionPA2013(recEvent)) return kFALSE;
  if(!fOldPileupSelection && evutils.IsPileUpEvent(recEvent)) return kFALSE;
  if(fOldPileupSelection && recEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return kFALSE;
  if(!fVertexCut.IsInRange(recEvent->GetPrimaryVertex()->GetZ())) return kFALSE;
  return true;
}

/**
 * Do vertex selection in the old buggy way
 * @param ev Event to check
 * @return True if the vertex is selected, false otherwise
 */
bool AliEMCalTriggerEventSelection::FalseVertexSelectionPA2013(const AliVEvent *const ev) const{
  const AliVVertex *trkVtx = ev->GetPrimaryVertex();
  if(!trkVtx || trkVtx->GetNContributors() < 1) return kFALSE;
  if (!TString(trkVtx->GetTitle()).Contains("VertexerTracks")) return kFALSE;

  Float_t zvtx = trkVtx->GetZ();
  const AliVVertex* spdVtx = ev->GetPrimaryVertexSPD();

  if (spdVtx->GetNContributors()< 1) return kFALSE;
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (spdVtx->IsFromVertexerZ() && (zRes>0.25)) return kFALSE; // doing this incorrectly on purpose
  if(TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return kFALSE;

  if (TMath::Abs(zvtx) > 10) return kFALSE;
  return kTRUE;
}
