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
 * Basic kinematic cuts for single track selection
 *
 *   Author: Markus Fasel
 */
#include <TMath.h>
#include "AliEMCalTriggerKineCuts.h"
#include "AliVParticle.h"

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerKineCuts::AliEMCalTriggerKineCuts():
  TObject(),
  fPtCut(0.1, 1000.),
  fEtaCut(-0.8, 0.8),
  fPhiCut()
{
}

//______________________________________________________________________________
bool AliEMCalTriggerKineCuts::IsSelected(const AliVParticle* const track) const {
  /*
   * Kinematic track selection
   */
  if(!fPtCut.IsInRange(TMath::Abs(track->Pt()))) return false;
  if(!fEtaCut.IsInRange(track->Eta())) return false;
  if(!fPhiCut.IsInRange(track->Phi())) return false;
  return true;
}

} /* namespace EMCalTriggerPtAnalysis */
