/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
 
//-------------------------------------------------------------------------
//                      Class AliRsnDaughterCut
//                     -------------------------
//           Implementation of track cuts for analysis.
//           These cuts must be added to the AliRsnAnalysis
//           object in order to make cut on single particles
//           or on pairs of particles.
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#include <Riostream.h>

#include "AliRsnDaughter.h"
#include "AliRsnDaughterCutPair.h"

ClassImp(AliRsnDaughterCutPair)
//
//--------------------------------------------------------------------------------------------------------
//
Bool_t AliRsnDaughterCutPair::Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const
{
// 
// Virtual method for cut passing.
// This function checks only that the two arguments are not NULL.
// 
   if (!track1 || !track2) return kFALSE;
	
   return kTRUE;
}
//
//--------------------------------------------------------------------------------------------------------
//
Bool_t AliRsnDaughterCutPairPt::Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const
{
// 
// Cut on single track momentum.
// 
   if (!AliRsnDaughterCutPair::Pass(track1, track2)) return kFALSE;
	
   AliRsnDaughter sum = AliRsnDaughter::Sum(*track1, *track2);
	
   if (sum.Pt() < fPtMin) return kFALSE;
   if (sum.Pt() > fPtMax) return kFALSE;
	
   return kTRUE;
}
//
//--------------------------------------------------------------------------------------------------------
//
Bool_t AliRsnDaughterCutPairAngle::Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const
{
// 
// Cut on single track momentum.
// 
	if (!AliRsnDaughterCutPair::Pass(track1, track2)) return kFALSE;
	
	AliRsnDaughter sum = AliRsnDaughter::Sum(*track1, *track2);
    
    Double_t num1 = sum.Px()*track1->Px() + sum.Py()*track1->Py() + sum.Pz()*track1->Pz();
    Double_t den1 = sum.P() * track1->P();
    Double_t cos1 = num1 / den1;
    Double_t num2 = sum.Px()*track2->Px() + sum.Py()*track2->Py() + sum.Pz()*track2->Pz();
    Double_t den2 = sum.P() * track2->P();
    Double_t cos2 = num2 / den2;
    if (cos1 < fAngleMin || cos2 < fAngleMin) return kFALSE;
    if (cos1 > fAngleMax || cos2 > fAngleMax) return kFALSE;
	
	return kTRUE;
}
//
//--------------------------------------------------------------------------------------------------------
//
void AliRsnDaughterCutPairArmenteros::Compute(AliRsnDaughter *track1, AliRsnDaughter *track2, Double_t &qt, Double_t &alpha) const
{
//
// Compute variables for Armenteros plot.
// Put as separate method to allow some monitoring
//
   
   if (!AliRsnDaughterCutPair::Pass(track1, track2)) return;
	
   AliRsnDaughter sum = AliRsnDaughter::Sum(*track1, *track2);
	
   // compute projection of both track momenta along mother's one
   Double_t qLpos = 0.0, qLneg = 0.0, qL1 = 0.0;
   if (track1->Charge() > 0) {
      qLpos = track1->Px()*sum.Px() + track1->Py()*sum.Py() + track1->Pz()*sum.Pz();
      qLneg = track2->Px()*sum.Px() + track2->Py()*sum.Py() + track2->Pz()*sum.Pz();
      qL1   = qLpos;
   } else {
      qLneg = track1->Px()*sum.Px() + track1->Py()*sum.Py() + track1->Pz()*sum.Pz();
      qLpos = track2->Px()*sum.Px() + track2->Py()*sum.Py() + track2->Pz()*sum.Pz();
      qL1   = qLneg;
   }
   qLpos /= sum.P();
   qLneg /= sum.P();
	
   // compute Qt
   qt = TMath::Sqrt(track1->P2() - qL1*qL1);
	
   // compute alpha
   alpha = (qLpos - qLneg) / (qLpos + qLneg);
}
//
//--------------------------------------------------------------------------------------------------------
//
Bool_t AliRsnDaughterCutPairArmenteros::Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const
{
// 
// Check Armenteros variables
//
   if (!AliRsnDaughterCutPair::Pass(track1, track2)) return kFALSE;
	
   Double_t qt, alpha;
   Compute(track1, track2, qt, alpha);
	
   if (qt < fQtMin || qt > fQtMax) return kFALSE;
   if (alpha < fAlphaMin || alpha > fAlphaMax) return kFALSE;
   return kTRUE;
}
