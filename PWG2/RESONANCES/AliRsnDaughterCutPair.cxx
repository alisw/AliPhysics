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

//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnDaughterCutPair::Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const
{
// 
// Virtual method for cut passing.
// This function checks only that the two arguments are not NULL.
// 
	if (!track1 || !track2) return kFALSE;
	
	return kTRUE;
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnDaughterCutPairPt::Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const
{
// 
// Cut on single track momentum.
// 
	if (!AliRsnDaughterCutPair::Pass(track1, track2)) return kFALSE;
	
	AliRsnDaughter sum = AliRsnDaughter::Sum(*track1, *track2);
	
	if (sum.GetPt() < fPtMin) return kFALSE;
	if (sum.GetPt() > fPtMax) return kFALSE;
	
	return kTRUE;
}
//--------------------------------------------------------------------------------------------------------
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
	if (track1->GetSign() > 0) {
		qLpos = track1->GetPx()*sum.GetPx() + track1->GetPy()*sum.GetPy() + track1->GetPz()*sum.GetPz();
		qLneg = track2->GetPx()*sum.GetPx() + track2->GetPy()*sum.GetPy() + track2->GetPz()*sum.GetPz();
		qL1   = qLpos;
	} else {
		qLneg = track1->GetPx()*sum.GetPx() + track1->GetPy()*sum.GetPy() + track1->GetPz()*sum.GetPz();
		qLpos = track2->GetPx()*sum.GetPx() + track2->GetPy()*sum.GetPy() + track2->GetPz()*sum.GetPz();
		qL1   = qLneg;
	}
	qLpos /= sum.GetP();
	qLneg /= sum.GetP();
	
	// compute Qt
	qt = TMath::Sqrt(track1->GetP2() - qL1*qL1);
	
	// compute alpha
	alpha = (qLpos - qLneg) / (qLpos + qLneg);
}
//--------------------------------------------------------------------------------------------------------
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
