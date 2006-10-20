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
#include "AliRsnDaughterCut.h"

ClassImp(AliRsnDaughterCut)

//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnDaughterCut::Pass(AliRsnDaughter* /*track1*/, AliRsnDaughter* /*track2*/) const
{
// 
// Virtual method for cut passing.
// This function must be overridden and return kTRUE when cut is passed.
// 
	TObject::Error("Pass", "This method must be overridden!");
	return kFALSE;
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnDaughterCutPtSingle::Pass(AliRsnDaughter *track1, AliRsnDaughter* /*track2*/) const
{
// 
// Cut on single track momentum.
//
	if (!track1) return kFALSE;
	if (track1->GetPt() < fPtMin) return kFALSE;
	if (track1->GetPt() > fPtMax) return kFALSE;
	
	return kTRUE;
}
//--------------------------------------------------------------------------------------------------------
Bool_t AliRsnDaughterCutPtPair::Pass(AliRsnDaughter *track1, AliRsnDaughter *track2) const
{
// 
// Cut on single track momentum.
// 
	if (!track1 || !track2) return kFALSE;
	
	AliRsnDaughter sum = AliRsnDaughter::Sum(*track1, *track2);
	
	if (sum.GetPt() < fPtMin) return kFALSE;
	if (sum.GetPt() > fPtMax) return kFALSE;
	
	return kTRUE;
}
