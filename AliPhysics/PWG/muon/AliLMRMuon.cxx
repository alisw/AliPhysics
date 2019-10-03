/**************************************************************************
 * Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
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

//========================================================================
//
//     Contact author: boris.teyssier@cern.ch | antonio.uras@cern.ch
//
//=========================================================================

// --- Standard libraries ---
#include <Riostream.h>

#include "TObject.h"

#include "TClonesArray.h"
#include "AliLMRMuon.h"

ClassImp(AliLMRMuon)

const Double_t AliLMRMuon::kMuonMass = 0.105658357;

AliLMRMuon::AliLMRMuon() : TObject(),fTriggerMatch(0),fLocalBoard(0),
				   fSelectionMask(0),fpDCA(0.), fRabs(0.),fCharge(0),
				   fPx(0.),fPy(0.), fPz(0.), fChi2(0.), fChi2Match(0.)
{
	
}

AliLMRMuon::AliLMRMuon(const AliLMRMuon& muon) : TObject(), 
							     fTriggerMatch(muon.fTriggerMatch),
							     fLocalBoard(muon.fLocalBoard),
							     fSelectionMask(muon.fSelectionMask),
							     fpDCA(muon.fpDCA),
							     fRabs(muon.fRabs),
							     fCharge(muon.fCharge),
							     fPx(muon.fPx),
							     fPy(muon.fPy),
							     fPz(muon.fPz),
							     fChi2(muon.fChi2),
							     fChi2Match(muon.fChi2Match)
{
  ///copy constructor
  
}

TLorentzVector AliLMRMuon::P4() const
{
  TLorentzVector v(fPx,fPy,fPz,Energy());
  return v;
}
