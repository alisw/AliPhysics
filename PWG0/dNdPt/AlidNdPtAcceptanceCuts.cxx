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
// last change: 2011-04-04 by M.Knichel

#include <iostream>
#include <TList.h>

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "TParticle.h"

#include "AlidNdPtAcceptanceCuts.h"

using namespace std;

ClassImp(AlidNdPtAcceptanceCuts)

//_____________________________________________________________________________
AlidNdPtAcceptanceCuts::AlidNdPtAcceptanceCuts(const Char_t* name,const Char_t *title) : 
AliAnalysisCuts(name, title)
, fMinEta(0)
, fMaxEta(0)
, fMinPhi(0)
, fMaxPhi(0)
, fMinPt(0)
, fMaxPt(0)
, fExcludeMinEta(0)
, fExcludeMaxEta(0)
, fExcludeMinPhi(0)
, fExcludeMaxPhi(0)
, fCheckRange(kFALSE)
, fMaxDCAr(0)
, fMaxDCAz(0)
{
  // default constructor 
  
  // init data members with defaults
  Init();
}

//_____________________________________________________________________________
AlidNdPtAcceptanceCuts::~AlidNdPtAcceptanceCuts()  
{
  // destructor
}

//_____________________________________________________________________________
void AlidNdPtAcceptanceCuts::Init()  
{
  // set default values
  SetEtaRange();
  SetPhiRange();
  SetPtRange();
  SetMaxDCAr();
  SetMaxDCAz();
}

//_____________________________________________________________________________
Bool_t AlidNdPtAcceptanceCuts::AcceptTrack(AliESDtrack *track)
{
  // check acceptance cuts for AliESDtrack
  if(!track) return kFALSE;

  Float_t eta = track->Eta();
  Float_t phi = track->Phi();
  Float_t pt = track->Pt();

  if(eta < fMinEta) return kFALSE;
  if(eta > fMaxEta) return kFALSE;
  if(phi < fMinPhi) return kFALSE;
  if(phi > fMaxPhi) return kFALSE;
  if(pt < fMinPt) return kFALSE;
  if(pt > fMaxPt) return kFALSE;
  
  if (fCheckRange) {
      if ((eta > fExcludeMinEta) && (eta < fExcludeMaxEta) && (phi > fExcludeMinPhi) && (phi < fExcludeMaxPhi)) { return kFALSE; }
  }

return kTRUE;
}

//_____________________________________________________________________________
Bool_t AlidNdPtAcceptanceCuts::AcceptTrack(AliExternalTrackParam *track)
{
  // check acceptance cuts for AliESDtrack
  if(!track) return kFALSE;

  Float_t eta = track->Eta();
  Float_t phi = track->Phi();
  Float_t pt = track->Pt();

  if(eta < fMinEta) return kFALSE;
  if(eta > fMaxEta) return kFALSE;
  if(phi < fMinPhi) return kFALSE;
  if(phi > fMaxPhi) return kFALSE;
  if(pt < fMinPt) return kFALSE;
  if(pt > fMaxPt) return kFALSE;

return kTRUE;
}

//_____________________________________________________________________________
Bool_t AlidNdPtAcceptanceCuts::AcceptTrack(TParticle *particle)
{
  // check acceptance cuts for TParticle
  if(!particle) return kFALSE;

  Float_t eta = particle->Eta();
  Float_t phi = particle->Phi();
  Float_t pt = particle->Pt();

  if(eta < fMinEta) return kFALSE;
  if(eta > fMaxEta) return kFALSE;
  if(phi < fMinPhi) return kFALSE;
  if(phi > fMaxPhi) return kFALSE;
  if(pt < fMinPt) return kFALSE;
  if(pt > fMaxPt) return kFALSE;

return kTRUE;
}

//_____________________________________________________________________________
Long64_t AlidNdPtAcceptanceCuts::Merge(TCollection* list) 
{
  // Merge list of objects (needed by PROOF)
  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AlidNdPtAcceptanceCuts* entry = dynamic_cast<AlidNdPtAcceptanceCuts*>(obj);
    if (entry == 0)  
      continue; 

  count++;
  }

return count;
}
