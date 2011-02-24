/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
//                Dielectron TrackRotator                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TObjArray.h>
#include <AliAODTrack.h>
#include <AliESDtrack.h>
#include <TRandom3.h>

#include "AliDielectronTrackRotator.h"

ClassImp(AliDielectronTrackRotator)

AliDielectronTrackRotator::AliDielectronTrackRotator() :
  TNamed(),
  fIterations(1),
  fRotationType(kRotateBothRandom),
  fStartAnglePhi(TMath::Pi()),
  fConeAnglePhi(TMath::Pi()/6.),
  fkArrTracksP(0x0),
  fkArrTracksN(0x0),
  fCurrentIteration(0),
  fCurrentTackP(0),
  fCurrentTackN(0),
  fTrackP(0x0),
  fTrackN(0x0)
{
  //
  // Default Constructor
  //
  gRandom->SetSeed();
}

//______________________________________________
AliDielectronTrackRotator::AliDielectronTrackRotator(const char* name, const char* title) :
  TNamed(name, title),
  fIterations(1),
  fRotationType(kRotateBothRandom),
  fStartAnglePhi(TMath::Pi()),
  fConeAnglePhi(TMath::Pi()/6.),
  fkArrTracksP(0x0),
  fkArrTracksN(0x0),
  fCurrentIteration(0),
  fCurrentTackP(0),
  fCurrentTackN(0),
  fTrackP(0x0),
  fTrackN(0x0)
{
  //
  // Named Constructor
  //
  gRandom->SetSeed();
}

//______________________________________________
AliDielectronTrackRotator::~AliDielectronTrackRotator()
{
  //
  // Default Destructor
  //
  
}

//______________________________________________
void AliDielectronTrackRotator::Reset()
{
  //
  // Reset the current iterators
  //
  fCurrentIteration=0;
  fCurrentTackP=0;
  fCurrentTackN=0;
}

//______________________________________________
Bool_t AliDielectronTrackRotator::NextCombination()
{
  //
  // Perform track rotation of the tracks in the track arrays as long as there are possible combinations
  //
  if (!fkArrTracksP || !fkArrTracksP) {
    Reset();
    return kFALSE;
  }

  Int_t nP=fkArrTracksP->GetEntriesFast();
  Int_t nN=fkArrTracksN->GetEntriesFast();
  if (nP==0||nN==0){
    Reset();
    return kFALSE;
  }
  
  if (fCurrentIteration==fIterations){
    fCurrentIteration=0;
    ++fCurrentTackP;
  }
  
  if (fCurrentTackP==nP){
    ++fCurrentTackN;
    fCurrentTackP=0;
  }
  
  if (fCurrentTackN==nN){
    Reset();
    return kFALSE;
  }
  
  if (!RotateTracks()){
    Reset();
    return kFALSE;
  }
  
  ++fCurrentIteration;
  return kTRUE;
}

//______________________________________________
Bool_t AliDielectronTrackRotator::RotateTracks()
{
  //
  // Actual track rotation
  // Find out particle type and perform the rotation
  //

  const AliVTrack *trackP=dynamic_cast<AliVTrack*>(fkArrTracksP->UncheckedAt(fCurrentTackP));
  const AliVTrack *trackN=dynamic_cast<AliVTrack*>(fkArrTracksN->UncheckedAt(fCurrentTackN));
  if (!trackP||!trackN) return kFALSE;

  
  Double_t angle  = fStartAnglePhi+(2*gRandom->Rndm()-1)*fConeAnglePhi;
  Int_t    charge = TMath::Nint(gRandom->Rndm());
  
  if (trackP->IsA()==AliESDtrack::Class()) {
    
    if (!fTrackP) {
      fTrackP=new AliESDtrack;
      fTrackN=new AliESDtrack;
    }
    
    trackP->Copy(*fTrackP);
    trackN->Copy(*fTrackN);
    
    if (fRotationType==kRotatePositive||(fRotationType==kRotateBothRandom&&charge==0)){
      ((AliESDtrack*)fTrackP)->Rotate(angle);
    }
    
    if (fRotationType==kRotateNegative||(fRotationType==kRotateBothRandom&&charge==1)){
      ((AliESDtrack*)fTrackN)->Rotate(angle);
    }
    
  } else if (trackP->IsA()==AliAODTrack::Class()) {
    
    if (!fTrackP) {
      fTrackP=new AliAODTrack;
      fTrackN=new AliAODTrack;
    }
    
    (*(AliAODTrack*)fTrackP)=(*(AliAODTrack*)trackP);
    (*(AliAODTrack*)fTrackN)=(*(AliAODTrack*)trackN);
        
    if (fRotationType==kRotatePositive||(fRotationType==kRotateBothRandom&&charge==0)){
      Double_t phi=fTrackP->Phi()+angle;
      if (phi>2*TMath::Pi()) phi-=2*TMath::Pi();
      ((AliAODTrack*)fTrackP)->SetPhi(phi);
    }
    
    if (fRotationType==kRotateNegative||(fRotationType==kRotateBothRandom&&charge==1)){
      Double_t phi=fTrackN->Phi()+angle;
      if (phi>2*TMath::Pi()) phi-=2*TMath::Pi();
      ((AliAODTrack*)fTrackN)->SetPhi(phi);
    }
    
  }

  return kTRUE;
}
