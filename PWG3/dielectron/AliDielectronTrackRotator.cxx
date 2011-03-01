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
#include <TRandom3.h>
#include <TObjArray.h>

#include <AliVTrack.h>

#include "AliDielectronHelper.h"

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
  fEvent(0x0),
  fTrackP(),
  fTrackN(),
  fVTrackP(0x0),
  fVTrackN(0x0),
  fPdgLeg1(-11),
  fPdgLeg2(11)
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
  fEvent(0x0),
  fTrackP(),
  fTrackN(),
  fVTrackP(0x0),
  fVTrackN(0x0),
  fPdgLeg1(-11),
  fPdgLeg2(11)
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

  AliVTrack *trackP=dynamic_cast<AliVTrack*>(fkArrTracksP->UncheckedAt(fCurrentTackP));
  AliVTrack *trackN=dynamic_cast<AliVTrack*>(fkArrTracksN->UncheckedAt(fCurrentTackN));
  fTrackP.Initialize();
  fTrackN.Initialize();
  fVTrackP=0x0;
  fVTrackN=0x0;
  if (!trackP||!trackN) return kFALSE;
  fTrackP+=AliKFParticle(*trackP,fPdgLeg1);
  fTrackN+=AliKFParticle(*trackN,fPdgLeg2);

  fVTrackP=trackP;
  fVTrackN=trackN;
  
  Double_t angle  = fStartAnglePhi+(2*gRandom->Rndm()-1)*fConeAnglePhi;
  Int_t    charge = TMath::Nint(gRandom->Rndm());
  
  if (fRotationType==kRotatePositive||(fRotationType==kRotateBothRandom&&charge==0)){
    AliDielectronHelper::RotateKFParticle(&fTrackP, angle, fEvent);
  }
  
  if (fRotationType==kRotateNegative||(fRotationType==kRotateBothRandom&&charge==1)){
    AliDielectronHelper::RotateKFParticle(&fTrackN, angle, fEvent);
  }

  return kTRUE;
}
