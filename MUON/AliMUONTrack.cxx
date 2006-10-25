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

/* $Id$ */

///////////////////////////////////////////////////
//
// Reconstructed track
// in
// ALICE
// dimuon
// spectrometer
//
///////////////////////////////////////////////////

#include "AliMUONTrack.h"

#include "AliMUONTrackParam.h" 
#include "AliMUONHitForRec.h" 
#include "AliMUONSegment.h" 
#include "AliMUONTriggerTrack.h"
#include "AliMUONConstants.h"

#include "AliLog.h"

#include <Riostream.h> // for cout
#include <TMath.h>
#include <TMatrixD.h>
#include <TObjArray.h>

#include <stdlib.h> // for exit()

ClassImp(AliMUONTrack) // Class implementation in ROOT context

//__________________________________________________________________________
AliMUONTrack::AliMUONTrack()
  : TObject(),
    fTrackParamAtVertex(),
    fTrackParamAtHit(0x0),
    fHitForRecAtHit(0x0),
    fNTrackHits(0),
    fFitFMin(-1.),
    fMatchTrigger(kFALSE),
    fChi2MatchTrigger(0.),
    fTrackID(0)
{
  // Default constructor
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack(AliMUONSegment* BegSegment, AliMUONSegment* EndSegment)
  : TObject(),
    fTrackParamAtVertex(),
    fTrackParamAtHit(0x0),
    fHitForRecAtHit(0x0),
    fNTrackHits(0),
    fFitFMin(-1.),
    fMatchTrigger(kFALSE),
    fChi2MatchTrigger(0.),
    fTrackID(0)
{
  // Constructor from two Segment's

  fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",10);
  fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",10);
  
  if (BegSegment) { //AZ
    AddTrackParamAtHit(0,BegSegment->GetHitForRec1());
    AddTrackParamAtHit(0,BegSegment->GetHitForRec2());
    AddTrackParamAtHit(0,EndSegment->GetHitForRec1());
    AddTrackParamAtHit(0,EndSegment->GetHitForRec2());
    fTrackParamAtHit->Sort(); // sort TrackParamAtHit according to increasing Z
  }
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack(AliMUONSegment* Segment, AliMUONHitForRec* HitForRec)
  : TObject(),
    fTrackParamAtVertex(),
    fTrackParamAtHit(0x0),
    fHitForRecAtHit(0x0),
    fNTrackHits(0),
    fFitFMin(-1.),
    fMatchTrigger(kFALSE),
    fChi2MatchTrigger(0.),
    fTrackID(0)
{
  // Constructor from one Segment and one HitForRec

  fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",10);
  fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",10);
  
  AddTrackParamAtHit(0,Segment->GetHitForRec1());
  AddTrackParamAtHit(0,Segment->GetHitForRec2());
  AddTrackParamAtHit(0,HitForRec);
  fTrackParamAtHit->Sort(); // sort TrackParamAtHit according to increasing Z
}

  //__________________________________________________________________________
AliMUONTrack::~AliMUONTrack()
{
  // Destructor
  if (fTrackParamAtHit) {
    // delete the TClonesArray of pointers to TrackParam
    delete fTrackParamAtHit;
    fTrackParamAtHit = NULL;
  }

  if (fHitForRecAtHit) {
    // delete the TClonesArray of pointers to HitForRec
    delete fHitForRecAtHit;
    fHitForRecAtHit = NULL;
  }
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack (const AliMUONTrack& theMUONTrack)
  : TObject(theMUONTrack),
    fTrackParamAtVertex(theMUONTrack.fTrackParamAtVertex),
    fTrackParamAtHit(0x0),
    fHitForRecAtHit(0x0),
    fNTrackHits(theMUONTrack.fNTrackHits),
    fFitFMin(theMUONTrack.fFitFMin),
    fMatchTrigger(theMUONTrack.fMatchTrigger),
    fChi2MatchTrigger(theMUONTrack.fChi2MatchTrigger),
    fTrackID(theMUONTrack.fTrackID)
{
  
  Int_t maxIndex = 0;
  
  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (theMUONTrack.fTrackParamAtHit) {
    maxIndex = (theMUONTrack.fTrackParamAtHit)->GetEntriesFast();
    fTrackParamAtHit  =  new TClonesArray("AliMUONTrackParam",maxIndex);
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fTrackParamAtHit)[index]) AliMUONTrackParam(*(AliMUONTrackParam*)theMUONTrack.fTrackParamAtHit->At(index));
    }
  }  

  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (theMUONTrack.fHitForRecAtHit) {
    maxIndex = (theMUONTrack.fHitForRecAtHit)->GetEntriesFast();
    fHitForRecAtHit  =  new TClonesArray("AliMUONHitForRec",maxIndex);
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fHitForRecAtHit)[index]) AliMUONHitForRec(*(AliMUONHitForRec*)theMUONTrack.fHitForRecAtHit->At(index));
    }
  }  

}

  //__________________________________________________________________________
AliMUONTrack & AliMUONTrack::operator=(const AliMUONTrack& theMUONTrack)
{

  // check assignement to self
  if (this == &theMUONTrack)
    return *this;

  // base class assignement
  TObject::operator=(theMUONTrack);

  fTrackParamAtVertex = theMUONTrack.fTrackParamAtVertex;

  Int_t maxIndex = 0;
  
  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  fTrackParamAtHit = 0;
  if (theMUONTrack.fTrackParamAtHit) {
    fTrackParamAtHit  =  new TClonesArray("AliMUONTrackParam",10);
    maxIndex = (theMUONTrack.fTrackParamAtHit)->GetEntriesFast();
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fTrackParamAtHit)[fTrackParamAtHit->GetEntriesFast()])
      	AliMUONTrackParam(*(AliMUONTrackParam*)(theMUONTrack.fTrackParamAtHit)->At(index));
    }
  }  

  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  fHitForRecAtHit = 0;
  if (theMUONTrack.fHitForRecAtHit) {
    fHitForRecAtHit  =  new TClonesArray("AliMUONHitForRec",10);
    maxIndex = (theMUONTrack.fHitForRecAtHit)->GetEntriesFast();
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fHitForRecAtHit)[fHitForRecAtHit->GetEntriesFast()])
      	AliMUONHitForRec(*(AliMUONHitForRec*)(theMUONTrack.fHitForRecAtHit)->At(index));
    }
  }  
  
  fNTrackHits         =  theMUONTrack.fNTrackHits;
  fFitFMin            =  theMUONTrack.fFitFMin;
  fMatchTrigger       =  theMUONTrack.fMatchTrigger;
  fChi2MatchTrigger   =  theMUONTrack.fChi2MatchTrigger;
  fTrackID            =  theMUONTrack.fTrackID;

  return *this;
}

  //__________________________________________________________________________
void AliMUONTrack::AddTrackParamAtHit(AliMUONTrackParam *trackParam, AliMUONHitForRec *hitForRec) 
{
  // Add TrackParamAtHit if "trackParam" != NULL else create empty TrackParamAtHit
  // Update link to HitForRec if "hitForRec" != NULL
  if (!fTrackParamAtHit) {
    fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",10);  
    fNTrackHits = 0;
  }
  AliMUONTrackParam* TrackParamAtHit;
  if (trackParam) TrackParamAtHit = new ((*fTrackParamAtHit)[fNTrackHits]) AliMUONTrackParam(*trackParam);
  else TrackParamAtHit = new ((*fTrackParamAtHit)[fNTrackHits]) AliMUONTrackParam();
  if (hitForRec) TrackParamAtHit->SetHitForRecPtr(hitForRec);
  fNTrackHits++;
}

  //__________________________________________________________________________
void AliMUONTrack::AddHitForRecAtHit(const AliMUONHitForRec *hitForRec) 
{
  if (!fHitForRecAtHit)
    fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",10); 
  
  if (!hitForRec)
    AliFatal("AliMUONTrack::AddHitForRecAtHit: hitForRec == NULL");
  
  new ((*fHitForRecAtHit)[fHitForRecAtHit->GetEntriesFast()]) AliMUONHitForRec(*hitForRec);
}

  //__________________________________________________________________________
Bool_t* AliMUONTrack::CompatibleTrack(AliMUONTrack * Track, Double_t Sigma2Cut) const
{
  // Return kTRUE/kFALSE for each chamber if hit is compatible or not 
  TClonesArray *hitArray, *thisHitArray;
  AliMUONHitForRec *hit, *thisHit;
  Int_t chamberNumber;
  Float_t deltaZ;
  Float_t deltaZMax = 1.; // 1 cm
  Float_t chi2 = 0;
  Bool_t *nCompHit = new Bool_t[AliMUONConstants::NTrackingCh()]; 

  for ( Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
    nCompHit[ch] = kFALSE;
  }

  thisHitArray = this->GetHitForRecAtHit();

  hitArray =  Track->GetHitForRecAtHit();

  for (Int_t iHthis = 0; iHthis < thisHitArray->GetEntriesFast(); iHthis++) {
    thisHit = (AliMUONHitForRec*) thisHitArray->At(iHthis);
    chamberNumber = thisHit->GetChamberNumber();
    if (chamberNumber < 0 || chamberNumber > AliMUONConstants::NTrackingCh()) continue; 
    nCompHit[chamberNumber] = kFALSE;
    for (Int_t iH = 0; iH < hitArray->GetEntriesFast(); iH++) {
      hit = (AliMUONHitForRec*) hitArray->At(iH);
      deltaZ = TMath::Abs(thisHit->GetZ() - hit->GetZ());
      chi2 = thisHit->NormalizedChi2WithHitForRec(hit,Sigma2Cut); // set cut to 4 sigmas
      if (chi2 < 3. && deltaZ < deltaZMax) {
	nCompHit[chamberNumber] = kTRUE;
	break;
      }
    }  
  }
  
  return nCompHit;
}

  //__________________________________________________________________________
Int_t AliMUONTrack::HitsInCommon(AliMUONTrack* Track) const
{
  // Returns the number of hits in common
  // between the current track ("this")
  // and the track pointed to by "Track".
  Int_t hitsInCommon = 0;
  AliMUONTrackParam *trackParamAtHit1, *trackParamAtHit2;
  // Loop over hits of first track
  trackParamAtHit1 = (AliMUONTrackParam*) this->fTrackParamAtHit->First();
  while (trackParamAtHit1) {
    // Loop over hits of second track
    trackParamAtHit2 = (AliMUONTrackParam*) Track->fTrackParamAtHit->First();
    while (trackParamAtHit2) {
      // Increment "hitsInCommon" if both TrackParamAtHits point to the same HitForRec
      if ((trackParamAtHit1->GetHitForRecPtr()) == (trackParamAtHit2->GetHitForRecPtr())) hitsInCommon++;
      trackParamAtHit2 = (AliMUONTrackParam*) Track->fTrackParamAtHit->After(trackParamAtHit2);
    } // trackParamAtHit2
    trackParamAtHit1 = (AliMUONTrackParam*) this->fTrackParamAtHit->After(trackParamAtHit1);
  } // trackParamAtHit1
  return hitsInCommon;
}

  //__________________________________________________________________________
void AliMUONTrack::RecursiveDump(void) const
{
  // Recursive dump of AliMUONTrack, i.e. with dump of TrackParamAtHit's and attached HitForRec's
  AliMUONTrackParam *trackParamAtHit;
  AliMUONHitForRec *hitForRec;
  cout << "Recursive dump of Track: " << this << endl;
  // Track
  this->Dump();
  for (Int_t trackHitIndex = 0; trackHitIndex < fNTrackHits; trackHitIndex++) {
    trackParamAtHit = (AliMUONTrackParam*) ((*fTrackParamAtHit)[trackHitIndex]);
    // TrackHit
    cout << "TrackParamAtHit: " << trackParamAtHit << " (index: " << trackHitIndex << ")" << endl;
    trackParamAtHit->Dump();
    hitForRec = trackParamAtHit->GetHitForRecPtr();
    // HitForRec
    cout << "HitForRec: " << hitForRec << endl;
    hitForRec->Dump();
  }
  return;
}
  
//_____________________________________________-
void AliMUONTrack::Print(Option_t* opt) const
{
//
  // Printing Track information 
  // "full" option for printing all the information about the track
  //
  TString sopt(opt);
  sopt.ToUpper();
 
  if ( sopt.Contains("FULL") ) { 
    cout << "<AliMUONTrack> No.Clusters=" << setw(2)   << GetNTrackHits() << 
      //      ", Bending P="<< setw(8) << setprecision(5)      << 1./GetInverseBendingMomentum() << 
      //", NonBendSlope=" << setw(8) << setprecision(5)  << GetNonBendingSlope()*180./TMath::Pi() <<
      //", BendSlope=" << setw(8) << setprecision(5)     << GetBendingSlope()*180./TMath::Pi() <<
      ", Match2Trig=" << setw(1) << GetMatchTrigger()  << 
      ", Chi2-tracking-trigger=" << setw(8) << setprecision(5) <<  GetChi2MatchTrigger() << endl ;
    GetTrackParamAtHit()->First()->Print("full");
  }
  else {
    cout << "<AliMUONTrack>";
    GetTrackParamAtHit()->First()->Print("");

  }
    
}
