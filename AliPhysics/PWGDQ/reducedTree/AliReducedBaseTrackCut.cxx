/*
***********************************************************
  Implementation of AliReducedBaseTrackCut class.
  Contact: iarsene@cern.ch
  2015/09/10
  *********************************************************
*/

#ifndef ALIREDUCEDBASETRACKCUT_H
#include "AliReducedBaseTrackCut.h"
#endif

#include "AliReducedBaseTrack.h"

ClassImp(AliReducedBaseTrackCut)

//____________________________________________________________________________
AliReducedBaseTrackCut::AliReducedBaseTrackCut() :
  AliReducedInfoCut(),
  fPtRange(),
  fCutOnPt(kFALSE),
  fEtaRange(),
  fCutOnEta(kFALSE),
  fPhiRange(),
  fCutOnPhi(kFALSE)
{
  //
  // default constructor
  //
}

//____________________________________________________________________________
AliReducedBaseTrackCut::AliReducedBaseTrackCut(const Char_t* name, const Char_t* title) :
  AliReducedInfoCut(name, title),
  fPtRange(),
  fCutOnPt(kFALSE),
  fEtaRange(),
  fCutOnEta(kFALSE),
  fPhiRange(),
  fCutOnPhi(kFALSE)
{
  //
  // named constructor
  //
}

//____________________________________________________________________________
AliReducedBaseTrackCut::~AliReducedBaseTrackCut() {
  //
  // destructor
  //
}

//____________________________________________________________________________
Bool_t AliReducedBaseTrackCut::IsSelected(TObject* obj) {
  //
  // apply cuts
  //
  if(!obj->InheritsFrom(AliReducedBaseTrack::Class())) return kFALSE;
  
  AliReducedBaseTrack* track = (AliReducedBaseTrack*)obj;
  
  if(fCutOnPt && (track->Pt()<fPtRange[0] || track->Pt()>fPtRange[1])) return kFALSE;
  if(fCutOnEta && (track->Eta()<fEtaRange[0] || track->Eta()>fEtaRange[1])) return kFALSE;
  if(fCutOnPhi && (track->Phi()<fPhiRange[0] || track->Phi()>fPhiRange[1])) return kFALSE;
  
  return kTRUE;
}
