/*
***********************************************************
  Implementation of AliReducedTrackCut class.
  Contact: iarsene@cern.ch
  2016/09/07
  *********************************************************
*/

#ifndef ALIREDUCEDTRACKCUT_H
#include "AliReducedTrackCut.h"
#endif

#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedVarManager.h"

ClassImp(AliReducedTrackCut)

//____________________________________________________________________________
AliReducedTrackCut::AliReducedTrackCut() :
  AliReducedVarCut(),
  fRejectKinks(kFALSE),
  fRejectTaggedGamma(kFALSE),
  fRejectTaggedPureGamma(kFALSE),
  fRequestITSrefit(kFALSE),
  fCutOnITShitMap(0),
  fUseANDonITShitMap(kFALSE),
  fRequestCutOnITShitMap(kFALSE),  
  fRequestTPCrefit(kFALSE),
  fRequestTOFout(kFALSE)
{
  //
  // default constructor
  //
}

//____________________________________________________________________________
AliReducedTrackCut::AliReducedTrackCut(const Char_t* name, const Char_t* title) :
  AliReducedVarCut(name, title),
  fRejectKinks(kFALSE),
  fRejectTaggedGamma(kFALSE),
  fRejectTaggedPureGamma(kFALSE),
  fRequestITSrefit(kFALSE),
  fCutOnITShitMap(0),
  fUseANDonITShitMap(kFALSE),
  fRequestCutOnITShitMap(kFALSE),  
  fRequestTPCrefit(kFALSE),
  fRequestTOFout(kFALSE)
{
  //
  // named constructor
  //
}

//____________________________________________________________________________
AliReducedTrackCut::~AliReducedTrackCut() {
  //
  // destructor
  //
}

//____________________________________________________________________________
Bool_t AliReducedTrackCut::IsSelected(TObject* obj) {
  //
  // apply cuts
  //
  if(!obj->InheritsFrom(AliReducedBaseTrack::Class())) return kFALSE;
  
  //Fill values
  Float_t values[AliReducedVarManager::kNVars];
  AliReducedVarManager::FillTrackInfo((AliReducedBaseTrack*)obj, values);
  
  return IsSelected(obj, values);
}


//____________________________________________________________________________
Bool_t AliReducedTrackCut::IsSelected(TObject* obj, Float_t* values) {
   //
   // apply cuts
   //      
   if(!obj->InheritsFrom(AliReducedBaseTrack::Class())) return kFALSE;
   
   if(obj->InheritsFrom(AliReducedTrackInfo::Class())) {
      AliReducedTrackInfo* track = (AliReducedTrackInfo*)obj;
      if(fRequestITSrefit && !track->CheckTrackStatus(AliReducedVarManager::kITSrefit)) return kFALSE;
      if(fRequestTPCrefit && !track->CheckTrackStatus(AliReducedVarManager::kTPCrefit)) return kFALSE;
      if(fRequestTOFout && !track->CheckTrackStatus(AliReducedVarManager::kTOFout)) return kFALSE;
      if(fRequestCutOnITShitMap) {
         UChar_t itsHitMap = track->ITSclusterMap();
         UChar_t eval = itsHitMap & fCutOnITShitMap;
         if(fUseANDonITShitMap && (eval!=fCutOnITShitMap)) return kFALSE;
         if(!fUseANDonITShitMap && (eval==0)) return kFALSE;
      }
   }
   //if(fRejectKinks && (((AliReducedBaseTrack*)obj)->IsKink(0) || ((AliReducedBaseTrack*)obj)->IsKink(1) || ((AliReducedBaseTrack*)obj)->IsKink(2))) return kFALSE;
   if(fRejectKinks && (((AliReducedBaseTrack*)obj)->IsKink(0))) return kFALSE;
   if(fRejectTaggedGamma && ((AliReducedBaseTrack*)obj)->IsGammaLeg()) return kFALSE;
   if(fRejectTaggedPureGamma && ((AliReducedBaseTrack*)obj)->IsPureGammaLeg()) return kFALSE;
   
   return AliReducedVarCut::IsSelected(values);   
}
