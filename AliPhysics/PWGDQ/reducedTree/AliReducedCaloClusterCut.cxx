/*
 ***********************************************************
 Implementation of AliReducedCaloClusterCut class.
 Contact: lucasaltenkamper@cern.ch
 13/04/2019
 *********************************************************
 */

#ifndef ALIREDUCEDCALOCLUSTERCUT_H
#include "AliReducedCaloClusterCut.h"
#endif

#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedVarManager.h"

ClassImp(AliReducedCaloClusterCut)

//____________________________________________________________________________
AliReducedCaloClusterCut::AliReducedCaloClusterCut() :
  AliReducedVarCut(),
  fDoTrackMatch(kFALSE),
  fMaxDistanceTrackMatch(-999.)
{
  //
  // default constructor
  //
}

//____________________________________________________________________________
AliReducedCaloClusterCut::AliReducedCaloClusterCut(const Char_t* name, const Char_t* title) :
  AliReducedVarCut(name, title),
  fDoTrackMatch(kFALSE),
  fMaxDistanceTrackMatch(-999.)
{
  //
  // named constructor
  //
}

//____________________________________________________________________________
AliReducedCaloClusterCut::~AliReducedCaloClusterCut() {
  //
  // destructor
  //
}

//____________________________________________________________________________
Bool_t AliReducedCaloClusterCut::IsSelected(TObject* obj) {
  //
  // apply cuts
  //
  if(!obj->InheritsFrom(AliReducedCaloClusterInfo::Class())) return kFALSE;
  
  //Fill values
  Float_t values[AliReducedVarManager::kNVars];
  AliReducedVarManager::FillCaloClusterInfo((AliReducedCaloClusterInfo*)obj, values);
  
  return IsSelected(obj, values);
}

//____________________________________________________________________________
Bool_t AliReducedCaloClusterCut::IsSelected(TObject* obj, Float_t* values) {
  //
  // apply cuts
  //
  if(!obj->InheritsFrom(AliReducedCaloClusterInfo::Class())) return kFALSE;
  
  // selection based on track match
  // NOTE: not sure about the distance cut, also could probably do something more sophisticated
  if (fDoTrackMatch) {
    Short_t nMatchedTracks = ((AliReducedCaloClusterInfo*)obj)->NMatchedTracks();
    if (nMatchedTracks<1) return kFALSE;
    Float_t deltaPhi = ((AliReducedCaloClusterInfo*)obj)->Dx();
    Float_t deltaEta = ((AliReducedCaloClusterInfo*)obj)->Dz();
    Float_t distance = TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
    if (distance>fMaxDistanceTrackMatch) return kFALSE;
  }
  
  return AliReducedVarCut::IsSelected(values);
}
