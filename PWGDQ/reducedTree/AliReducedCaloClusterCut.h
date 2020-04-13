// Class for cutting on ALICE Var manager and other cluster specific information
// Author: Lucas Altenkamper (lucas.altenkamper@cern.ch)
//   13/04/2019

#ifndef ALIREDUCEDCALOCLUSTERCUT_H
#define ALIREDUCEDCALOCLUSTERCUT_H

#include "AliReducedVarCut.h"

//_________________________________________________________________________
class AliReducedCaloClusterCut : public AliReducedVarCut {

public:
  AliReducedCaloClusterCut();
  AliReducedCaloClusterCut(const Char_t* name, const Char_t* title);
  virtual ~AliReducedCaloClusterCut();

  // setters
  void SetMaxTrackMatchDistance(Float_t maxDist) {fDoTrackMatch=kTRUE; fMaxDistanceTrackMatch=maxDist;}
  
  // getters
  Float_t GetMaxTrackMatchDistance() const {return fMaxDistanceTrackMatch;}

  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(TObject* obj, Float_t* values);
  
protected:
  
  Bool_t  fDoTrackMatch;            // if true, cluster w/o track match within defined distance are rejected
  Float_t fMaxDistanceTrackMatch;   // max distance for cluster-track matching
  
  AliReducedCaloClusterCut(const AliReducedCaloClusterCut &c);
  AliReducedCaloClusterCut& operator= (const AliReducedCaloClusterCut &c);
  
  ClassDef(AliReducedCaloClusterCut, 1);
};

#endif
