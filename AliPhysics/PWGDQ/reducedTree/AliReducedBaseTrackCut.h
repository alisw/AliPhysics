// Class for cutting on base track (kinematics) information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
//   10/09/2015

#ifndef ALIREDUCEDBASETRACKCUT_H
#define ALIREDUCEDBASETRACKCUT_H

#include "AliReducedInfoCut.h"

//_________________________________________________________________________
class AliReducedBaseTrackCut : public AliReducedInfoCut {

 public:
  AliReducedBaseTrackCut();
  AliReducedBaseTrackCut(const Char_t* name, const Char_t* title);
  virtual ~AliReducedBaseTrackCut();
  
  virtual Bool_t IsSelected(TObject* obj);
  
  void SetPtRange(Float_t min, Float_t max) {fPtRange[0] = min; fPtRange[1]=max; fCutOnPt=kTRUE;}
  void SetEtaRange(Float_t min, Float_t max) {fEtaRange[0] = min; fEtaRange[1]=max; fCutOnEta=kTRUE;}
  void SetPhiRange(Float_t min, Float_t max) {fPhiRange[0] = min; fPhiRange[1]=max; fCutOnPhi=kTRUE;}
  
 protected: 
   
  Float_t  fPtRange[2];
  Bool_t   fCutOnPt;
  
  Float_t  fEtaRange[2];
  Bool_t   fCutOnEta;
  
  Float_t  fPhiRange[2];
  Bool_t   fCutOnPhi;
  
  AliReducedBaseTrackCut(const AliReducedBaseTrackCut &c);
  AliReducedBaseTrackCut& operator= (const AliReducedBaseTrackCut &c);
  
  ClassDef(AliReducedBaseTrackCut,1);
};

#endif