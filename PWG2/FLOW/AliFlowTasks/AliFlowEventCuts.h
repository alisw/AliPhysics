/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliFlowEventCuts:
// An event cut class
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#ifndef ALIFLOWEVENTCUTS_H
#define ALIFLOWEVENTCUTS_H

#include <float.h>
#include <limits.h>
#include "TNamed.h"

class AliVEvent;
class AliFlowTrackCuts;

class AliFlowEventCuts : public TNamed {

 public:
  enum refMultMethod { kTPConly, kSPDtracklets };

  AliFlowEventCuts();
  AliFlowEventCuts(const char* name, const char* title = "AliFlowEventCuts");
  AliFlowEventCuts(const AliFlowEventCuts& someCuts);
  AliFlowEventCuts& operator=(const AliFlowEventCuts& someCuts);
  virtual  ~AliFlowEventCuts() {}
  
  virtual Bool_t IsSelected(const TObject* obj);

  Bool_t PassesCuts(const AliVEvent* event);
  
  static AliFlowEventCuts* StandardCuts();
  
  void SetNumberOfTracksMax(Int_t value) {fNumberOfTracksMax=value;fCutNumberOfTracks=kTRUE;}
  void SetNumberOfTracksMin(Int_t value) {fNumberOfTracksMin=value;fCutNumberOfTracks=kTRUE;}
  void SetNumberOfTracksRange(Int_t min, Int_t max) {fNumberOfTracksMin=min;fNumberOfTracksMax=max;fCutNumberOfTracks=kTRUE;}
  void SetRefMultMax(Int_t value) {fRefMultMax=value;fCutRefMult=kTRUE;}
  void SetRefMultMin(Int_t value) {fRefMultMin=value;fCutRefMult=kTRUE;}
  void SetRefMultRange(Int_t min, Int_t max) {fRefMultMin=min;fRefMultMax=max;fCutRefMult=kTRUE;}
  void SetPrimaryVertexXrange(Double_t min, Double_t max)
       {fCutPrimaryVertexX=kTRUE; fPrimaryVertexXmin=min; fPrimaryVertexXmax=max;}
  void SetPrimaryVertexYrange(Double_t min, Double_t max)
       {fCutPrimaryVertexY=kTRUE; fPrimaryVertexYmin=min; fPrimaryVertexYmax=max;}
  void SetPrimaryVertexZrange(Double_t min, Double_t max)
       {fCutPrimaryVertexZ=kTRUE; fPrimaryVertexZmin=min; fPrimaryVertexZmax=max;}
  void SetNContributorsRange(Int_t min, Int_t max=INT_MAX) 
       {fCutNContributors=kTRUE; fNContributorsMin=min; fNContributorsMax=max;}

  Int_t GetNumberOfTracksMax() const {return fNumberOfTracksMax;}
  Int_t GetNumberOfTracksMin() const {return fNumberOfTracksMin;}
  Int_t GetRefMultMax() const {return fRefMultMax;}
  Int_t GetRefMultMin() const {return fRefMultMin;}
  void SetRefMultMethod(refMultMethod m) {fRefMultMethod=m;}
  refMultMethod GetRefMultMethod() const {return fRefMultMethod;}
  void SetRefMultCuts( AliFlowTrackCuts* cuts ) {fRefMultCuts=cuts;}
  AliFlowTrackCuts* GetRefMultCuts() const {return fRefMultCuts;}

  Int_t RefMult(const AliVEvent* event);
  //Int_t GetRefMult() {return fRefMult;}
  Int_t GetReferenceMultiplicity(const AliVEvent* event) {return RefMult(event);}

 private:
  Bool_t fCutNumberOfTracks;//cut on # of tracks
  Int_t fNumberOfTracksMax;  //limits
  Int_t fNumberOfTracksMin;  //limits
  Bool_t fCutRefMult; //cut on refmult
  refMultMethod fRefMultMethod; //how do we calculate refmult?
  Int_t fRefMultMax; //max refmult
  Int_t fRefMultMin; //min refmult
  AliFlowTrackCuts* fRefMultCuts; //cuts
  Bool_t fCutPrimaryVertexX; //cut on x of prim vtx
  Double_t fPrimaryVertexXmax; //max x prim vtx
  Double_t fPrimaryVertexXmin; //min x prim vtx
  Bool_t fCutPrimaryVertexY; //cut on y of prim vtx
  Double_t fPrimaryVertexYmax; //max y prim vtx
  Double_t fPrimaryVertexYmin; //min y prim vtx
  Bool_t fCutPrimaryVertexZ; //cut on z of prim vtx
  Double_t fPrimaryVertexZmax; //max z prim vtx
  Double_t fPrimaryVertexZmin; //min z prim vtx
  Bool_t fCutNContributors; //cut on number of contributors
  Int_t fNContributorsMax; //maximal number of contrib
  Int_t fNContributorsMin; //minimal number of contrib


  ClassDef(AliFlowEventCuts,2)
};

#endif


