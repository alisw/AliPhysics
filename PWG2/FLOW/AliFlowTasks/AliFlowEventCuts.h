/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliFlowEventCuts:
// An event cut class
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#ifndef ALIFLOWEVENTCUTS_H
#define ALIFLOWEVENTCUTS_H

#include <float.h>
#include "TNamed.h"

class AliVEvent;;

class AliFlowEventCuts : public TNamed {

 public:
  AliFlowEventCuts();
  AliFlowEventCuts(const char* name, const char* title = "AliFlowEventCuts");
  //AliFlowEventCuts(const AliFlowEventCuts& someCuts);
  //AliFlowEventCuts& operator=(const AliFlowEventCuts& someCuts);
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

  Int_t GetNumberOfTracksMax() const {return fNumberOfTracksMax;}
  Int_t GetNumberOfTracksMin() const {return fNumberOfTracksMin;}
  Int_t GetRefMultMax() const {return fRefMultMax;}
  Int_t GetRefMultMin() const {return fRefMultMin;}

  static Int_t ReferenceMultiplicity(const AliVEvent* event);
  Int_t GetReferenceMultiplicity() {return fReferenceMultiplicity;}

 private:
  Bool_t fCutNumberOfTracks;//cut on # of tracks
  Int_t fNumberOfTracksMax;  //limits
  Int_t fNumberOfTracksMin;  //limits
  Bool_t fCutRefMult; //cut on refmult
  Int_t fRefMultMax; //max refmult
  Int_t fRefMultMin; //min refmult

  Int_t fReferenceMultiplicity; //store the reference multiplicity

  ClassDef(AliFlowEventCuts,1)
};

#endif


