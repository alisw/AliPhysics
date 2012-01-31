/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliStarEventCuts:
// An event cut class for AliStarEvent
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#ifndef ALISTAREVENTCUTS_H
#define ALISTAREVENTCUTS_H

#include <float.h>
#include "TNamed.h"

class AliStarTrack;

class AliStarEventCuts : public TNamed {

 public:
  AliStarEventCuts();
  //AliStarEventCuts(const AliStarEventCuts& someCuts);
  //AliStarEventCuts& operator=(const AliStarEventCuts& someCuts);
  virtual  ~AliStarEventCuts() {}
  
  Bool_t PassesCuts(const AliStarEvent* event) const;
  static AliStarEventCuts* StandardCuts();
  
  void SetRunIDMax(const Int_t value) {fRunIDMax=value;fCutRunID=kTRUE;}
  void SetRunIDMin(const Int_t value) {fRunIDMin=value;fCutRunID=kTRUE;}
  void SetEventNumberMax(const Int_t value) {fEventNumberMax=value;fCutEventNumber=kTRUE;}
  void SetEventNumberMin(const Int_t value) {fEventNumberMin=value;fCutEventNumber=kTRUE;}
  void SetVtxXMax(const Float_t value) {fVtxXMax=value;fCutVtxX=kTRUE;}
  void SetVtxXMin(const Float_t value) {fVtxXMin=value;fCutVtxX=kTRUE;}
  void SetVtxYMax(const Float_t value) {fVtxYMax=value;fCutVtxY=kTRUE;}
  void SetVtxYMin(const Float_t value) {fVtxYMin=value;fCutVtxY=kTRUE;}
  void SetVtxZMax(const Float_t value) {fVtxZMax=value;fCutVtxZ=kTRUE;}
  void SetVtxZMin(const Float_t value) {fVtxZMin=value;fCutVtxZ=kTRUE;}
  void SetBFieldMax(const Float_t value) {fBFieldMax=value;fCutBField=kTRUE;}
  void SetBFieldMin(const Float_t value) {fBFieldMin=value;fCutBField=kTRUE;}
  void SetRefMultMax(const Int_t value) {fRefMultMax=value;fCutRefMult=kTRUE;}
  void SetRefMultMin(const Int_t value) {fRefMultMin=value;fCutRefMult=kTRUE;}
  void SetCentralityIDMax(const Int_t value) {fCentralityIDMax=value;fCutCentralityID=kTRUE;}
  void SetCentralityIDMin(const Int_t value) {fCentralityIDMin=value;fCutCentralityID=kTRUE;}
  void SetNumberOfPrimaryTracksMax(const Int_t value) {fNumberOfPrimaryTracksMax=value;fCutNumberOfPrimaryTracks=kTRUE;}
  void SetNumberOfPrimaryTracksMin(const Int_t value) {fNumberOfPrimaryTracksMin=value;fCutNumberOfPrimaryTracks=kTRUE;}
  void SetNumberOfTracksMax(const Int_t value) {fNumberOfTracksMax=value;fCutNumberOfTracks=kTRUE;}
  void SetNumberOfTracksMin(const Int_t value) {fNumberOfTracksMin=value;fCutNumberOfTracks=kTRUE;}

  Int_t GetRunIDMax() const {return fRunIDMax;}
  Int_t GetRunIDMin() const {return fRunIDMin;}
  Int_t GetEventNumberMax() const {return fEventNumberMax;}
  Int_t GetEventNumberMin() const {return fEventNumberMin;}
  Float_t GetVtxXMax() const {return fVtxXMax;}
  Float_t GetVtxXMin() const {return fVtxXMin;}
  Float_t GetVtxYMax() const {return fVtxYMax;}
  Float_t GetVtxYMin() const {return fVtxYMin;}
  Float_t GetVtxZMax() const {return fVtxZMax;}
  Float_t GetVtxZMin() const {return fVtxZMin;}
  Float_t GetBFieldMax() const {return fBFieldMax;}
  Float_t GetBFieldMin() const {return fBFieldMin;}
  Int_t GetRefMultMax() const {return fRefMultMax;}
  Int_t GetRefMultMin() const {return fRefMultMin;}
  Int_t GetCentralityIDMax() const {return fCentralityIDMax;}
  Int_t GetCentralityIDMin() const {return fCentralityIDMin;}
  Int_t GetNumberOfPrimaryTracksMax() const {return fNumberOfPrimaryTracksMax;}
  Int_t GetNumberOfPrimaryTracksMin() const {return fNumberOfPrimaryTracksMin;}
  Int_t GetNumberOfTracksMax() const {return fNumberOfTracksMax;}
  Int_t GetNumberOfTracksMin() const {return fNumberOfTracksMin;}

 private:
  Bool_t   fCutRunID; //cut on run id
  Int_t fRunIDMax;  //limits
  Int_t fRunIDMin;  //limits
  Bool_t   fCutEventNumber; //cut on event number
  Int_t fEventNumberMax;  //limits
  Int_t fEventNumberMin;  //limits
  Bool_t   fCutVtxX;//cut on vertex
  Float_t fVtxXMax;  //limits
  Float_t fVtxXMin;  //limits
  Bool_t   fCutVtxY;//cut on vertex
  Float_t fVtxYMax;  //limits
  Float_t fVtxYMin;  //limits
  Bool_t   fCutVtxZ;//cut on vertex
  Float_t fVtxZMax;  //limits
  Float_t fVtxZMin;  //limits
  Bool_t   fCutBField;//cut on bfield
  Float_t fBFieldMax;  //limits
  Float_t fBFieldMin;  //limits
  Bool_t   fCutRefMult;//cut on reference multiplicity
  Int_t fRefMultMax;  //limits
  Int_t fRefMultMin;  //limits
  Bool_t   fCutCentralityID;//cut on centrality id
  Int_t fCentralityIDMax;  //limits
  Int_t fCentralityIDMin;  //limits
  Bool_t   fCutNumberOfPrimaryTracks;//cut on # prim tracks
  Int_t fNumberOfPrimaryTracksMax;  //limits
  Int_t fNumberOfPrimaryTracksMin;  //limits
  Bool_t   fCutNumberOfTracks;//cut on # of tracks
  Int_t fNumberOfTracksMax;  //limits
  Int_t fNumberOfTracksMin;  //limits

  ClassDef(AliStarEventCuts,1)
};

#endif


