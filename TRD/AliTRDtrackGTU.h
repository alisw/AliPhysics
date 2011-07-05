#ifndef ALITRDTRACKGTU_H
#define ALITRDTRACKGTU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackGTU.h 27496 2008-07-22 08:35:45Z cblume $ */

//---------------------------------------------------------------
//
// TRD track as calculated in the GTU (for L1 contribution)
//
//---------------------------------------------------------------

#include "TClonesArray.h"

#include "AliTRDtrackletGTU.h"
class AliESDTrdTrack;

class AliTRDtrackGTU : public TObject {
 public:
  AliTRDtrackGTU();
  ~AliTRDtrackGTU();

// ----- Track properties
  Int_t    GetPtInt() const { return fPt; }
  Float_t  GetPt() const { return (Float_t) fPt / 128.; }
  Int_t    GetPID() const { return fPID; }
  Int_t    GetSector() const { return fSector; }
  Int_t    GetStack() const { return fStack; }
  Int_t GetLabel() const { return fLabel; }

  AliESDTrdTrack* CreateTrdTrack() const;

// ----- compositing tracklets
  Int_t    GetNTracklets() const;
  Int_t    GetTrackletMask() const { return fTrackletMask; }
  Bool_t   IsTrackletInLayer(Int_t layer) const;
  Int_t    GetTrackletIndex(Int_t layer) { return IsTrackletInLayer(layer) ? ((AliTRDtrackletGTU*) (*fTracklets)[layer])->GetIndex() : -1; }
  AliTRDtrackletGTU* GetTracklet(Int_t layer);

// ----- Quantities used internally for the calculation
  Float_t GetA() const { return fA; }
  Float_t GetB() const { return fB; }
  Float_t GetC() const { return fC; }
  Int_t GetZChannel() const { return fZChannel; }
  Int_t GetZSubChannel();
  Int_t GetRefLayerIdx() const { return fRefLayerIdx; }
  Int_t GetYapprox();


  void AddTracklet(const AliTRDtrackletGTU * const tracklet, Int_t layer);

  void SetStack(Int_t stack) { fStack = stack; }
  void SetSector(Int_t sector) { fSector = sector; }
  void SetPtInt(Int_t pt) { fPt = pt; }
  void SetPID(Int_t pid) { fPID = pid; }

  void SetZChannel(Int_t zch) { fZChannel = zch; }
  void SetRefLayerIdx(Int_t reflayer) { fRefLayerIdx = reflayer; }
//  void SetInnerIntPoint(Float_t *x);
//  void SetOuterIntPoint(Float_t *x);
  void SetFitParams(Float_t a, Float_t b, Float_t c);

  Bool_t CookLabel();

 protected:

  Int_t fStack; // TRD stack to which this track belongs
  Int_t fSector; // sector in which the track was found

  Int_t fPt; // pt in integer representation
  Int_t fPID; // PID calculated from tracklet PID

  TClonesArray *fTracklets; // array holding the tracklets composing this track
  Int_t fTrackletMask; // mask in which layers tracklets have been assigned
  Int_t fNTracklets; // number of tracklets in this track

  Int_t fRefLayerIdx; // index of the reference layer in which this track was found
  Int_t fZChannel; // z-channel unit in which this track was found
  Int_t fZSubChannel; // z-subchannel of the assigned tracklets

  Float_t fA; // fit parameter of y' = a + b*x + c*z
  Float_t fB; // fit parameter of y' = a + b*x + c*z
  Float_t fC; // fit parameter of y' = a + b*x + c*z

  Int_t fLabel; // MC label

 private:
  AliTRDtrackGTU(const AliTRDtrackGTU &rhs); // not implemented
  AliTRDtrackGTU& operator=(const AliTRDtrackGTU &rhs); // not implemented

  ClassDef(AliTRDtrackGTU, 1);
};

#endif
