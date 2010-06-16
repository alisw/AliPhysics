#ifndef ALITRDTRACKLETMCM_H
#define ALITRDTRACKLETMCM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackletMCM.h 27496 2008-07-22 08:35:45Z cblume $ */

//-----------------------------------
//
// TRD tracklet word (as from FEE)
// only 32-bit of information + detector ID
//
//----------------------------------

#include "AliTRDtrackletBase.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

class AliTRDtrackletMCM : public AliTRDtrackletBase {
 public:
  AliTRDtrackletMCM(UInt_t trackletWord = 0);
  AliTRDtrackletMCM(UInt_t trackletWword, Int_t hcid);
  AliTRDtrackletMCM(UInt_t trackletWword, Int_t hcid, Int_t rob, Int_t mcm);
  AliTRDtrackletMCM(const AliTRDtrackletMCM &rhs);
  ~AliTRDtrackletMCM();

  // ----- Getters for contents of tracklet word -----
  Int_t GetYbin() const; 
  Int_t GetdY() const; 
  Int_t GetZbin() const { return ((fTrackletWord >> 20) & 0xf); }
  Int_t GetPID() const { return ((fTrackletWord >> 24) & 0xff); }

  // ----- Getters for MCM-tracklet information -----
  Int_t GetMCM() const { return fMCM; }
  Int_t GetROB() const { return fROB; }
  Int_t GetLabel() const { return fLabel; }

  // ----- Getters for offline corresponding values -----
  Bool_t CookPID() { return kFALSE; }
  Double_t GetPID(Int_t /* is */) const { return GetPID()/255.; }
  Int_t GetDetector() const { return fHCId / 2; }
  Int_t GetHCId() const { return fHCId; }
  Float_t GetdYdX() const { return (GetdY() * 140e-4 / 3.); }
  Float_t GetX() const { return fGeo->GetTime0((fHCId % 12) / 2); }
  Float_t GetY() const { return (GetYbin() * 160e-4); }
  Float_t GetZ() const { return fGeo->GetPadPlane((fHCId % 12) / 2, (fHCId / 12) % 5)->GetRowPos( 4 * (fROB / 2) + fMCM / 4) - 
      fGeo->GetPadPlane((fHCId % 12) / 2, (fHCId /12) % 5)->GetRowSize(4 * (fROB / 2) + fMCM / 4) * .5; }
  Float_t GetLocalZ() const { return GetZ() - fGeo->GetPadPlane((fHCId % 12) / 2, (fHCId / 12) % 5)->GetRowPos(8); }

  Int_t GetQ0() const { return fQ0; }
  Int_t GetQ1() const { return fQ1; }
  Int_t GetNHits() const { return fNHits; }
  Int_t GetNHits0() const { return fNHits0; }
  Int_t GetNHits1() const { return fNHits1; }

  UInt_t GetTrackletWord() const { return fTrackletWord; }
  void SetTrackletWord(UInt_t trackletWord) { fTrackletWord = trackletWord; }

  void SetDetector(Int_t id) { fHCId = 2 * id + (GetYbin() < 0 ? 0 : 1); }
  void SetHCId(Int_t id) { fHCId = id; }
  void SetMCM(Int_t mcm) { fMCM = mcm; }
  void SetROB(Int_t rob) { fROB = rob; }
  void SetLabel(Int_t label) { fLabel = label; }
  void SetQ0(Int_t charge) { fQ0 = charge; }
  void SetQ1(Int_t charge) { fQ1 = charge; }
  void SetNHits(Int_t nhits) { fNHits = nhits; }
  void SetNHits0(Int_t nhits) { fNHits0 = nhits; }
  void SetNHits1(Int_t nhits) { fNHits1 = nhits; }

  void SetSlope(Float_t slope) { fSlope = slope; }
  void SetOffset(Float_t offset) { fOffset = offset; }
  void SetError(Float_t error) { fError = error; }
  void SetClusters(Float_t *res, Float_t *q, Int_t n);

  Float_t GetSlope() const { return fSlope; }
  Float_t GetOffset() const { return fOffset; }
  Float_t GetError() const { return fError; }
  Int_t   GetNClusters() const { return fNClusters; }
  Float_t *GetResiduals() const { return fResiduals; }
  Float_t *GetClsCharges() const { return fClsCharges; }

 protected:
  AliTRDgeometry *fGeo; //! TRD geometry

  Int_t fHCId;                  // half-chamber ID (only transient)
  UInt_t fTrackletWord;		// tracklet word: PID | Z | deflection length | Y 
				//          bits:  12   4            7          13
  Int_t fMCM; // MCM no. in which the tracklet was found
  Int_t fROB; // ROB no. on which the tracklet was found

  Int_t fQ0; // accumulated charge in the first time window
  Int_t fQ1; // accumulated charge in the second time window

  Int_t fNHits;  // no. of contributing clusters
  Int_t fNHits0; // no. of contributing clusters in window 0
  Int_t fNHits1; // no. of contributing clusters in window 1

  Int_t fLabel; // label for MC track
  
  Float_t  fSlope;	      // tracklet slope
  Float_t  fOffset;	      // tracklet offset
  Float_t  fError;            // tracklet error
  Int_t    fNClusters;	      // no. of clusters
  Float_t *fResiduals;	      //[fNClusters] cluster to tracklet residuals
  Float_t *fClsCharges;	      //[fNClusters] cluster charge

 private:
  AliTRDtrackletMCM& operator=(const AliTRDtrackletMCM &rhs);   // not implemented

  ClassDef(AliTRDtrackletMCM, 2);
};

#endif
