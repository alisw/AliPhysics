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
  Double_t GetPID(Int_t /* is */) const { return 0; }
  Int_t GetDetector() const { return fHCId / 2; }
  Int_t GetHCId() const { return fHCId; }
  Float_t GetdYdX() const { return (GetdY() * 140e-4 / 3.); }
  Float_t GetX() const { return fGeo->GetTime0((fHCId % 12) / 2); }
  Float_t GetY() const { return (GetYbin() * 160e-4); }
  Float_t GetZ() const { return fGeo->GetPadPlane((fHCId % 12) / 2, (fHCId / 12) % 5)->GetRowPos( 4 * (fROB / 2) + fMCM / 4); }

  UInt_t GetTrackletWord() const { return fTrackletWord; }
  void SetTrackletWord(UInt_t trackletWord) { fTrackletWord = trackletWord; }

  void SetDetector(Int_t id) { fHCId = 2 * id + (GetYbin() < 0 ? 0 : 1); }
  void SetHCId(Int_t id) { fHCId = id; }
  void SetMCM(Int_t mcm) { fMCM = mcm; }
  void SetROB(Int_t rob) { fROB = rob; }
  void SetLabel(Int_t label) { fLabel = label; }

 protected:
  AliTRDgeometry *fGeo; //! TRD geometry

  Int_t fHCId;                  // half-chamber ID (only transient)
  UInt_t fTrackletWord;		// tracklet word: PID | Z | deflection length | Y 
				//          bits:  12   4            7          13
  Int_t fMCM; // MCM no. in which the tracklet was found
  Int_t fROB; // ROB no. on which the tracklet was found

  Int_t fLabel; // label for MC track

 private:
  AliTRDtrackletMCM& operator=(const AliTRDtrackletMCM &rhs);   // not implemented

  ClassDef(AliTRDtrackletMCM, 1);
};

#endif
