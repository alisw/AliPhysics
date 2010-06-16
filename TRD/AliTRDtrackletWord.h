#ifndef ALITRDTRACKLETWORD_H
#define ALITRDTRACKLETWORD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackletWord.h 27496 2008-07-22 08:35:45Z cblume $ */

//-----------------------------------
//
// TRD tracklet word (as from FEE)
// only 32-bit of information + detector ID
//
//----------------------------------

#include "AliTRDtrackletBase.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

class AliTRDtrackletWord : public AliTRDtrackletBase {
 public:
  AliTRDtrackletWord(UInt_t trackletWord = 0);
  AliTRDtrackletWord(UInt_t trackletWord, Int_t hcid);
  AliTRDtrackletWord(const AliTRDtrackletWord &rhs);
  ~AliTRDtrackletWord();

  // ----- Getters for contents of tracklet word -----
  Int_t GetYbin() const; 
  Int_t GetdY() const; 
  Int_t GetZbin() const { return ((fTrackletWord >> 20) & 0xf); }
  Int_t GetPID() const { return ((fTrackletWord >> 24) & 0xff); }

  Int_t GetROB() const;
  Int_t GetMCM() const; 

  // ----- Getters for offline corresponding values -----
  Bool_t CookPID() { return kFALSE; }
  Double_t GetPID(Int_t /* is */) const { return 0; }
  Int_t GetDetector() const { return fHCId / 2; }
  Int_t GetHCId() const { return fHCId; }
  Float_t GetdYdX() const { return (GetdY() * 140e-4 / 3.); }
  Float_t GetX() const { return fgGeo->GetTime0((fHCId%12)/2); }
  Float_t GetY() const { return (GetYbin() * 160e-4); }
  Float_t GetZ() const { return fgGeo->GetPadPlane((fHCId % 12) / 2, (fHCId/12) % 5)->GetRowPos(GetZbin()) -
      fgGeo->GetPadPlane((fHCId % 12) / 2, (fHCId/12) % 5)->GetRowSize(GetZbin()); }
  Float_t GetLocalZ() const { return GetZ() - fgGeo->GetPadPlane((fHCId % 12) / 2, (fHCId/12) % 5)->GetRowPos(8); }

  UInt_t GetTrackletWord() const { return fTrackletWord; }
  void SetTrackletWord(UInt_t trackletWord) { fTrackletWord = trackletWord; }

  void SetDetector(Int_t id) { fHCId = 2 * id + (GetYbin() < 0 ? 0 : 1); }
  void SetHCId(Int_t id) { fHCId = id; }

 protected:
  Int_t fHCId;                  // half-chamber ID
  UInt_t fTrackletWord;		// tracklet word: PID | Z | deflection length | Y 
				//          bits:  12   4            7          13
  static AliTRDgeometry *fgGeo;  // pointer to TRD geometry for coordinate calculations

 private:
  AliTRDtrackletWord& operator=(const AliTRDtrackletWord &rhs);   // not implemented

  ClassDef(AliTRDtrackletWord, 2);
};

#endif
