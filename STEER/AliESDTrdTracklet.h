#ifndef ALIESDTRDTRACKLET_H
#define ALIESDTRDTRACKLET_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// ESD format for TRD tracklet from FEE used for triggering

#include "TObject.h"

class AliESDTrdTracklet : public TObject
{
 public:
  AliESDTrdTracklet();
  AliESDTrdTracklet(UInt_t trackletWord, Short_t hcid, Int_t label = -1);
  AliESDTrdTracklet(const AliESDTrdTracklet &trkl);
  AliESDTrdTracklet& operator=(const AliESDTrdTracklet &trkl);
  ~AliESDTrdTracklet();

  void SetTrackletWord(UInt_t trklWord) { fTrackletWord = trklWord; }
  void SetHCId(Short_t hcid) { fHCId = hcid; }

  // ----- tracklet information -----
  UInt_t GetTrackletWord() const { return fTrackletWord; }
  Int_t  GetBinY()  const;
  Int_t  GetBinDy() const;
  Int_t  GetBinZ()  const { return ((fTrackletWord >> 20) & 0xf);  }
  Int_t  GetPID()   const { return ((fTrackletWord >> 24) & 0xff); }

  // ----- position information (chamber-local) -----
  Float_t GetLocalX() const { return fgkX0; }
  Float_t GetLocalY() const { return GetBinY() * fgkBinWidthY; }
  Float_t GetDyDx() const;

  // ----- geometrical information -----
  Int_t GetHCId() const { return fHCId; }
  Int_t GetDetector() const { return fHCId / 2; }
  Int_t GetROB() const { return -1; }
  Int_t GetMCM() const { return -1; }

  // ----- MC information -----
  Int_t GetLabel() const { return fLabel; }

 protected:
  Short_t fHCId;		// half-chamber ID

  UInt_t fTrackletWord;		// tracklet word (as from FEE)
				// pppp : pppp : zzzz : dddd : dddy : yyyy : yyyy : yyyy
  Int_t  fLabel;		// MC label

  static const Float_t fgkBinWidthY;   // bin width y-position
  static const Float_t fgkBinWidthDy;  // bin width deflection length

  static const Float_t fgkX0;          // reference position in x-direction
  static const Float_t fgkDriftLength; // drift length to which the deflection length is scaled

  ClassDef(AliESDTrdTracklet, 1);
};

#endif
