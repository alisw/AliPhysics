#ifndef ALITRDTRACKLETGTU_H
#define ALITRDTRACKLETGTU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackletGTU.h 27496 2008-07-22 08:35:45Z cblume $ */

// --------------------------------------------------------
//
// GTU tracklet
//
//
// --------------------------------------------------------

#include "AliTRDtrackletBase.h"
#include "AliESDTrdTracklet.h"
#include "AliTRDgtuParam.h"
#include "AliLog.h"

class AliTRDtrackletGTU : public AliTRDtrackletBase {
 public:
  AliTRDtrackletGTU();
  AliTRDtrackletGTU(AliTRDtrackletBase *tracklet);
  AliTRDtrackletGTU(AliESDTrdTracklet *tracklet);
  AliTRDtrackletGTU(const AliTRDtrackletGTU& trk);

  ~AliTRDtrackletGTU();

  AliTRDtrackletGTU& operator=(const AliTRDtrackletGTU &rhs);

  // ----- Getters for information from the tracklet word -----
  Int_t GetYbin() const { return fTrackletESD ? fTrackletESD->GetBinY() : fTracklet->GetYbin(); }
  Int_t GetdY() const { return fTrackletESD ? fTrackletESD->GetBinDy() : fTracklet->GetdY(); }
  Int_t GetZbin() const { return fTrackletESD ? fTrackletESD->GetBinZ() : fTracklet->GetZbin(); }
  Int_t GetPID() const { return fTrackletESD ? fTrackletESD->GetPID() : ((Int_t) (256 * fTracklet->GetPID())); }
  Double_t GetPID(Int_t is) const { return fTracklet ? fTracklet->GetPID(is) : 0.; }

  // ----- Getters for calculated properties -----
  Int_t GetYProj() const { return fYProj; }
  Int_t GetAlpha() const { return fAlpha; }
  Int_t GetYPrime() const { return fYPrime; }
  Int_t GetSubChannel(Int_t zch) const;

  // ----- Getters for offline corresponding values -----
  Bool_t CookPID() { return kFALSE; }
  Int_t GetDetector() const { return fTrackletESD ? fTrackletESD->GetDetector() : fTracklet->GetDetector(); }
  Int_t GetIndex() const { return fIndex; }

  Float_t GetdYdX() const { return (GetdY() * 140e-4 / 3.); }
  Float_t GetX() const { return 0; }
  Float_t GetY() const { return (GetYbin() * 160e-4); }
  Float_t GetZ() const { return 0; }

  AliTRDtrackletBase* GetTracklet() const { return fTracklet; }
  AliESDTrdTracklet* GetTrackletESD() const { return fTrackletESD; }
  UInt_t GetTrackletWord() const { return fTrackletESD ? fTrackletESD->GetTrackletWord() : fTracklet->GetTrackletWord(); }
  Int_t GetMCMtrackletIndex() const { return fMCMtrackletIndex; }
  void  SetMCMtrackletIndex(Int_t val) { fMCMtrackletIndex=val; }

  Int_t GetSide() const { return GetYbin() < 0 ? 0 : 1; }

  Int_t GetLabel() const; // { return fLabel; }

  // ----- Setters -----
  void SetAlpha(Int_t alpha) { fAlpha = alpha; }
  void SetYProj(Int_t yproj) { fYProj = yproj; }
  void SetYPrime(Int_t yprime) { fYPrime = yprime; }

  void SetSubChannel(Int_t zch, Int_t subch);
  void SetDetector(Int_t /* id */ ) { AliError("Cannot change base tracklet"); }
  void SetIndex(Int_t idx) { fIndex = idx; }

  void RemoveTracklet() { fTracklet = fgkDummyTracklet; }

 protected:
  AliTRDgtuParam *fGtuParam;	 //!
  AliTRDtrackletBase *fTracklet;    //! pointer to the underlying tracklet
  AliESDTrdTracklet  *fTrackletESD; //! pointer to the underlying ESD tracklet

  Int_t fMCMtrackletIndex;      // Index number of the original tracklet in the TrackletTree

  Int_t  fSubChannel[AliTRDgtuParam::fgkNZChannels]; // z-channel assignments
  Bool_t fAssignedZ;		// tracklet assigned to a Z-channel

  Int_t fAlpha;			// calculated value for alpha
  Int_t fYProj;			// calculated value for y_proj
  Int_t fYPrime;		// calculated value for y'
  Int_t fIndex;                 // index of tracklet in the sequence after the input units

  static AliTRDtrackletBase* fgkDummyTracklet; // dummy tracklet, used in case no tracklet is given

 private:

  ClassDef(AliTRDtrackletGTU, 1);
};

#endif
