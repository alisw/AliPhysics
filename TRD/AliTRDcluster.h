#ifndef ALITRDCLUSTER_H
#define ALITRDCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

class AliTRDrecPoint;

class AliTRDcluster : public TObject {

 public:

  AliTRDcluster();
  AliTRDcluster(AliTRDrecPoint *rp);
  
  virtual Int_t   GetDetector() const           { return fDetector; };
  virtual Int_t   GetLocalTimeBin() const       { return fTimeBin; }

  virtual Float_t GetSigmaY2() const            { return fSigmaY2; }
  virtual Float_t GetSigmaZ2() const            { return fSigmaZ2; }
  virtual Float_t GetY() const                  { return fY; }
  virtual Float_t GetZ() const                  { return fZ; }
  virtual Float_t GetQ() const                  { return fQ; }

  Int_t   IsUsed() const                        { return (fQ<0) ? 1 : 0; }
  void    Use()                                 { fQ=-fQ; }
  Int_t   GetTrackIndex(Int_t i) const          { return fTracks[i]; }


 protected:

  Int_t    fDetector;    // TRD detector number
  Int_t    fTimeBin;     // Time bin number within the detector

  Int_t    fTracks[3];   // labels of overlapped tracks
  Float_t  fQ;           // amplitude 
  Float_t  fY;           // local Rphi coordinate (cm) within tracking sector
  Float_t  fZ;           // local Z coordinate (cm) within tracking sector
  Float_t  fSigmaY2;     // Y variance (cm)
  Float_t  fSigmaZ2;     // Z variance (cm)  


  ClassDef(AliTRDcluster,1) // Reconstructed point for the TRD
 
};

#endif
