#ifndef ALITRDRECPOINT_H
#define ALITRDRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliRecPoint.h"

class AliTRDrecPoint : public AliRecPoint {

 public:

  AliTRDrecPoint();
  virtual ~AliTRDrecPoint();

  virtual void    Print(Option_t* opt) {};
  virtual void    AddDigit(Int_t digit);
  virtual void    AddDigit(AliDigitNew &digit) {};

  virtual void    SetEnergy(Float_t amp)          { fAmp      = amp; };
  virtual void    SetDetector(Int_t det)          { fDetector = det; };
  virtual void    SetLocalPosition(TVector3 &pos);
  virtual void    SetLocalRow(Float_t r)          { fLocPos.SetX(r); };
  virtual void    SetLocalCol(Float_t c)          { fLocPos.SetY(c); };
  virtual void    SetLocalTime(Float_t t)         { fLocPos.SetZ(t); };

  virtual void    SetLocalTimeBin(Int_t tb)       { fTimeBin  = tb;  }
  virtual void    SetTrackingYZ(Float_t fSigmaY = 0.0, Float_t fSigmaZ = 0.0);  

  virtual Int_t   GetDetector() const             { return fDetector; };
  virtual Int_t   GetDigit(Int_t i = 0) const     { if (i < fMulDigit)
                                                      return fDigitsList[i]; 
                                                    else
                                                      return -1;};
  virtual Float_t GetLocalRow() const             { return fLocPos(0); };
  virtual Float_t GetLocalCol() const             { return fLocPos(1); };
  virtual Float_t GetLocalTime() const            { return fLocPos(2); };

  virtual Int_t   GetLocalTimeBin() const         { return Int_t(fLocPos(2)); }
  virtual Float_t GetSigmaY2() const              { return fSigmaY2; }
  virtual Float_t GetSigmaZ2() const              { return fSigmaZ2; }
  virtual Float_t GetY() const                    { return fY; }
  virtual Float_t GetZ() const                    { return fZ; }
          Int_t   IsUsed() const                  { return fUsed; }
          void    Use()                           { fUsed++; }
          Int_t   GetTrackIndex(Int_t i) const    { return fTracks[i]; }
          void    AddTrackIndex(Int_t *i);  

 protected:

  Int_t    fDetector;        // TRD detector number
  Int_t    fTimeBin;         // Time bin number within the detector
  Int_t    fUsed;            // 0 initially and incremented if the point is "used"
  Int_t    fTracks[3];       // labels of overlapped tracks
  Float_t  fY;               // local Rphi coordinate (cm) within tracking sector
  Float_t  fZ;               // local Z coordinate (cm) within tracking sector
  Float_t  fSigmaY2;         // Y variance (cm)
  Float_t  fSigmaZ2;         // Z variance (cm)  

  ClassDef(AliTRDrecPoint,1) // Reconstructed point for the TRD

};

#endif
