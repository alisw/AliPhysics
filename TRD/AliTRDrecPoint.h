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

  virtual Int_t   GetDetector()                   { return fDetector; };
  virtual Int_t   GetDigit(Int_t i = 0)           { if (i < fMulDigit)
                                                      return fDigitsList[i]; 
                                                    else
                                                      return -1;};
  virtual Float_t GetLocalRow()                   { return fLocPos(0); };
  virtual Float_t GetLocalCol()                   { return fLocPos(1); };
  virtual Float_t GetLocalTime()                  { return fLocPos(2); };

 protected:

  Int_t    fDetector;        // TRD detector number

  ClassDef(AliTRDrecPoint,1) // Reconstructed point for the TRD

};

#endif
