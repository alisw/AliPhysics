#ifndef TRDrecPoint_h
#define TRDrecPoint_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliRecPoint.h"

class AliTRDrecPoint : public AliRecPoint {

 public:

  AliTRDrecPoint();
  virtual ~AliTRDrecPoint() {};
  virtual void  Print(Option_t * opt = "void") {};
  virtual void  AddDigit(Int_t digit);

  virtual void  SetAmplitude(Float_t amp)       { fAmp      = amp; };
  virtual void  SetDetector(Int_t det)          { fDetector = det; };
  virtual void  SetLocalPosition(TVector3 &pos);

  virtual Int_t GetDetector()                   { return fDetector; };

 protected:

  Int_t    fDetector;        // TRD detector number

  ClassDef(AliTRDrecPoint,1) // Reconstructed point for the TRD

};

#endif
