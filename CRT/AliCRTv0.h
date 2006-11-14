#ifndef ALICRTV0_H
#define ALICRTV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager class for detector: CRTv0         //
////////////////////////////////////////////////

#include "AliCRT.h"

class AliCRTv0 : public AliCRT {
public:
  AliCRTv0();
  AliCRTv0(const char *name, const char *title);
  virtual ~AliCRTv0();

  virtual void CreateGeometry();
  virtual void BuildGeometry();
  virtual void DrawDetector() const;

protected:
  virtual void CreateMolasse() {}
  virtual void CreateShafts() {}

private: 
  AliCRTv0(const AliCRTv0& crt);
  AliCRTv0& operator=(const AliCRTv0& crt);

  ClassDef(AliCRTv0,1) // Cosmic Ray Trigger (ACORDE).
};
#endif // ALICRTV0_H
