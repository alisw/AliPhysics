#ifndef ALIACORDEV0_H
#define ALIACORDEV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager class for detector: ACORDEv0         //
////////////////////////////////////////////////

#include "AliACORDE.h"

class AliACORDEv0 : public AliACORDE {
public:
  AliACORDEv0();
  AliACORDEv0(const char *name, const char *title);
  virtual ~AliACORDEv0();

  virtual void CreateGeometry();
  virtual void BuildGeometry();
  virtual void DrawDetector() const;

protected:
  virtual void CreateMolasse() {}
  virtual void CreateShafts() {}

private: 
  AliACORDEv0(const AliACORDEv0& crt);
  AliACORDEv0& operator=(const AliACORDEv0& crt);

  ClassDef(AliACORDEv0,1) // Cosmic Ray Trigger (ACORDE).
};
#endif // ALIACORDEV0_H
