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

  virtual TString Version() { return TString("v0"); }
  virtual Int_t IsVersion() const { return 1; }
  virtual void AddHit(Int_t track, Int_t *vol, Float_t *hits);

  virtual void CreateGeometry();

  virtual void Init();
  virtual void StepManager();


protected:
  virtual void CreateCavern();
  virtual void CreateShafts();
  virtual void CreateMolasse();
  virtual void CreateAcorde();

private: 
  AliACORDEv0(const AliACORDEv0& crt);
  AliACORDEv0& operator=(const AliACORDEv0& crt);

  ClassDef(AliACORDEv0,1) // Cosmic Ray Trigger (ACORDE).
};

typedef AliACORDEv0 AliCRTv0; // for backward compatibility

#endif // ALIACORDEV0_H
