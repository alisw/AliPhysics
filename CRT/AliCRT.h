#ifndef ALICRT_H
#define ALICRT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: ACORDE        //
////////////////////////////////////////////////

#include "AliDetector.h"

class AliCRTModule;

class AliCRT : public AliDetector {
public:
  AliCRT();
  AliCRT(const char* name, const char* title);
  virtual ~AliCRT();

  virtual void CreateMaterials();

  virtual Int_t IsVersion() const { return -1; }

  virtual TString Version() { return TString(""); }

  virtual void SetTreeAddress();
  virtual void SetModule(AliCRTModule* module) {fModule = module;}
  virtual const AliCRTModule* GetModule() const {return fModule; }
  virtual void MakeBranch(Option_t* opt = "");

protected:
  AliCRTModule* fModule;
private:
  AliCRT(const AliCRT& crt);
  AliCRT& operator=(const AliCRT& crt);

  ClassDef(AliCRT, 1) // Cosmic Ray Trigger (ACORDE) base class
};
#endif // ALICRT_H
