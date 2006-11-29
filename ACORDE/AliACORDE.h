#ifndef ALIACORDE_H
#define ALIACORDE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: ACORDE        //
////////////////////////////////////////////////

#include "AliDetector.h"

class AliACORDEModule;

class AliACORDE : public AliDetector {
public:
  AliACORDE();
  AliACORDE(const char* name, const char* title);
  virtual ~AliACORDE();

  virtual void CreateMaterials();

  virtual Int_t IsVersion() const { return -1; }

  virtual TString Version() { return TString(""); }

  virtual void SetTreeAddress();
  virtual void SetModule(AliACORDEModule* module) {fModule = module;}
  virtual const AliACORDEModule* GetModule() const {return fModule; }
  virtual void MakeBranch(Option_t* opt = "");

protected:
  AliACORDEModule* fModule;
private:
  AliACORDE(const AliACORDE& crt);
  AliACORDE& operator=(const AliACORDE& crt);

  ClassDef(AliACORDE, 1) // Cosmic Ray Trigger (ACORDE) base class
};

typedef AliACORDE AliCRT; // for backward compatibility

#endif // ALIACORDE_H
