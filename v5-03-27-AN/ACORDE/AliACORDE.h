#ifndef ALIACORDE_H
#define ALIACORDE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: ACORDE        //
////////////////////////////////////////////////

#include "AliDetector.h"
#include "AliACORDELoader.h"
#include "AliACORDEDigitizer.h"
#include "AliACORDETrigger.h"

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
  virtual void MakeBranch(Option_t* opt = "");

  virtual AliLoader* MakeLoader(const char* topfoldername);

  AliDigitizer*  CreateDigitizer(AliDigitizationInput* digInput) const;

  virtual AliTriggerDetector* CreateTriggerDetector() const
  { return new AliACORDETrigger(); }

  void  Digits2Raw ();
  virtual Bool_t Raw2SDigits(AliRawReader*);

  virtual void SetCreateCavern(Bool_t b) {fCreateCavern = b;}
  virtual void Set4CentralModulesGeometry(Bool_t b) {f4CentralModulesGeometry = b;}
  virtual Bool_t GetCreateCavern() const {return fCreateCavern;}
  virtual Bool_t Get4CentralModulesGeometry() const {return f4CentralModulesGeometry;}

private:
  AliACORDE(const AliACORDE& crt);
  AliACORDE& operator=(const AliACORDE& crt);

  Bool_t fCreateCavern;
  Bool_t f4CentralModulesGeometry;

  ClassDef(AliACORDE, 1) // Cosmic Ray Trigger (ACORDE) base class
};

typedef AliACORDE AliCRT; // for backward compatibility

#endif // ALIACORDE_H
