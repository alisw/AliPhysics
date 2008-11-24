#ifndef ALIACORDEV1_H
#define ALIACORDEV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliACORDEv1.h,v 1.3 2007/08/18 08:40:00 hristov Exp $ */
////////////////////////////////////////////////
//  Manager class for detector: ACORDEv1         //
////////////////////////////////////////////////

#include "AliACORDE.h"

class AliACORDEv1 : public AliACORDE {
public:
  AliACORDEv1();
  AliACORDEv1(const char *name, const char *title);
  virtual void AddAlignableVolumes() const;

  virtual ~AliACORDEv1();

  virtual TString Version() { return TString("v1"); }
  virtual Int_t IsVersion() const { return 1; }
  virtual void AddHit(Int_t track, Int_t *vol, Float_t *hits);
  virtual void AddDigits(Int_t* track, Int_t module, Float_t time);
  virtual void   MakeBranch(Option_t *option);
  virtual void CreateGeometry();

  virtual void Init();
  virtual void StepManager();


protected:

  virtual void CreateAcorde();

private: 
  AliACORDEv1(const AliACORDEv1& crt);
  AliACORDEv1& operator=(const AliACORDEv1& crt);

  ClassDef(AliACORDEv1,2) // Cosmic Ray Trigger (ACORDE).
};

typedef AliACORDEv1 AliCRTv1; // for backward compatibility

#endif // ALIACORDEV1_H
