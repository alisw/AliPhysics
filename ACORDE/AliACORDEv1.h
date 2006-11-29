#ifndef ALIACORDEV1_H
#define ALIACORDEV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager class for detector: ACORDEv1         //
////////////////////////////////////////////////

#include "AliACORDE.h"

class AliACORDEv1 : public AliACORDE {
public:
  AliACORDEv1();
  AliACORDEv1(const char *name, const char *title);
  virtual ~AliACORDEv1();

  virtual TString Version() { return TString("v1"); }
  virtual Int_t IsVersion() const { return 1; }

  virtual void AddHit(Int_t track, Int_t *vol, Float_t *hits);
  //virtual void    FinishEvent();
  //virtual void    ResetHits();
  //virtual void    ResetDigits();

  virtual void CreateMaterials();
  virtual void CreateGeometry();

  virtual void Init();
  virtual void DrawDetector() const;
  virtual void StepManager();

protected:
  virtual void CreateMolasse();
  virtual void CreateShafts();

private: 
  AliACORDEv1(const AliACORDEv1& crt);
  AliACORDEv1& operator=(const AliACORDEv1& crt);

  ClassDef(AliACORDEv1, 1)  //Class for ACORDE, version 1, Shafts outside of AliHALL
};

typedef AliACORDEv1 AliCRTv1; // for backward compatibility

#endif // ALIACORDEV1_H
