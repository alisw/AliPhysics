#ifndef ALICRT_H
#define ALICRT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: ACORDE        //
////////////////////////////////////////////////

#include "AliDetector.h"

class TFile;
class TDirectory;
class TString ;  
class TTask ;
class TFolder ;

class AliCRT : public AliDetector {
 
public:
  
                AliCRT();
                AliCRT(const char *name, const char *title);
                AliCRT(const AliCRT& crt);
                AliCRT& operator= (const AliCRT& crt);
  virtual       ~AliCRT();

  virtual void  AddHit(Int_t track, Int_t *vol, Float_t *hits);
  virtual void  AddDigit(Int_t *tracks,Int_t *digits);

  virtual void  BuildGeometry();
  virtual void  CreateGeometry();
  virtual void  Init() const;
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 0;}
  virtual TString Version(void) {return TString("");}
  virtual void  DrawDetector() const {}
  virtual void  DrawModule() const {}
  virtual void  StepManager() {}

  virtual void  FinishEvent();
  virtual void  ResetHits();
  virtual void  ResetDigits();
  virtual void  SetTreeAddress();

protected:

private:
  ClassDef(AliCRT,1)  //Class manager for CRT(ACORDE)
};

#endif // ALICRT_H
