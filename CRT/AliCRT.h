#ifndef ALICRT_H
#define ALICRT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: ACORDE        //
////////////////////////////////////////////////

#include <TObject.h>
#include <TTree.h>

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
  virtual       ~AliCRT();

  virtual void  AddHit(Int_t track, Int_t *vol, Float_t *hits);
  virtual void  AddDigit( Int_t* tracks, Int_t* digits);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry();
  virtual void  Init() const;
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  DrawDetector() const {};
  virtual void  DrawModule() const {};
  virtual void  StepManager() = 0;
  virtual void  MakeBranch(Option_t *opt=" ", const char *file=0);

  virtual void  FinishEvent();
  virtual void  ResetHits();
  virtual void  ResetDigits();

private: 
   ClassDef(AliCRT,1)  //Class manager for CRT(ACORDE)
};

#endif // ALICRT_H
