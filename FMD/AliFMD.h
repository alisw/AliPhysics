#ifndef ALIFMD_H
#define ALIFMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager and hits classes for set:Si-FMD     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "TString.h"
 
 class TFile;
 class TTree;
 class AliFMD : public AliDetector {
 
public:
  AliFMD();
  AliFMD(const char *name, const char *title);
  virtual       ~AliFMD(); 
  virtual void   AddHit(Int_t, Int_t*, Float_t*);
  virtual void   AddDigit(Int_t*);
  virtual void   BuildGeometry();
  virtual void   CreateGeometry() {}
  virtual void   CreateMaterials()=0; 
  virtual Int_t  DistanceToPrimitive(Int_t px, Int_t py);
  virtual Int_t  IsVersion() const =0;
  virtual void   Init();
  virtual void   MakeBranch(Option_t *opt=" ",char *file=0);
  virtual void   SetTreeAddress();
  virtual void   ResetHits();
  virtual void   ResetDigits();
  virtual void   DrawDetector()=0;
  virtual void   StepManager() {}
  void  Eta2Radius(Float_t, Float_t, Float_t*);
  void Hits2SDigits();//
   // Digitisation
  TClonesArray *SDigits() const {return fSDigits;}
//  virtual void   SDigits2Digits();
	
 protected:
  Int_t fIdSens1;     //Si sensetive volume
  TClonesArray *fSDigits      ; // List of summable digits
  ClassDef(AliFMD,2)  //Class for the FMD detector
};
#endif // AliFMD_H
