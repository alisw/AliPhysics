#ifndef ALIFMD_H
#define ALIFMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager and hits classes for set:Si-FMD     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
 
class AliFMD : public AliDetector {
 
public:
  AliFMD();
  AliFMD(const char *name, const char *title);
  virtual       ~AliFMD(); 
  virtual void   AddHit(Int_t, Int_t*, Float_t*);
  virtual void   BuildGeometry();
  virtual void   CreateGeometry()=0;
  virtual void   CreateMaterials()=0; 
  virtual Int_t  DistanceToPrimitive(Int_t px, Int_t py);
  virtual Int_t  IsVersion() const =0;
  virtual void   Init();
  virtual void   MakeBranch(Option_t *opt=" ");
  virtual void   DrawDetector()=0;
  virtual void   StepManager()=0;
  void  Eta2Radius(Float_t, Float_t, Float_t*);
  
 protected:
  Int_t fIdSens1;     //Si sensetive volume
  ClassDef(AliFMD,1)  //Class for the FMD detector
};
#endif
