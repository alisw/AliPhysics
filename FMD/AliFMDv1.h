#ifndef ALIFMDV1_H
#define ALIFMDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////
//  Manager and hits classes for set:FMD     //
////////////////////////////////////////////////
 
#include "AliFMD.h"
 
class AliFMDv1 : public AliFMD {
  
public:
  AliFMDv1() {};
  AliFMDv1(const char *name, const char *title);
  virtual       ~AliFMDv1() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawDetector();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();
  //  virtual void Hit2Digits(Int_t bgrEvent, Option_t *opt1=" ",
  // 	Option_t *opt2=" ",Text_t *name=" "); // hit to digit for v1 :test  
 virtual void  Response( Float_t Edep);
//private:
 //Int_t fCharge; 


protected:
   Int_t fIdSens1; // Sensetive volume  in FMD
   
// Background event for event mixing
  
  ClassDef(AliFMDv1,2)  //Class for FMD version 0
};

#endif


