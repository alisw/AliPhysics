#ifndef ALIFMDV2_H
#define ALIFMDV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////
//  Manager and hits classes for set:FMD     //
////////////////////////////////////////////////
 
#include "AliFMD.h"
 
class AliFMDv2 : public AliFMD {
  
public:
  AliFMDv2() {};
  AliFMDv2(const char *name, const char *title);
  virtual       ~AliFMDv2() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawDetector();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 1;}
  virtual void   StepManager();
  //  virtual void Hit2Digits(Int_t bgrEvent, Option_t *opt1=" ",
  // 	Option_t *opt2=" ",Text_t *name=" "); // hit to digit for v1 :test  
 virtual void  Response( Float_t Edep);
//private:
 //Int_t fCharge; 


protected:
   Int_t fIdSens1; // Sensetive volume  in FMD
   Int_t fIdSens2; // Sensetive volume  in FMD
   Int_t fIdSens3; // Sensetive volume  in FMD
   Int_t fIdSens4; // Sensetive volume  in FMD
   Int_t fIdSens5; // Sensetive volume  in FMD
   
// Background event for event mixing
  
  ClassDef(AliFMDv2,2)  //Class for FMD version 0
};

#endif


