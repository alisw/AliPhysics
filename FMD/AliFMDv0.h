#ifndef ALIFMDV0_H
#define ALIFMDV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////
//  Manager and hits classes for set:FMD     //
////////////////////////////////////////////////
 
#include "AliFMD.h"
 
class AliFMDv0 : public AliFMD {
  
public:
  AliFMDv0() {};
  AliFMDv0(const char *name, const char *title);
  virtual       ~AliFMDv0() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawDetector();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void   StepManager();
  
protected:
  Int_t fIdSens1; // Sensetive volume  in FMD
  Int_t fIdSens2; // Sensetive volume  in FMD
  Int_t fIdSens3; // Sensetive volume  in FMD
  Int_t fIdSens4; // Sensetive volume  in FMD
  Int_t fIdSens5; // Sensetive volume  in FMD
   

   
  ClassDef(AliFMDv0,2)  //Class for FMD version 0
};

#endif


