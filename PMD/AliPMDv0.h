#ifndef PMDV0_H
#define PMDV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:PMD      //
////////////////////////////////////////////////
 
#include "AliPMD.h"

//___________________________________________
 
class AliPMDv0 : public AliPMD {

private:
  Int_t fMedSens;
  
public:
  AliPMDv0();
  AliPMDv0(const char *name, const char *title);
  virtual      ~AliPMDv0() {}
  virtual void  CreateGeometry();
  virtual void  CreatePMD();
  virtual void  CreateSupermodule();
  virtual void  GetParameters();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  StepManager();
  virtual void  DrawModule();
 
   ClassDef(AliPMDv0,1)  //Hits manager for set:PMD
};
 
#endif


