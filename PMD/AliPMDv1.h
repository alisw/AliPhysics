#ifndef PMDV1_H
#define PMDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:PMD      //
////////////////////////////////////////////////
 
#include "AliPMD.h"

//___________________________________________
 
class AliPMDv1 : public AliPMD {

private:
  Int_t fMedSens;
  
public:
  AliPMDv1();
  AliPMDv1(const char *name, const char *title);
  virtual      ~AliPMDv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateInside();
  virtual void  CreatePads();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  StepManager();
  virtual void  DrawModule();
 
   ClassDef(AliPMDv1,1)  //Hits manager for set:PMD
};
 
#endif
