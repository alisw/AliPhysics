#ifndef PMDV1_H
#define PMDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Rectangular geometry - Bedanga Mohanty - Spetember 2003

////////////////////////////////////////////////
//  Manager and hits classes for set:PMD      //
////////////////////////////////////////////////
 
#include "AliPMD.h"

//___________________________________________
 
class AliPMDv1 : public AliPMD {

private:
  Int_t fMedSens;
  Int_t fMedSens1;
  Float_t dbox_mm1[3];
  Float_t dbox_mm12[3];
  Float_t dbox_mm2[3];
  Float_t dbox_mm22[3];
  
public:
  AliPMDv1();
  AliPMDv1(const char *name, const char *title);
  virtual      ~AliPMDv1() {}
  virtual void  CreateGeometry();
  virtual void  CreatePMD();
  virtual void  CreateSupermodule();
  virtual void  GetParameters();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  StepManager();
  virtual void  DrawModule();
 
   ClassDef(AliPMDv1,1)  //Hits manager for set:PMD
};
 
#endif


