#ifndef ALIPMDV1_H
#define ALIPMDV1_H
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
  virtual void  DrawModule() const;

private:
  Int_t   fMedSens;        // Sensitive Medium Ar+CO2
  Float_t fDboxmm1[3];     // Master MODULE EMPA of aluminum for PMD
  Float_t fDboxmm12[3];    // Master MODULE EMCA of aluminum for CPV
  Float_t fDboxmm2[3];     // Master MODULE EMPB of aluminum for PMD
  Float_t fDboxmm22[3];    // Master MODULE EMCB of aluminum for CPV
 
  ClassDef(AliPMDv1,1)     //Hits manager for set:PMD
};
 
#endif
