#ifndef PMDv2_H
#define PMDv2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////
//  Manager and hits classes for set:PMD  version 0     //
//////////////////////////////////////////////////////////
 
#include "AliPMD.h"
 
class AliPMDv2 : public AliPMD {

private:
  Int_t fMedSens;
  
public:
  AliPMDv2();
  AliPMDv2(const char *name, const char *title);
  virtual      ~AliPMDv2() {}
  virtual void  Coordinates();
  virtual void  Coordnew();
  virtual void  CreateGeometry();
  virtual void  CreateInside();
  virtual void  CreatePads();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 2;}
  virtual void  StepManager();
  virtual void  DrawModule();
 
   ClassDef(AliPMDv2,1)  //Hits manager for set:PMD
};
 
#endif
