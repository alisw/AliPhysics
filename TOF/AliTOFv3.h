#ifndef TOFv3_H
#define TOFv3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////
//  Manager and hits classes for set:TOF  version 3  //
///////////////////////////////////////////////////////
 
#include "AliTOF.h"
#include "AliHit.h"
 
 
class AliTOFv3 : public AliTOF {

private:
  Int_t fIdFTOA;
  Int_t fIdFTOB;
  Int_t fIdFTOC;
  Int_t fIdFLTA;
  Int_t fIdFLTB;
  Int_t fIdFLTC;
 
public:
  AliTOFv3();
  AliTOFv3(const char *name, const char *title);
  virtual       ~AliTOFv3() {}
  virtual void   BuildGeometry();
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 3;}
  virtual void   TOFpc(Float_t,Float_t,Float_t,Float_t,Float_t,Float_t);
  virtual void   StepManager();
  virtual void   DrawModule();
 
   ClassDef(AliTOFv3,1)  //Time Of Flight version 3
};
 
#endif
