#ifndef TOFv4_H
#define TOFv4_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////
//  Manager and hits classes for set:TOF  version 4  //
///////////////////////////////////////////////////////
 
#include "AliTOF.h"
#include "AliHit.h"
 
 
class AliTOFv4 : public AliTOF {

private:
  Int_t fIdFTOA;
  Int_t fIdFTOB;
  Int_t fIdFTOC;
  Int_t fIdFLTA;
  Int_t fIdFLTB;
  Int_t fIdFLTC;
 
public:
  AliTOFv4();
  AliTOFv4(const char *name, const char *title);
  virtual       ~AliTOFv4() {}
  virtual void   BuildGeometry();
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 4;}
  virtual void   TOFpc(Float_t,Float_t,Float_t,Float_t,Float_t,Float_t);
  virtual void   StepManager();
  virtual void   DrawModule();
 
   ClassDef(AliTOFv4,1)  //Time Of Flight version 4
};
 
#endif
