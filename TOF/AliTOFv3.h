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
  Int_t fIdFTO2; // First sensitive volume identifier
  Int_t fIdFTO3; // Second sensitive volume identifier
  Int_t fIdFLT1; // Third sensitive volume identifier
  Int_t fIdFLT2; // Fourth sensitive volume identifier
  Int_t fIdFLT3; // Fifth sensitive volume identifier
 
public:
  AliTOFv3();
  AliTOFv3(const char *name, const char *title);
  virtual       ~AliTOFv3() {}
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
