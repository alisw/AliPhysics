#ifndef TOFv5_H
#define TOFv5_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////
//  Manager and hits classes for set:TOF  version 5  //
///////////////////////////////////////////////////////
 
#include "AliTOF.h"
#include "AliHit.h"
 
 
class AliTOFv5 : public AliTOF {

private:
  Int_t fIdFTO2; // First sensitive volume identifier
  Int_t fIdFTO3; // Second sensitive volume identifier
  Int_t fIdFLT1; // Third sensitive volume identifier
  Int_t fIdFLT2; // Fourth sensitive volume identifier
  Int_t fIdFLT3; // Fifth sensitive volume identifier
 
public:
  AliTOFv5();
  AliTOFv5(const char *name, const char *title);
  virtual       ~AliTOFv5() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 5;}
  virtual void   TOFpc(Float_t,Float_t,Float_t,Float_t,Float_t,Float_t);
  virtual void   StepManager();
  virtual void   DrawModule();
 
   ClassDef(AliTOFv5,1)  //Time Of Flight version 1
};
 
#endif
