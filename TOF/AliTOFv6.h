#ifndef TOFv6_H
#define TOFv6_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////
//  Manager and hits classes for set:TOF  version 1  //
///////////////////////////////////////////////////////
 
#include "AliTOF.h"
#include "AliHit.h"
 
 
class AliTOFv6 : public AliTOF {

private:
  Int_t fIdFTO2; // First sensitive volume identifier
  Int_t fIdFTO3; // Second sensitive volume identifier
  Int_t fIdFLT1; // Third sensitive volume identifier
  Int_t fIdFLT2; // Fourth sensitive volume identifier
  Int_t fIdFLT3; // Fifth sensitive volume identifier
 
public:
  AliTOFv6();
  AliTOFv6(const char *name, const char *title);
  virtual       ~AliTOFv6() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 6;}
  virtual void   TOFpc(Float_t,Float_t,Float_t,Float_t,Float_t,Float_t);
  virtual void   StepManager();
  virtual void   DrawModule();
 
   ClassDef(AliTOFv6,1)  //Time Of Flight version 6
};
 
#endif
