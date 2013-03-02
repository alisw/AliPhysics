#ifndef ALITPCV1_H
#define ALITPCV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Version 1 for TPC                         //
////////////////////////////////////////////////

 
#include "AliTPC.h"

class AliTPCv1 : public AliTPC {

public:
  AliTPCv1(); 
  AliTPCv1(const char *name, const char *title);
  virtual      ~AliTPCv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  StepManager();

protected:
  Int_t fIdSens;    //Sensitive volume identifier

  
  ClassDef(AliTPCv1,1)  // Time Projection Chamber version 1
};

#endif

