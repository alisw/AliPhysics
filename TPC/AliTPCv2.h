#ifndef ALITPCV2_H
#define ALITPCV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Version 2 for TPC                         //
////////////////////////////////////////////////

 
#include "AliTPC.h"

class AliTPCv2 : public AliTPC {

public:
  AliTPCv2() {}
  AliTPCv2(const char *name, const char *title);
  virtual      ~AliTPCv2() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 2;}
  virtual void  StepManager();
  virtual void  DrawDetector();

protected:
  Int_t fIdSens1;    //First  sensitive volume identifier - lower sector
  Int_t fIdSens2;    //Second sensitive volume identifier - upper sector
  Int_t fIdSens3;    //Sensitive strip - lower sector
  Int_t fIdSens4;    //Sensitive strip - upper sector     

private:

  Float_t BetheBloch(Float_t bg);
  
  ClassDef(AliTPCv2,1)  // Time Projection Chamber version 2
};

#endif
