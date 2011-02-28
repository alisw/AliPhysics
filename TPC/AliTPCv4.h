#ifndef ALITPCV4_H
#define ALITPCV4_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////
//  Version 4 for TPC                         //
////////////////////////////////////////////////

 
#include "AliTPC.h"
#include <stdlib.h>
#include <TMath.h>
#include "AliMC.h"
#include "AliConst.h"
#include <TVirtualMC.h>
#include <TSystem.h>
#include "AliTPCParamSR.h"
#include "AliRun.h"
#include "AliTPCDigitsArray.h"
#include "TGeoManager.h"
class AliTPCv4 : public AliTPC {

public:
  AliTPCv4():AliTPC(),
  fIdSens(0),
  fIDrift(0),
  fSecOld(0){}
  AliTPCv4(const char *name, const char *title);
  virtual      ~AliTPCv4() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual void  AddAlignableVolumes() const;
  void SetInnerChambersAlignable() const;
  void SetOuterChambersAlignable() const;
  virtual void  Init();
  virtual Int_t IsVersion() const {return 4;}
  virtual void  StepManager();

protected:
  Int_t fIdSens;    // sensitive strip
  Int_t fIDrift;    // drift gas
  Int_t fSecOld;    // indicate the previous sector - for reference points    

  ClassDef(AliTPCv4,1)  // Time Projection Chamber version 1
};

#endif
