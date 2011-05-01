#ifndef ALITPCCALIBMATERIAL_H
#define ALITPCCALIBMATERIAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
#include "TObjArray.h"

class THnSparse;
class TList;
class AliESDEvent;
class AliESDtrack;

#include "TTreeStream.h"
#include "TMap.h"
 

class AliTPCcalibMaterial:public AliTPCcalibBase {
public:
  AliTPCcalibMaterial(); 
  AliTPCcalibMaterial(const Text_t *name, const Text_t *title); 
  virtual ~AliTPCcalibMaterial();
  virtual void           Process(AliESDEvent *event);
  virtual void           ProcessPairs(AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *li);
  virtual void           Analyze(){;}
 public:
  Bool_t CheckLooper(Int_t index, AliESDEvent *event);
  Bool_t CheckV0(Int_t index, AliESDEvent *event);
 
  THnSparse * MakeHisto();
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}
  THnSparse * fHisMaterial;
  THnSparse * fHisMaterialRPhi;

private:
   AliTPCcalibMaterial(const AliTPCcalibMaterial&); // Not implemented
   AliTPCcalibMaterial& operator=(const AliTPCcalibMaterial&); // Not implemented
 
  ClassDef(AliTPCcalibMaterial, 1); 
};

#endif


