#ifndef ALITPCCALIBTRIGGER_H
#define ALITPCCALIBTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
#include "TH2F.h"
#include "TF1.h"
#include "TArrayD.h"
#include "TObjArray.h"

class TH1F;
class TH3F;
class TH2F;
class THnSparse;
class TList;
class AliESDEvent;
class AliESDtrack;
class AliTPCcalibLaser;
class TGraphErrors;

#include "TTreeStream.h"
#include "TMap.h"
 

class AliTPCcalibTrigger:public AliTPCcalibBase {
public:
  AliTPCcalibTrigger(); 
  AliTPCcalibTrigger(const Text_t *name, const Text_t *title); 
  virtual ~AliTPCcalibTrigger(){;}
  virtual void           Process(AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *li);
  virtual void           Analyze(){;}
  THnSparse * GetHisto(const char *trigger);
  void   AddHisto(const char *trigger, THnSparse *his); 
  THnSparse *MakeHisto(const char* trigger);
  //
public:
  TMap *fHisMap;      // map of the histogram per trigger class 
  ClassDef(AliTPCcalibTrigger, 1); 
};

#endif


