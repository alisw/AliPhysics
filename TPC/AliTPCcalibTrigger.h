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
  virtual ~AliTPCcalibTrigger();
  virtual void           Process(AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *li);
  virtual void           Analyze(){;}
  THnSparse * GetHisto(const char *trigger);
  void   AddHisto(const char *trigger, THnSparse *his); 
  THnSparse *MakeHisto(const char* trigger);
  //
  TTree * MakeTree(const char *fname);
  void   MakeTree(TTreeStream &pcstream, const char *tname);
  Bool_t HasTOF(TObjString *tname);
  Bool_t HasACORDE(TObjString *tname);
  Bool_t HasPIXEL(TObjString *tname);
  Int_t HasTRD(TObjString *tname);
public:
  TMap *fHisMap;      // map of the histogram per trigger class 
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}
private:
   AliTPCcalibTrigger(const AliTPCcalibTrigger&); // Not implemented
   AliTPCcalibTrigger& operator=(const AliTPCcalibTrigger&); // Not implemented
 
  ClassDef(AliTPCcalibTrigger, 1); 
};

#endif


