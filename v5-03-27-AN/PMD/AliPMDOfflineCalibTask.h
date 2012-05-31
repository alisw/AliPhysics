#ifndef AliPMDOfflineCalibTask_cxx
#define AliPMDOfflineCalibTask_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/************************************
 *
 * Satyajit Jena, IIT Bombay
 * sjena@cern.ch
 * Fri Feb 12 13:30:19 IST 2010
 *
 ************************************/

class TH1F;
class TList;
class TObjArray;

#include "TObjString.h"
#include "AliAnalysisTaskSE.h"

class AliPMDOfflineCalibTask : public AliAnalysisTaskSE {
 public:
  AliPMDOfflineCalibTask(const char *name = "AliPMDOfflineCalibTask");
  virtual ~AliPMDOfflineCalibTask() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void AddSelectedTriggerClass(const char*name) {fSelectedTrigger->Add(new TObjString(name));};
  void SetReject(Bool_t rejected) {fRejected = rejected;};
  
 private:

  TList *fListOfHistos;
  
  TH1F *fPmdCalibAdcP; 
  TH1F *fPmdCalibAdcC; 
  TH1F *fPmdCalibEntP; 
  TH1F *fPmdCalibEntC; 

  TH1I *fNEvents; 
  
  TObjArray *fSelectedTrigger;   
  Bool_t fRejected;

  AliPMDOfflineCalibTask(const AliPMDOfflineCalibTask&); 
  AliPMDOfflineCalibTask& operator=(const AliPMDOfflineCalibTask&); 
  
  ClassDef(AliPMDOfflineCalibTask, 1); 
};

#endif

