#ifndef ALITPCCALIBKRTASK_H
#define ALITPCCALIBKRTASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TList.h>

#include "AliTPCclusterKr.h"
#include "AliAnalysisTask.h"

#include "AliTPCCalibKr.h"

class AliTPCCalibKrTask : public AliAnalysisTask {

public:
  AliTPCCalibKrTask(const char *name = "AliTPCCalibKrTask");
  virtual ~AliTPCCalibKrTask();

  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual bool   Notify() { return kTRUE;}

  Bool_t ReadEntry(Int_t evt);
  void SetInputChain(TChain *inChain) {fTree = (TTree*) inChain;}
  void SetTPCCalibKr(AliTPCCalibKr *calibKr) {fTPCCalibKr = calibKr;}

private:

  AliTPCCalibKrTask(const AliTPCCalibKrTask&); // not implemented
  AliTPCCalibKrTask operator=(const AliTPCCalibKrTask&); // not implemented

  static Int_t  fEvtNumber;     //! event number
  AliTPCclusterKr *fClustKr;  //! input AliTPCclusterKr objects
  AliTPCCalibKr *fTPCCalibKr; // output AliTPCCalibKr objects

  TTree *fTree;               //! input tree
  TList *fOutput;             //! output list of objects

  ClassDef(AliTPCCalibKrTask, 1)  // TPC task
};

#endif

