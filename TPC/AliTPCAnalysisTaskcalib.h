#ifndef ALIANALYSISTASKTPCCALIB_H
#define ALIANALYSISTASKTPCCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

#include "AliAnalysisTask.h"
#include "TObjArray.h"
#include "AliTPCcalibBase.h"
class AliESDEvent;
class AliESDfriend;
class AliTPCseed;

class AliTPCAnalysisTaskcalib:public AliAnalysisTask {
public:
  AliTPCAnalysisTaskcalib(const char *name);
  virtual ~AliTPCAnalysisTaskcalib();
  void AddJob(AliTPCcalibBase *job) {fCalibJobs.Add(job);}
  TObjArray* GetJobs() {return &fCalibJobs;}

  virtual void ConnectInputData(Option_t *option);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
protected:
  virtual void     Process(AliESDEvent *event);
  virtual void     Process(AliTPCseed *track);
  virtual Long64_t Merge(TCollection *li);
  virtual void     Analyze();
private:
  TObjArray fCalibJobs;
  AliESDEvent *fESD;
  AliESDfriend *fESDfriend;
  AliTPCAnalysisTaskcalib(const AliTPCAnalysisTaskcalib&);
  AliTPCAnalysisTaskcalib& operator=(const AliTPCAnalysisTaskcalib&);
  ClassDef(AliTPCAnalysisTaskcalib,1)
};

#endif
