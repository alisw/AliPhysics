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
//class AliESDEvent;
//class AliESDtrack;
//class AliESDfriend;
class AliVEvent;
class AliVTrack;
class AliVfriendEvent;
class AliTPCseed;

class AliTPCAnalysisTaskcalib:public AliAnalysisTask {
public:
  AliTPCAnalysisTaskcalib();
  AliTPCAnalysisTaskcalib(const char *name);
  virtual ~AliTPCAnalysisTaskcalib();
  void AddJob(AliTPCcalibBase *job) {fCalibJobs->Add(job);}
  TObjArray* GetJobs() {return fCalibJobs;}

  virtual void ConnectInputData(Option_t *option);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void FinishTaskOutput();
  void         SetDebugOuputhPath(const char * name){fDebugOutputPath=name;}
protected:
  virtual void     Process(AliVEvent *event);
  virtual void     Process(AliTPCseed *track);
  virtual void     Process(AliVTrack *track, Int_t run);
  virtual Long64_t Merge(TCollection *li);
  virtual void     Analyze();
  void             RegisterDebugOutput();
private:
  TObjArray *fCalibJobs;      // array of calibration objects - WE ARE NOT OWNER?
  AliVEvent *fV;         //! current esd
  AliVfriendEvent *fVfriend;  //! current esd friend
  TString      fDebugOutputPath; // debug output path   
  AliTPCAnalysisTaskcalib(const AliTPCAnalysisTaskcalib&);
  AliTPCAnalysisTaskcalib& operator=(const AliTPCAnalysisTaskcalib&);
  ClassDef(AliTPCAnalysisTaskcalib,1)
};

#endif
