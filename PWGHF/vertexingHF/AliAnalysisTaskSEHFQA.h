#ifndef ALIANALYSISTASKSEHFQUALITYASSURANCE_H
#define ALIANALYSISTASKSEHFQUALITYASSURANCE_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSEHFQA
// AliAnalysisTaskSE for HF quality assurance
// Authors: C.Bianchin, chiara.bianchin@pd.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>

#include "AliAnalysisTaskSE.h"

class AliRDHFCuts;
class TH1F;

class AliAnalysisTaskSEHFQA : public AliAnalysisTaskSE
{

 public:

  enum DecChannel {kDplustoKpipi,kD0toKpi,kDstartoKpipi,kDstoKKpi,kD0toKpipipi,kLambdactopKpi};

  AliAnalysisTaskSEHFQA();
  AliAnalysisTaskSEHFQA(const char *name, DecChannel ch, AliRDHFCuts* cuts);
  virtual ~AliAnalysisTaskSEHFQA();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  //setters
  void SetReadMC(Bool_t mcflag){fReadMC=mcflag;}
  void SetSimpleMode(Bool_t flag){fSimpleMode=flag;}
  void SetTrackOn(Bool_t trackon=kTRUE){fOnOff[0]=trackon;}
  void SetPIDOn(Bool_t pidon=kTRUE){fOnOff[1]=pidon;}
  void SetCentralityOn(Bool_t centron=kTRUE){fOnOff[2]=centron;}
  void SetEvSelectionOn(Bool_t evselon=kTRUE){fOnOff[3]=evselon;}

  //getters
  AliRDHFCuts* GetCutObject() const {return fCuts;}
  DecChannel GetDecayChannel()const {return fDecayChannel;}
  Bool_t GetTrackStatus() const {return fOnOff[0];}
  Bool_t GetPIDStatus() const {return fOnOff[1];}
  Bool_t GetCentralityStatus() const {return fOnOff[2];}
  Bool_t GetEvSelStatus() const {return fOnOff[3];}

 private:
  AliAnalysisTaskSEHFQA(const AliAnalysisTaskSEHFQA &source);
  AliAnalysisTaskSEHFQA operator=(const AliAnalysisTaskSEHFQA &source);

 TH1F*  fNEntries;         //! histogram with number of events on output slot 1
 TList* fOutputPID;        //! list sent on output slot 2
 TList* fOutputTrack;      //! list sent on output slot 3
 TList* fOutputCounters;   //! list sent on output slot 5
 TList* fOutputCheckCentrality;   //! list sent on output slot 6
 TList* fOutputEvSelection; //! list sent on output slot 7
 DecChannel fDecayChannel; //identify the decay channel
 AliRDHFCuts* fCuts;       // object containing cuts
 AliRDHFCuts::ECentrality fEstimator; //2nd estimator for centrality
 Bool_t fReadMC;           // flag to read MC
 Bool_t fSimpleMode;       // if true, don't do candidates (much faster in PbPb)
 Bool_t fOnOff[4];         // on-off the QA on tracks (0), PID (1), centrality (2), event selection -- default is {kTRUE,kTRUE,kTRUE,kTRUE}
 ClassDef(AliAnalysisTaskSEHFQA,6); //AnalysisTaskSE for the quality assurance of HF in hadrons

};

#endif

