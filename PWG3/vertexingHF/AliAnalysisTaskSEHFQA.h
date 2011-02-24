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

  //getters
  AliRDHFCuts* GetCutObject() const {return fCuts;}
  DecChannel GetDecayChannel()const {return fDecayChannel;}

 private:
  AliAnalysisTaskSEHFQA(const AliAnalysisTaskSEHFQA &source);
  AliAnalysisTaskSEHFQA operator=(const AliAnalysisTaskSEHFQA &source);

 TH1F*  fNEntries;         //! histogram with number of events on output slot 1
 TList* fOutputPID;        //! list sent on output slot 2
 TList* fOutputTrack;      //! list sent on output slot 3
 TList* fOutputCounters;   //! list sent on output slot 5
 TList* fOutputCheckCentrality;   //! list sent on output slot 6
 DecChannel fDecayChannel; //identify the decay channel
 AliRDHFCuts* fCuts;       // object containing cuts
 AliRDHFCuts::ECentrality fEstimator; //2nd estimator for centrality
 Bool_t fReadMC;           // flag to read MC
 Bool_t fSimpleMode;       // if true, don't do candidates (much faster in PbPb) 
 ClassDef(AliAnalysisTaskSEHFQA,4); //AnalysisTaskSE for the quality assurance of HF in hadrons

};

#endif
