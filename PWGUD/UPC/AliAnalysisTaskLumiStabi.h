/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKLUMISTABI_H
const Int_t STARTRUN = 240000;
const Int_t ENDRUN = 300000;
#define ALIANALYSISTASKLUMISTABI_H

class TH1;
class TH2;
class TTree;
class TList;
class TFile;
class TBits;

#include "AliAnalysisTaskSE.h"
class AliMultSelection;

class AliAnalysisTaskLumiStabi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskLumiStabi();
  AliAnalysisTaskLumiStabi(const char *name);
  virtual ~AliAnalysisTaskLumiStabi();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
 private:
 

  TList *fOutputList;
   
  TTree *tOutput;

  TH1I *hTriggerClassesCounter;
  TH1I *hTriggerInputsCounter;
  TH1D *hCentralityV0M;
  Bool_t fTrgClassCINTZAC, fTrgInputV0M;
  Int_t fRunNumber;

  Float_t fCentralityPercentile = 300;
  AliMultSelection *fCentrality = 0x0;
  
  AliAnalysisTaskLumiStabi(const AliAnalysisTaskLumiStabi&); //not implemented
  AliAnalysisTaskLumiStabi& operator =(const AliAnalysisTaskLumiStabi&); //not implemented
  
  ClassDef(AliAnalysisTaskLumiStabi, 37);
};

#endif
