/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKLUMISTABI_H
#define ALIANALYSISTASKLUMISTABI_H

class TH1D;
class TH1I;
class TTree;
class TList;
class TFile;
class AliMultSelection;
class AliAnalysisCuts;
class AliAODZDC;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskLumiStabi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskLumiStabi();
  AliAnalysisTaskLumiStabi(const char *name);
  virtual ~AliAnalysisTaskLumiStabi();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  Bool_t IsSatellite(AliAODZDC *ZDCdata);

 private:
  TList *fOutputList;
   
  TTree *tOutput;

  TH1I *hDummyCounter;
  TH1I *hTriggerClassesCounter;
  TH1I *hTriggerInputsCounter;
  TH1D *hCentralityV0M;
  TH1D *hCentralityV0MandPS;
  TH1D *hCentralityV0MandSat;
  Bool_t fTrgClassCINTZAC, fTrgInputV0M, fSelectPhysics, fIsSatellite;
  Int_t fRunNumber;
  UInt_t fL0inputs;
  Float_t fV0McentPercentile = 300;
  Float_t fZNATDCm[4], fZNCTDCm[4];
//  AliMultSelection *fCentrality = 0x0;
  
  AliAnalysisTaskLumiStabi(const AliAnalysisTaskLumiStabi&); //not implemented
  AliAnalysisTaskLumiStabi& operator =(const AliAnalysisTaskLumiStabi&); //not implemented
  
  ClassDef(AliAnalysisTaskLumiStabi, 1);
};

#endif
