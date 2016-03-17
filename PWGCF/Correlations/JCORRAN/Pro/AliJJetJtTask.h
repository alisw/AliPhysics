#ifndef ALIJJETJTTASK_H
#define ALIJJETJTTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////


#include "AliAnalysisTaskSE.h"
#include "AliJJetJtAnalysis.h"
#include "AliAnalysisUtils.h"
#include "AliJJetTask.h"
#include "AliJCard.h"
#include "AliJRunTable.h"

class AliJJetJtAnalysis;
class AliJEfficiency;
//==============================================================

using namespace std;

class AliJJetJtTask : public AliAnalysisTaskSE {

 public:
  AliJJetJtTask();
  AliJJetJtTask(const char *name,  TString inputformat);
  AliJJetJtTask(const AliJJetJtTask& ap);   
  AliJJetJtTask& operator = (const AliJJetJtTask& ap);
  virtual ~AliJJetJtTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
          bool IsGoodEvent(AliVEvent *event);
  virtual void Terminate(Option_t *);
  virtual void FinishTaskOutput();
  virtual Bool_t UserNotify() { std::cout<<"DEBUG UserNotify"<<std::endl; return kTRUE;}

  void SetJetTaskName(TString name){ fJetTaskName=name; }
  void SetMCJetTaskName(TString name){ fMCJetTaskName=name; }
  void SetJJetJtAnalysis(AliJJetJtAnalysis * jco){ fJJetJtAnalysis=jco; }
  void SetCard( AliJCard * card ){ fCard=card; }
  void SetMC(int mc) {fDoMC = mc;};
  void SetNrandom( int Nrand) { NRandom = Nrand;}
  void SetMoveJet( int move) { moveJet = move;}

 private:
  
  // TODO new Task - AliJJetTask?
  AliJJetTask           * fJetTask;
  AliJJetTask           * fMCJetTask;
  TString                 fJetTaskName;
  TString                 fMCJetTaskName;
  AliJJetJtAnalysis     * fJJetJtAnalysis;
  TDirectory     * fOutput;
  AliJCard              * fCard;
  Bool_t fFirstEvent;
  int cBin;
  int zBin;
  int NRandom;
  int moveJet;
  int fDoMC;
  double zVert;
  AliAnalysisUtils *fAnaUtils;
  AliJRunTable *fRunTable;
  TH1D * fEventHist;

  ClassDef(AliJJetJtTask, 1); 
};
#endif // ALIJJETJTTASK_H
