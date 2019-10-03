#ifndef ALIJJETCORRTASK_H
#define ALIJJETCORRTASK_H

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
#include "AliJJetCORRAnalysis.h"
#include "AliJJetTask.h"
#include "AliJCard.h"

class AliJJetCORRAnalysis;
class AliJRunTable;
//==============================================================

using namespace std;

class AliJJetCORRTask : public AliAnalysisTaskSE {

 public:
  AliJJetCORRTask();
  AliJJetCORRTask(const char *name,  TString inputformat);
  AliJJetCORRTask(const AliJJetCORRTask& ap);   
  AliJJetCORRTask& operator = (const AliJJetCORRTask& ap);
  virtual ~AliJJetCORRTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  virtual Bool_t UserNotify() { std::cout<<"DEBUG UserNotify"<<std::endl; return kTRUE;}

  void SetJetTaskName(TString name){ fJetTaskName=name; }
  void SetJJetCORRAnalysis(AliJJetCORRAnalysis * jco){ fJJetCORRAnalysis=jco; }
  void SetCard( AliJCard * card ){ fCard=card; }
  void SetTargetJetIndex( int jfindex ) { fTargetJetIndex = jfindex; } // PlayCorrelation only for selected jetfinder

  bool IsGoodEvent(AliVEvent *event);
  void SetDebugMode( int debug) { fDebugMode = debug; };

 private:
  AliJJetTask           * fJetTask;
  TString                 fJetTaskName;
  AliJJetCORRAnalysis     * fJJetCORRAnalysis;	//!
  TDirectory     * fOutput;
  AliJCard              * fCard;
  AliAnalysisUtils *fAnaUtils;
  AliJRunTable *fRunTable; //
  Bool_t fFirstEvent; //
  int cBin;
  int zBin;
  double zVert;
  int fDebugMode;
  int fTargetJetIndex;
  int fevt;

  ClassDef(AliJJetCORRTask, 1); 
};
#endif // AliJJetCORRTask_H
