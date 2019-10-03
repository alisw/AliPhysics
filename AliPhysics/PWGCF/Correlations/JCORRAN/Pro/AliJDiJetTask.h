#ifndef ALIJDIJETTASK_H
#define ALIJDIJETTASK_H

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
#include "AliJDiJetAnalysis.h"
#include "AliJJetTask.h"
#include "AliJCard.h"

class AliJRunTable;
class AliJDiJetAnalysisTask;
//==============================================================

using namespace std;

class AliJDiJetTask : public AliAnalysisTaskSE {

public:
  AliJDiJetTask();
  AliJDiJetTask(const char *name,  TString inputformat);
  AliJDiJetTask(const AliJDiJetTask& ap);   
  AliJDiJetTask& operator = (const AliJDiJetTask& ap);
  virtual ~AliJDiJetTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  virtual Bool_t UserNotify() { return kTRUE;}

  bool IsGoodEvent( AliVEvent * event );

  void SetJetTaskName(TString name){ fJetTaskName=name; }
  void SetDiJetAnalysis(AliJDiJetAnalysis * jco){ fJDiJetAnalysis=jco; }
  void SetCard( AliJCard * card ){ fCard=card; }

private:

  // TODO new Task - AliJJetTask?
  AliJJetTask           * fJetTask;
  TString                 fJetTaskName;
  AliJDiJetAnalysis     * fJDiJetAnalysis;	//!
  TDirectory     * fOutput;
  Bool_t fFirstEvent; //
  AliAnalysisUtils *fAnaUtils;
  AliJRunTable *fRunTable; //
  AliJCard              * fCard;

  ClassDef(AliJDiJetTask, 1); 
};
#endif // AliJDiJetTask_H
