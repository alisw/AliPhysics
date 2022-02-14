#ifndef ALIJJETQATASK_H
#define ALIJJETQATASK_H

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
#include "AliJJetQAAna.h"
#include "AliJJetTask.h"
#include "AliJCard.h"

class AliJJetQAAna;
class AliJRunTable;
//==============================================================

using namespace std;

class AliJJetQATask : public AliAnalysisTaskSE {

 public:
  AliJJetQATask();
  AliJJetQATask(const char *name,  TString inputformat);
  AliJJetQATask(const AliJJetQATask& ap);   
  AliJJetQATask& operator = (const AliJJetQATask& ap);
  virtual ~AliJJetQATask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  virtual Bool_t UserNotify() { std::cout<<"DEBUG UserNotify"<<std::endl; return kTRUE;}

  void SetJetTaskName(TString name){ fJetTaskName=name; }
  void SetJJetQAAna(AliJJetQAAna * jco){ fJJetQAAna=jco; }
  void SetCard( AliJCard * card ){ fCard=card; }
  void SetTargetJetIndex( int jfindex ) { fTargetJetIndex = jfindex; } // PlayCorrelation only for selected jetfinder

  bool IsGoodEvent(AliVEvent *event);
  void SetDebugMode( int debug) { fDebugMode = debug; };
  void AddFlags(UInt_t nflags){flags |= nflags;}
  enum{
    FLUC_MC = 0x1,
    FLUC_EXCLUDEWDECAY = 0x2,
    FLUC_KINEONLY = 0x4,
    FLUC_CUT_OUTLIERS = 0x200,
    FLUC_ALICE_IPINFO = 0x400,
  };

 private:
  AliJJetTask           * fJetTask;
  TString                 fJetTaskName;
  AliJJetQAAna     * fJJetQAAna;	//!
  TDirectory     * fOutput;
  AliJCard              * fCard;
  AliAnalysisUtils *fAnaUtils;
  AliJRunTable *fRunTable; //
  Bool_t fFirstEvent; //
  float fcent; 
  int cBin;
  int zBin;
  double zVert;
  int fDebugMode;
  int fTargetJetIndex;
  int fevt;
  UInt_t flags; //

  ClassDef(AliJJetQATask, 1); 
};
#endif // AliJJetQATask_H
