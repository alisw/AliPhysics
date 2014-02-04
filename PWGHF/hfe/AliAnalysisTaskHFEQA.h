/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Flow task class for the ALICE HFE group
//
//
#ifndef ALIANALYSISTASKHFEQA_H
#define ALIANALYSISTASKHFEQA_H


#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#ifndef ROOT_TString
#include <TString.h>
#endif

#ifndef ROOT_TBits
#include <TBits.h>
#endif


class TList;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliAODMCHeader;
class TClonesArray;
class TH1F;


class AliAnalysisTaskHFEQA: public AliAnalysisTaskSE {
public:

  typedef enum{
    kpp = 0,
    kpPb = 1,
    kPbPb = 2
  } ESystem_t;
  

  AliAnalysisTaskHFEQA();
  AliAnalysisTaskHFEQA(const char *name);
  AliAnalysisTaskHFEQA(const AliAnalysisTaskHFEQA &ref);
  AliAnalysisTaskHFEQA& operator=(const AliAnalysisTaskHFEQA &ref);
  virtual void Copy(TObject &o) const;
  virtual ~AliAnalysisTaskHFEQA();
  
  virtual void  UserExec(Option_t */*option*/);
  virtual void  UserCreateOutputObjects();

  void SetDoTPConly(Bool_t tpconlydo)   { fTPConlydo = tpconlydo; };
  void SetDoTRDonly(Bool_t trdonlydo)   { fTRDonlydo = trdonlydo; };
  void SetDoTOFTPC(Bool_t toftpcdo)   { fTOFTPCdo = toftpcdo; };
  void SetDoTPCTRD(Bool_t tpctrddo)   { fTPCTRDdo = tpctrddo; };
  void SetDoTPCEMCal(Bool_t tpcemcaldo)   { fTPCEMCaldo = tpcemcaldo; };
  void SetAODAnalysis(Bool_t aodAnalysis)   { fAODAnalysis = aodAnalysis; };
  void SetCentralityEstimator(const char *estimator) { fCentralityEstimator = estimator; }
  void SetppAnalysis(){
    fCollisionSystem.SetBitNumber(kpPb, kFALSE); 
    fCollisionSystem.SetBitNumber(kPbPb, kFALSE); 
    fCollisionSystem.SetBitNumber(kpp, kTRUE); 
  }
  void SetpPbAnalysis() {
    fCollisionSystem.SetBitNumber(kpp, kFALSE); 
    fCollisionSystem.SetBitNumber(kPbPb, kFALSE); 
    fCollisionSystem.SetBitNumber(kpPb, kTRUE); 
  }
  void SetPbPbAnalysis() { 
    fCollisionSystem.SetBitNumber(kpp, kFALSE); 
    fCollisionSystem.SetBitNumber(kpPb, kFALSE); 
    fCollisionSystem.SetBitNumber(kPbPb, kTRUE); 
  };
  Bool_t Ispp() const { return fCollisionSystem.TestBitNumber(kpp); }
  Bool_t IsPbPb() const { return fCollisionSystem.TestBitNumber(kPbPb); }
  Bool_t IspPb() const { return fCollisionSystem.TestBitNumber(kpPb); }
  
  AliHFEpid *GetPIDTPConly() const { return fPIDTPConly; }
  AliHFEpid *GetPIDTRDonly() const { return fPIDTRDonly; }
  AliHFEpid *GetPIDTOFTPC() const { return fPIDTOFTPC; }
  AliHFEpid *GetPIDTPCTRD() const { return fPIDTPCTRD; }
  AliHFEpid *GetPIDTPCEMCal() const { return fPIDTPCEMCal; }
  AliHFEpidQAmanager *GetPIDQAManagerTRDonly() const { return fPIDqaTRDonly; }
  AliHFEpidQAmanager *GetPIDQAManagerTOFTPC() const { return fPIDqaTOFTPC; }
  AliHFEpidQAmanager *GetPIDQAManagerTPCTRD() const { return fPIDqaTPCTRD; }
  AliHFEpidQAmanager *GetPIDQAManagerTPCEMCal() const { return fPIDqaTPCEMCal; }
 

  void SetHFECuts(AliHFEcuts * const cuts) { fHFECuts = cuts; };
 
 private:
  TList     *fListHist;         //! TH list
  Bool_t    fAODAnalysis;       // AOD analysis
  AliAODMCHeader *fAODMCHeader;         // ! MC info AOD
  TClonesArray *fAODArrayMCInfo;        // ! MC info particle AOD
 
  // Cuts for HFE
  AliHFEcuts *fHFECuts;                   // HFE cuts
  AliHFEpid  *fPIDTPConly;                // PID cuts 
  AliHFEpid  *fPIDTRDonly;                // PID cuts 
  AliHFEpid  *fPIDTOFTPC;                 // PID cuts TOF-TPC only
  AliHFEpid  *fPIDTPCTRD;                 // PID cuts TPC-TRD 
  AliHFEpid  *fPIDTPCEMCal;                 // PID cuts TPC-EMCal 
  AliHFEpidQAmanager *fPIDqaTRDonly;       // QA Manager TOF TPC
  AliHFEpidQAmanager *fPIDqaTOFTPC;       // QA Manager TOF TPC
  AliHFEpidQAmanager *fPIDqaTPCTRD;       // QA Manager TPC TRD
  AliHFEpidQAmanager *fPIDqaTPCEMCal;       // QA Manager TPC EMCal
  TString fCentralityEstimator;         // Centrality Estimator
  TBits fCollisionSystem;              // Collision System;

  // Histo yields
  TH1F *fNbEvent;                      // Number of events
  TH1F *fTPConly;                      // TPC only electron yield
  TH1F *fTOFTPC;                      // TOF TPC  electron yield
  TH1F *fTPCTRD;                      // TPC TRD  electron yield
  TH1F *fTPCEMCal;                      // TPC EMCal electron yield

  // Do PID or not
  Bool_t fTPConlydo;                   // Do TPC only PID
  Bool_t fTRDonlydo;                   // Do TRD only PID
  Bool_t fTOFTPCdo;                    // Do TOF TPC 
  Bool_t fTPCTRDdo;                    // Do TPC TRD 
  Bool_t fTPCEMCaldo;                    // Do TPC EMCal 

  
  // Debuging Cuts step by step all centrality together: pt, step (6)
  //THnSparseF *fTrackingCuts; //! Tracking Cuts

  
  
  ClassDef(AliAnalysisTaskHFEQA, 3); // analysisclass
};

#endif
