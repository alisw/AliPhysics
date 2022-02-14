#ifndef ALIANALYSISTASKPIDPERFORMCOMBPIDPTDEP_H
#define ALIANALYSISTASKPIDPERFORMCOMBPIDPTDEP_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#########################################################
//#                                                       # 
//#        Task for testing the combined PID              #
//#                                                       #
//#  Pietro Antonioli, INFN / Pietro.Antonioli@bo.infn.it #
//#  Jens Wiechula, Uni Tübingen / Jens.Wiechula@cern.ch  #
//#                                                       #
//#########################################################

#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliPIDCombined.h>

#include <AliESDtrackCuts.h>
#include <AliAnalysisFilter.h>

#include "AliAnalysisTaskSE.h"

class TH1F;
class TH2D;
class TH3D;
class fPIDResponse;
class fPIDCombined;
class TList;

class AliAnalysisTaskPIDPerformCombIDPtDep : public AliAnalysisTaskSE {
  
 public:
  static const Int_t kPtBins = 9;
  
  
  AliAnalysisTaskPIDPerformCombIDPtDep();
  AliAnalysisTaskPIDPerformCombIDPtDep(const char *name);
  virtual ~AliAnalysisTaskPIDPerformCombIDPtDep(){;}
  
  virtual void  UserExec(Option_t *option);
  virtual void  UserCreateOutputObjects();
  virtual void  Terminate(Option_t *option);
  
  void SetFilterBit(Int_t fb){fFB = fb;}
  Int_t GetMomBin(Float_t mom);

  //Centrality
  void UseCentrality() { fUseCentrality = kTRUE;}
  void SetCentralityEstimator(const char* centralityEstimator) {
    fCentralityEstimator = centralityEstimator;}
  void SetCentralityPercentileRange(Float_t min, Float_t max) { 
    fCentralityPercentileMin=min;
    fCentralityPercentileMax=max;
  }
  
  void SetBayesThreshold1(Float_t bth1) {fbayesth1 = bth1;}
  void SetBayesThreshold2(Float_t bth2) {fbayesth2 = bth2;}
  void SetBayesThreshold3(Float_t bth3) {fbayesth3 = bth3;}

  void SetPIDMomCut(Float_t pidMomCut)  {fPIDMomCut = pidMomCut;}
    
  void SetParticleOfInterest(AliPID::EParticleType partOfInterest) { fpartOfInterest = partOfInterest;}
    
  void SetRejectInjectedSignalsGenName(TString genToBeKept) {
        fGenToBeKept = genToBeKept;
        fRejectCheckGenName=kTRUE;
  }
  
  void SetExtraTPCCutsAOD(Double_t maxTPCchi2, Int_t minNClustersTPC, Int_t minNTPCCrossedRows, Float_t minNTPCFindableCls){
    fTPCchi2Cut      = maxTPCchi2;
    fNClustersTPCCut = minNClustersTPC;
    fMinTPCCrossedRows = minNTPCCrossedRows;
    fMinTPCRowsOverFindableCls =  minNTPCFindableCls;
  }
  
  
  private:

  TList *fHistList;                   //! list of histograms

  const AliPIDResponse *fPIDResponse;     //! PID response object
  AliPIDCombined       *fPIDCombined;     //! combined PID object

  AliVEvent *fAOD;

  Bool_t  fUseCentrality;// Bool_t use centrality or not
  TString fCentralityEstimator;//"V0M","TRK","TKL","ZDC","FMD"
  Float_t fCentralityPercentileMin, fCentralityPercentileMax; //min-max centrality percentile
  TH1F *fHistCentrality;
  
  TH1F *hTrue[AliPID::kSPECIES][2]; //!
  TH1F *hTrueInAccTPC[AliPID::kSPECIES][2];  //!
  TH1F *hTrueInAccTOF[AliPID::kSPECIES][2]; //!
  TH1F *hTrueInAccTPCTOF[AliPID::kSPECIES][2]; //!
  TH1F *hTrueInAccTPCTOFBayes[AliPID::kSPECIES][2]; //!
  
  TH1F *hIdTPConly2s[2]; //!
  TH1F *hIdTPConly3s[2]; //!

  TH1F *hIdTOFonly2s[2]; //!
  TH1F *hIdTOFonly3s[2]; //!

  TH1F *hIDnSigmaComb1[2]; //!
  TH1F *hIDnSigmaComb2[2]; //!
  TH1F *hIDnSigmaComb3[2]; //!

  TH1F *hIDBayes1[2]; //!
  TH1F *hIDBayes2[2]; //!
  TH1F *hIDBayes3[2]; //!
  
  TH1F *hEffPlotsTPConly2s[AliPID::kSPECIES][2]; //!
  TH1F *hEffPlotsTPConly3s[AliPID::kSPECIES][2]; //!
  TH1F *hEffPlotsTOFonly2s[AliPID::kSPECIES][2]; //!
  TH1F *hEffPlotsTOFonly3s[AliPID::kSPECIES][2]; //!
  
  TH1F* hEffPlotsnSigmaComb1[AliPID::kSPECIES][2]; //!
  TH1F* hEffPlotsnSigmaComb2[AliPID::kSPECIES][2]; //!
  TH1F* hEffPlotsnSigmaComb3[AliPID::kSPECIES][2]; //!

  TH1F* hEffPlotsBayes1[AliPID::kSPECIES][2]; //!
  TH1F* hEffPlotsBayes2[AliPID::kSPECIES][2]; //!
  TH1F* hEffPlotsBayes3[AliPID::kSPECIES][2]; //!
    

  TH2D *fPriorsUsed[AliPID::kSPECIES];//! priors used

  static const char* fgkBinMomDesc[kPtBins];

  AliPID::EParticleType fpartOfInterest; 
  
  Int_t fFB; //track filterbit
  Float_t fbayesth1;
  Float_t fbayesth2;
  Float_t fbayesth3;
    
  TString fGenToBeKept;
  Bool_t fRejectCheckGenName;
  Float_t fPIDMomCut;

  Double_t fTPCchi2Cut;//only used for AODs
  Int_t fNClustersTPCCut;//only used for AODs
  Int_t fMinTPCCrossedRows; //only used for AODs
  Float_t fMinTPCRowsOverFindableCls; //only used for AODs

  ClassDef(AliAnalysisTaskPIDPerformCombIDPtDep, 2);

};
#endif
