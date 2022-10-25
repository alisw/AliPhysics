/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskSpectraFlatenicity_H
#define AliAnalysisTaskSpectraFlatenicity_H

class AliESDtrackCuts;
class AliESDEvent;
class TList;
class TH1D;
class TH2D;
class TH1I;
class TProfile;
class THnSparse;

#include "AliAnalysisTaskSE.h"

#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"

class AliAnalysisTaskSpectraFlatenicity : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSpectraFlatenicity();
  AliAnalysisTaskSpectraFlatenicity(const char *name);
  virtual ~AliAnalysisTaskSpectraFlatenicity();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void  SetPeriod(const char* dataset) { fDataSet = dataset; }  

  Double_t GetFlatenicityTPC();
  Double_t GetFlatenicityTPCMC();
  Double_t GetFlatenicity();
  Double_t GetFlatenicityV0A();
  Double_t GetFlatenicityV0C();
  Double_t GetFlatenicityMC();
  Double_t GetFlatenicityV0AMC();
  Double_t GetFlatenicityV0CMC();
  Double_t GetFlatenicityCombinedMC();
  Double_t GetFlatenicityCombined();
  void CheckMultiplicities();
  void CheckMultiplicitiesMC();
  void MakeMCanalysis();
  void MakeDataanalysis();

  void SetPtMin(Double_t val) {
    fPtMin = val;
  } // Set pT cut for associated particles
  void SetUseMC(Bool_t flat_flag = kFALSE) {
    fUseMC = flat_flag;
  } // use to analyse MC data
  void SetMCclosureTest(Bool_t flat_flag = kFALSE) { fIsMCclosure = flat_flag; }
  void SetDetectorForFlatenicity(TString det = "V0") { fDetFlat = det; }
  void SetRemoveTrivialScaling(Bool_t flat_flag = kFALSE) { fRemoveTrivialScaling = flat_flag; }
  void SetUseCalibration(Bool_t flat_flag = kFALSE) { fUseCalib = flat_flag; }  
  bool HasRecVertex();

protected:
    //
private:
    
  Double_t V0AmplCalibration(const Int_t &chnl);    
  Double_t V0AmplCalibrationTruth(const Int_t &chnl);    
    
  AliESDEvent *fESD; //! input ESD event
  AliEventCuts fEventCuts;
  AliStack *fMCStack; //! MC stack
  AliMCEvent *fMC;    //! MC Event
  Bool_t fUseMC;      // analyze MC events
  Int_t fV0Mindex;
  Float_t fmultV0A;
  Float_t fmultV0C;
  Float_t fmultTPC;
  Float_t fmultV0Amc;
  Float_t fmultV0Cmc;
  Float_t fmultTPCmc;
  TString fDetFlat;
  Bool_t fIsMCclosure;
  Bool_t fRemoveTrivialScaling;
  Int_t fnGen;
  AliPIDResponse *fPIDResponse;
  AliAnalysisFilter *fTrackFilter;
  TList *fOutputList; //! output list in the root file
  Double_t fEtaCut;
  Double_t fPtMin;
  Double_t ftrackmult08;
  Double_t fv0mpercentile;
  Float_t fFlat;
  Float_t fFlatV0ATPC;
  Float_t fFlatV0CTPC;
  Float_t fFlatV0TPC1;
  Float_t fFlatV0TPC2;
  Float_t fFlatV0AV0CTPC;
  Float_t fFlatV0AV0CTPC2;
  Float_t fFlatMC;
  Float_t fFlatMC2;
  Float_t fFlatMC3;
  Float_t fFlatMC4;
  Float_t fFlatMC5;
  Float_t fFlatMC6;
  AliMultSelection *fMultSelection;
  TH1D *hPtPrimIn;
  TH1D *hPtPrimOut;
  TH1D *hPtSecOut;
  TH1D *hPtOut;
  TH2D *hFlatV0vsFlatTPC;
  TH2D *hFlatV0AvsFlatTPC;
  TH2D *hFlatV0CvsFlatTPC;
  TH2D *hFlatV0vsFlatTPCmc;
  TH2D *hFlatV0AvsFlatTPCmc;
  TH2D *hFlatV0CvsFlatTPCmc;
  TH1D *hFlatenicity;
  TH1D *hFlatenicityMC;
  TH2D *hFlatCominedMC;
  TH2D *hFlat2CominedMC;
  TH2D *hFlat3CominedMC;
  TH2D *hFlat4CominedMC;
  TH2D *hFlat5CominedMC;
  TH2D *hFlat6CominedMC;
  TH2D *hFlatComined;
  TH2D *hFlatComined2;
  TH2D *hFlatComined3;
  TH2D *hFlatComined4;
  TH2D *hFlatComined5;
  TH2D *hFlatComined6;
  TH2D *hFlatResponse;
  TH2D *hFlatVsPt;
  TH2D *hFlatVsPtMC;
  Bool_t fUseCalib;
  TString fDataSet;
  TF1 *fV0Camp;
  TF1 *fV0Aamp;
  TProfile *pActivityV0DataSect;
  TProfile *pActivityV0ADataSect;
  TProfile *pActivityV0CDataSect;
  TProfile *pActivityV0multData;
  TProfile *pActivityV0AmultData;
  TProfile *pActivityV0CmultData;
  TProfile *pActivityV0McSect;
  TProfile *pActivityV0multMc;
  TH2D *hFlatVsNch;
  TH2D *hFlatVsNchTPC;
  TH2D *hFlatVsNchMC;
  TH2D *hFlatVsNchTPCmc;
  TH2D *hFlatVsNchCombinedMC;
  TH1D *hNchV0M;
  TH1D *hNchV0MMC;
  TH1D *hNchTPC;
  TH1D *hNchTPCmc;
  TH1D *hNchMidRapMC;
  TH1D *hNchV0a;
  TH1D *hNchV0c;
  TH1D *hNchCombinedmc;
  TH1D *hNchV0aMC;
  TH1D *hNchV0cMC;
  TH2D *hFlatVsV0M;
  TH2D* hFlatVsV0ATPC;
  TH2D* hFlatVsV0CTPC;
  TH2D* hFlatVsV0TPC1;
  TH2D* hFlatVsV0TPC2;
  TH2D* hFlatVsV0AV0CTPC;
  TH2D* hFlatVsV0AV0CTPC2;
  TH1D *hEta;
  TH1D *hEtamc;
  TH1D *hCounter;
  TH2D *hFlatVsPtV0M[9];
  TH2D *hFlatVsNchTPCV0M[9];
  TH2D *hFlatVsPtV0MMC[9];
  TH2D *hFlatVsNchTPCV0MMC[9];
  
  AliAnalysisTaskSpectraFlatenicity(const AliAnalysisTaskSpectraFlatenicity &); // not implemented
  AliAnalysisTaskSpectraFlatenicity &operator=(const AliAnalysisTaskSpectraFlatenicity &); // not implemented

  ClassDef(AliAnalysisTaskSpectraFlatenicity, 3);
};

#endif
