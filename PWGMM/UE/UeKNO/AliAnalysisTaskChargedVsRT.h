/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskChargedVsRT_H
#define AliAnalysisTaskChargedVsRT_H

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

class AliAnalysisTaskChargedVsRT : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskChargedVsRT();
  AliAnalysisTaskChargedVsRT(const char *name);
  virtual ~AliAnalysisTaskChargedVsRT();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void GetLeadingObjectFromArray(const std::vector<Float_t> &pt,
                                 const std::vector<Float_t> &phi,
                                 const std::vector<Int_t> &id, Int_t multPart,
                                 Bool_t isMC);
  void GetDetectorResponse(const std::vector<Float_t> &phiGen, Int_t multGen,
                           const std::vector<Float_t> &phiRec, Int_t multRec,
                           const std::vector<Int_t> &idGen);
  void GetBinByBinCorrections(Int_t multGen, Int_t multRec,
                              const std::vector<Float_t> &ptGen,
                              const std::vector<Float_t> &ptRec,
                              const std::vector<Int_t> &idGen,
                              const std::vector<Int_t> &idRec,
                              const std::vector<Int_t> &isprimRec);
  void GetMultiplicityDistributionsTrue(const std::vector<Float_t> &phiGen,
                                        const std::vector<Float_t> &ptGen,
                                        Int_t multGen,
                                        const std::vector<Int_t> &idGen);

  void GetMultiplicityDistributions(const std::vector<Float_t> &phiRec,
                                    const std::vector<Float_t> &ptRec,
                                    Int_t multRec,
                                    const std::vector<Float_t> &ptRecwodca,
                                    const std::vector<Float_t> &dcaxyRecwodca,
                                    const std::vector<Int_t> &isprimwodca,
                                    Int_t multRecwodca);

  void GetMultiplicityDistributionsData(
      const std::vector<Float_t> &phiRec, const std::vector<Float_t> &ptRec,
      Int_t multRec, const std::vector<Float_t> &ptRecwodca,
      const std::vector<Float_t> &dcaxyRecwodca, Int_t multRecwodca);
  void SetPtMin(Double_t val) {
    fPtMin = val;
  } // Set pT cut for associated particles
  void SetLeadingPtMin(Double_t PtLmin) {
    fLeadPtCutMin = PtLmin;
  } // use differnet ptcuts
  void SetLeadingPtMax(Double_t PtLmax) {
    fLeadPtCutMax = PtLmax;
  } // use differnet ptcuts
  void SetNchNbin(Int_t NchNbins) {
    fNchNbin = NchNbins;
  } // use different bining
  void SetNchBinMax(Double_t maxbinNch) {
    fNchBinMax = maxbinNch;
  } // use different bin max

  void SetUseMC(Bool_t mc = kFALSE) { fUseMC = mc; } // use to analyse MC data
  void SetMCclosureTest(Bool_t mcc = kFALSE) { fIsMCclosure = mcc; }
  void SetIsHybridAnalysis(Bool_t isHy = kFALSE) { fIsHybAna = isHy; }
  void SetMultPercenV0(Bool_t multV0 = kFALSE) { fMultPercenV0 = multV0; }
  bool HasRecVertex();
  void SetCutsHybrid0WoDCA(AliESDtrackCuts *name);
  void SetCutsHybrid1WoDCA(AliESDtrackCuts *name);
  void SetCutsFilterWoDCA(AliESDtrackCuts *name);
  // Systematic ============================
  void SetTPCclustersVar1(Bool_t TPCclustersVar1 = kFALSE) {
    fTPCclustersVar1 = TPCclustersVar1;
  }
  void SetTPCclustersVar2(Bool_t TPCclustersVar2 = kFALSE) {
    fTPCclustersVar2 = TPCclustersVar2;
  }
  void SetNcrVar1(Bool_t NcrVar1 = kFALSE) { fNcrVar1 = NcrVar1; }
  void SetNcrVar2(Bool_t NcrVar2 = kFALSE) { fNcrVar2 = NcrVar2; }
  void SetChisqTPCVar1(Bool_t ChisqTPCVar1 = kFALSE) {
    fChisqTPCVar1 = ChisqTPCVar1;
  }
  void SetChisqTPCVar2(Bool_t ChisqTPCVar2 = kFALSE) {
    fChisqTPCVar2 = ChisqTPCVar2;
  }
  void SetChisqITSVar1(Bool_t ChisqITSVar1 = kFALSE) {
    fChisqITSVar1 = ChisqITSVar1;
  }
  void SetChisqITSVar2(Bool_t ChisqITSVar2 = kFALSE) {
    fChisqITSVar2 = ChisqITSVar2;
  }
  //         void       SetChisqITSmTPCVar1(Bool_t ChisqITSmTPCVar1 = kFALSE)
  //         {fChisqITSmTPCVar1 = ChisqITSmTPCVar1;} void
  //         SetChisqITSmTPCVar2(Bool_t ChisqITSmTPCVar2 = kFALSE)
  //         {fChisqITSmTPCVar2 = ChisqITSmTPCVar2;}
  void SetDcazVar1(Bool_t DcazVar1 = kFALSE) { fDcazVar1 = DcazVar1; }
  void SetDcazVar2(Bool_t DcazVar2 = kFALSE) { fDcazVar2 = DcazVar2; }
  void SetGeoTPCVar1(Bool_t GeoTPCVar1 = kFALSE) { fGeoTPCVar1 = GeoTPCVar1; }
  void SetGeoTPCVar2(Bool_t GeoTPCVar2 = kFALSE) { fGeoTPCVar2 = GeoTPCVar2; }
  void SetGeoTPCVar3(Bool_t GeoTPCVar3 = kFALSE) { fGeoTPCVar3 = GeoTPCVar3; }
  void SetGeoTPCVar4(Bool_t GeoTPCVar4 = kFALSE) { fGeoTPCVar4 = GeoTPCVar4; }
  //         void       SetSPDreqVar1(Bool_t SPDreqVar1 = kFALSE) {fSPDreqVar1 =
  //         SPDreqVar1;}
  // Systematic ============================
  virtual Double_t DeltaPhi(Double_t phia, Double_t phib,
                            Double_t rangeMin = -TMath::Pi() / 2,
                            Double_t rangeMax = 3 * TMath::Pi() / 2);
  Int_t FillArrayMC(std::vector<Float_t> &pt, std::vector<Float_t> &phi,
                    std::vector<Int_t> &id);
  Int_t FillArray(std::vector<Float_t> &pt, std::vector<Float_t> &phi,
                  std::vector<Float_t> &dcaxy, std::vector<Float_t> &dcaz,
                  std::vector<Int_t> &isprim, std::vector<Int_t> &id,
                  const bool wDCA, const bool useHy);

protected:
private:
  AliESDEvent *fESD; //! input ESD event
  AliEventCuts fEventCuts;
  AliStack *fMCStack; //! MC stack
  AliMCEvent *fMC;    //! MC Event
  Bool_t fUseMC;      // analyze MC events
  Bool_t fIsMCclosure;
  Bool_t fIsHybAna;
  Bool_t fMultPercenV0;
  Int_t fnRecHy;
  Int_t fnRecHyWoDCA;
  Int_t fnGen;
  // Systematic------------------------------------
  Bool_t fNcrVar1;
  Bool_t fNcrVar2;
  Bool_t fTPCclustersVar1;
  Bool_t fTPCclustersVar2;
  Bool_t fGeoTPCVar1;
  Bool_t fGeoTPCVar2;
  Bool_t fGeoTPCVar3;
  Bool_t fGeoTPCVar4;
  Bool_t fChisqTPCVar1;
  Bool_t fChisqTPCVar2;
  Bool_t fChisqITSVar1;
  Bool_t fChisqITSVar2;
  Bool_t fDcazVar1;
  Bool_t fDcazVar2;
  // Systematic------------------------------------
  AliPIDResponse *fPIDResponse;
  AliAnalysisFilter *fTrackFilter;
  AliAnalysisFilter *fTrackFilterwoDCA;
  AliAnalysisFilter *fTrackFilterHybrid0;
  AliAnalysisFilter *fTrackFilterHybrid1;
  AliAnalysisFilter *fTrackFilterHybrid0woDCA;
  AliAnalysisFilter *fTrackFilterHybrid1woDCA;
  TList *fOutputList; //! output list in the root file
  Double_t fEtaCut;
  Double_t fPtMin;
  Double_t fLeadPtCutMin;
  Double_t fLeadPtCutMax;
  Int_t fNchNbin;
  Double_t fNchBinMax;
  Double_t fGenLeadPhi;
  Double_t fGenLeadPt;
  Int_t fGenLeadIn;
  Double_t fRecLeadPhi;
  Double_t fRecLeadPt;
  Int_t fRecLeadIn;
  Double_t ftrackmult08;
  Double_t fv0mpercentile;
  Float_t fdcaxy;
  Float_t fdcaz;
  AliMultSelection *fMultSelection;
  AliMultSelection *fMultSelectionbefvtx;
  // KNO
  TH1D *hNchTSGen;
  TH1D *hNchTSGenTest;
  TH1D *hNchGen;
  TH1D *hNchGenTest;
  TH1D *hNchTSRec;
  TH1D *hNchTSRecTest;
  TH1D *hNchData;
  TH1D *hNchTSData;
  TH1F *hPhiTotal;
  TH1F *hPhiStandard;
  TH1F *hPhiHybrid1;
  TH2D *hNchResponse;
  TH1D *hNchRec;
  TH1D *hNchRecTest;
  TH1D *hPtInPrim;
  TH1D *hPtInPrim_lambda;
  TH1D *hPtInPrim_pion;
  TH1D *hPtInPrim_kaon;
  TH1D *hPtInPrim_proton;
  TH1D *hPtInPrim_sigmap;
  TH1D *hPtInPrim_sigmam;
  TH1D *hPtInPrim_omega;
  TH1D *hPtInPrim_xi;
  TH1D *hPtInPrim_rest;
  TH1D *hPtOut;
  TH1D *hPtOutPrim;
  TH1D *hPtOutPrim_lambda;
  TH1D *hPtOutPrim_pion;
  TH1D *hPtOutPrim_kaon;
  TH1D *hPtOutPrim_proton;
  TH1D *hPtOutPrim_sigmap;
  TH1D *hPtOutPrim_sigmam;
  TH1D *hPtOutPrim_omega;
  TH1D *hPtOutPrim_xi;
  TH1D *hPtOutPrim_rest;
  TH1D *hPtOutSec;
  TH1D *hCounter;
  TH2D *hPtVsUEGenTest[3];
  TH2D *hPtVsUERecTest[3];
  TH2D *hPtVsUEData[3];
  TH2D *hPtVsNchGenTest[3];
  TH2D *hPtVsNchRecTest[3];
  TH2D *hPtVsNchData[3];
  TH1D *hPhiGen[3];
  TH1D *hPhiRec[3];
  TH2D *hPTVsDCAData;
  TH2F *hptvsdcaPrim;
  TH2F *hptvsdcaDecs;
  TH2F *hptvsdcaMatl;
  TH2F *hptvsdcaAll;
  TProfile *hV0MvsNchT;

  AliAnalysisTaskChargedVsRT(
      const AliAnalysisTaskChargedVsRT &); // not implemented
  AliAnalysisTaskChargedVsRT &
  operator=(const AliAnalysisTaskChargedVsRT &); // not implemented

  ClassDef(AliAnalysisTaskChargedVsRT, 3);
};

#endif
