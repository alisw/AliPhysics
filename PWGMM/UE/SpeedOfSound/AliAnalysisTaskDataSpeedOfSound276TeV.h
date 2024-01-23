#ifndef ALIANALYSISTASKDATASPEEDOFSOUND276TEV_H
#define ALIANALYSISTASKDATASPEEDOFSOUND276TEV_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id$ */

// ROOT includes
#include <TH1.h>
#include <TList.h>
#include <TObject.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TTreeStream.h>

// AliRoot includes
#include <AliAODEvent.h>
#include <AliAODMCParticle.h>
#include <AliAnalysisFilter.h>
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliESDtrackCuts.h>
#include <AliGenEventHeader.h>
#include <AliMCEvent.h>
#include <AliStack.h>
#include <AliVHeader.h>

#include <cmath>

class AliAnalysisTaskDataSpeedOfSound276TeV : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDataSpeedOfSound276TeV();
  AliAnalysisTaskDataSpeedOfSound276TeV(const char* name);
  virtual ~AliAnalysisTaskDataSpeedOfSound276TeV();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);

  Bool_t GetAnalysisMC() { return fAnalysisMC; }
  Double_t GetVtxCut() { return fVtxCut; }
  Double_t GetEtaCut() { return fEtaCut; }

  virtual void SetTrigger(UInt_t ktriggerInt) { ftrigBit = ktriggerInt; }
  virtual void SetCentralityEstimator(const char* centEst) {
    fCentEst = centEst;
  }
  virtual void SetAnalysisType(const char* analysisType) {
    fAnalysisType = analysisType;
  }
  virtual void SetAnalysisMC(Bool_t isMC) { fAnalysisMC = isMC; }
  virtual void SetVtxCut(Double_t vtxCut) { fVtxCut = vtxCut; }
  virtual void SetEtaCut(Double_t etaCut) { fEtaCut = etaCut; }
  void SetEtaMinCut(const double& etamin) { fEtaMin = etamin; }
  void SetEtaMaxCut(const double& etamax) { fEtaMax = etamax; }
  void SetEtaGappT(const double& eta) { fEtaGappT = eta; }
  void SetEtaGapNch(const double& etamin, const double& etamax) {
    fEtaGapNchMin = etamin;
    fEtaGapNchMax = etamax;
  }
  void SetPtMin(const double& ptmin) { fPtMin = ptmin; }

  virtual void SetPileUpRej(Bool_t isrej) { fPileUpRej = isrej; }
  virtual void SetMinCent(Float_t minvalc) { fMinCent = minvalc; }
  virtual void SetMaxCent(Float_t maxvalc) { fMaxCent = maxvalc; }
  virtual void SetStoreMcIn(Bool_t value) { fStoreMcIn = value; }
  virtual void SetAnalysisPbPb(Bool_t isanaPbPb) { fAnalysisPbPb = isanaPbPb; }

 private:
  virtual Float_t GetVertex(const AliVEvent* event) const;
  virtual void AnalyzeESD(AliESDEvent* esd);
  virtual void AnalyzeAOD(AliAODEvent* aod);
  Short_t GetPidCode(Int_t pdgCode) const;
  void ProcessMCTruthESD();
  void ProcessMCTruthAOD();

  Short_t GetPythiaEventProcessType(Int_t pythiaType);
  Short_t GetDPMjetEventProcessType(Int_t dpmJetType);
  ULong64_t GetEventIdAsLong(AliVHeader* header) const;

  TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
  Int_t FindPrimaryMotherLabel(AliStack* stack, Int_t label);

  AliAODMCParticle* FindPrimaryMotherAOD(AliAODMCParticle* startParticle);

  TParticle* FindPrimaryMotherV0(AliStack* stack, Int_t label);
  Int_t FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps);
  Bool_t PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t mag,
                TF1* phiCutLow, TF1* phiCutHigh);

  AliAODMCParticle* FindPrimaryMotherAODV0(AliAODMCParticle* startParticle,
                                           Int_t& nSteps);

  void GetSPDMultiplicity();
  void GetCalibratedV0Amplitude();
  void MultiplicityDistributions();

  AliESDEvent* fESD;                //! ESD object
  AliAODEvent* fAOD;                //! AOD object
  AliMCEvent* fMC;                  //! MC object
  AliStack* fMCStack;               //! MC ESD stack
  TClonesArray* fMCArray;           //! MC array for AOD
  AliAnalysisFilter* fTrackFilter;  // Track filter, same as PbPb 5.02 TeV
  TString fCentEst;                 // V0A , V0M,
  TString fAnalysisType;            //  "ESD" or "AOD"
  Bool_t fAnalysisMC;               //  Real(kFALSE) or MC(kTRUE) flag
  Bool_t fAnalysisPbPb;  //  true you want to analyze PbPb data, false for pp
  UInt_t ftrigBit;
  TRandom* fRandom;   //! random number generator
  Bool_t fPileUpRej;  // kTRUE is pile-up is rejected

  //
  // Cuts and options
  //

  Double_t fVtxCut;    // Vtx cut on z position in cm
  Double_t fEtaCut;    // Eta cut used to select particles
  Double_t fEtaMin;    // Min eta cut for TPC estimator
  Double_t fEtaMax;    // Max eta cut for TPC estimator
  Double_t fPtMin;     // Min pT cut to select particles
  Double_t fEtaGappT;  // pT cut to for eta gap estimator
  Float_t fMinCent;    // minimum centrality
  Float_t fMaxCent;    // maximum centrality
  Bool_t fStoreMcIn;   // Store MC input tracks
  Float_t fEtaGapNchMin;
  Float_t fEtaGapNchMax;
  //
  // Help variables
  //
  Short_t fMcProcessType;     // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
  Short_t fTriggeredEventMB;  // 1 = triggered, 0 = not trigged (MC only)
  Short_t fVtxStatus;         // -1 = no vtx, 0 = outside cut, 1 = inside cut
  Float_t fZvtx;              // z vertex
  Float_t fZvtxMC;            // z vertex MC (truth)
  Int_t fRun;                 // run no
  ULong64_t fEventId;         // unique event id
  Float_t fv0mamplitude;
  Float_t fv0mpercentile;

  int fTracklets14;      // Tracklets estimator |eta|<1.4
  int fTracklets10;      // Tracklets estimator |eta|<1.0
  int fTrackletsEtaGap;  // Tracklets estimator 0.7<|eta|<1.4

  //
  // Output objects
  //
  TList* fListOfObjects;  //! Output list of objects
  TH1I* fEvents;          //! No of accepted events
  TH1I* fVtx;             //! Event vertex info
  TH1F* fVtxMC;           //! Event vertex info for ALL MC events
  TH1F* fVtxBeforeCuts;   //! Vertex z dist before cuts
  TH1F* fVtxAfterCuts;    //! Vertex z dist after cuts
  TH1F* fn1;
  TH1F* fcent;

  TH2D* hPhi;

  TH2F* hPhiEtaSPD;
  TH2F* hPhiEtaGapSPD;
  TH2F* hVtxZvsTracklets;
  TProfile* pV0MAmpChannel;
  TH2D* hPtEtaNegvsNchEtaPos;  // From here
  TProfile* pPtEtaNegvsNchEtaPos;
  TH2D* hPtEtaPosvsNchEtaNeg;
  TProfile* pPtEtaPosvsNchEtaNeg;
  TProfile* pPtvsNch;
  TProfile* pPtvsV0MAmp;
  TH2D* hPtvsV0MAmp;
  TH2D* hPtvsTracklets10;
  TProfile* pPtvsTracklets10;
  TH2D* hPtvsTracklets14;
  TProfile* pPtvsTracklets14;
  TH2D* hPtvsTrackletsEtaGap;
  TProfile* pPtvsTrackletsEtaGap;
  TH1F* hV0Mmult;
  TH1F* hV0MAmplitude;
  TH2F* hNchvsV0M;
  TH2F* hNchvsV0MAmp;
  TH2F* hV0MvsV0MAmp;
  TH2D* hNchEtaPosvsNchEtaNeg;
  TH2F* hTrackletsvsV0MAmp10;
  TH2F* hTrackletsvsV0MAmp14;
  TH1F* hTrackletsEtaGap;
  TH2F* hTrackletsvsV0MAmpEtaGap;

  TH1D* hMcIn[7][9];
  TH1D* hMcOut[7][9];

  AliAnalysisTaskDataSpeedOfSound276TeV(
      const AliAnalysisTaskDataSpeedOfSound276TeV&);  // not implemented
  AliAnalysisTaskDataSpeedOfSound276TeV& operator=(
      const AliAnalysisTaskDataSpeedOfSound276TeV&);  // not implemented

  ClassDef(AliAnalysisTaskDataSpeedOfSound276TeV,
           1);  // Analysis task for high pt analysis
};

#endif
