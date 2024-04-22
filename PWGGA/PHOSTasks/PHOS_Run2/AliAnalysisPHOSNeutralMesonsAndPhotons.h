#ifndef AliAnalysisPHOSNeutralMesonsAndPhotons_cxx
#define AliAnalysisPHOSNeutralMesonsAndPhotons_cxx

#include "AliCaloPhoton.h"
#include "AliAODCaloCluster.h"
#include "AliPHOSGeometry.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliCaloTriggerMimicHelper.h"
#include "AliAnalysisUtils.h"
// #include "AliMCAnalysisUtils.h"

// Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// See cxx source for full Copyright notice

// Analysis task for neutral meson and inclusive/direct photon measurements (AODs only)
// + supproting QA studies
// Author: Vladislav Kuskov
// since march 2024

#include "AliAnalysisTaskSE.h"

class AliAnalysisPHOSNeutralMesonsAndPhotons : public AliAnalysisTaskSE
{
 public:
  AliAnalysisPHOSNeutralMesonsAndPhotons(const char* name = "AliAnalysisPHOSNeutralMesonsAndPhotons");
  virtual ~AliAnalysisPHOSNeutralMesonsAndPhotons();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  void SetMBTrigger(Bool_t isMB) { fIsMBTriggerAnalysis = isMB; };
  void SetPHOSTrigger(Bool_t isPHS, AliCaloTriggerMimicHelper::phosTriggerType PHSTrigType)
  {
    fIsPHOSTriggerAnalysis = isPHS;
    fPHOSTrigType = PHSTrigType;
  };
  void SetMC(Bool_t isMC) { fIsMC = isMC; };
  void SetCoreEnergyFlag(Bool_t useCoreE) { fUseCoreEnergy = useCoreE; };
  void SetClusterCuts(Bool_t useCoreE, Double_t Emin, Double_t DispCut, Double_t CPVCut, Double_t TOFCut, Bool_t MinBiasLowECut)
  {
    fUseCoreEnergy = useCoreE;
    fEminCut = Emin;
    fDispSigma = DispCut;
    fCPVSigma = CPVCut;
    fTOFCut = TOFCut;
    fMinBiasLowEClustCut = MinBiasLowECut;
  };
  void SetQAStudy(Bool_t DoClustQA, Bool_t DoCellQA)
  {
    fDoClustQA = DoClustQA;
    fDoCellQA = DoCellQA;
  };
  void SetCentralityStudy(Bool_t doCentralityStudy, TString CentEstim, Float_t CentMin, Float_t CentMax)
  {
    fDoCentralityStudy = doCentralityStudy;
    fCentralityEstimator = CentEstim;
    fCentralityMin = CentMin;
    fCentralityMax = CentMax;
  }

  void SetTOFEfficiency(TF1* tofeff) { fUserTOFEff = tofeff; };
  void SetUserNonlinearityFunction(TF1* nonlin) { fUserNonlinFunc = nonlin; };
  void SetPHOSModEnergyCorr(Double_t PHOSWeights[4]){
    for (Int_t imod = 0; imod < 4; imod++)
      fPHOSModEnCorr[imod] = PHOSWeights[imod];
  };

  void SetNEventsForMixing(Int_t NMix) { fNMixEvents = NMix; };

 private:
  AliAnalysisPHOSNeutralMesonsAndPhotons(const AliAnalysisPHOSNeutralMesonsAndPhotons&);
  AliAnalysisPHOSNeutralMesonsAndPhotons& operator=(const AliAnalysisPHOSNeutralMesonsAndPhotons&);

  void FillHistogram(const char* key, Double_t x) const;                         // Fill 1D histogram witn name key
  void FillHistogram(const char* key, Double_t x, Double_t y) const;             // Fill 2D histogram witn name key
  void FillHistogram(const char* key, Double_t x, Double_t y, Double_t z) const; // Fill 3D histogram witn name key

  void SelectCentrality();
  void InitPHOSGeometry();

  void ProcessCaloPhotons();
  void FillMgg();
  void EstimatePIDCutEfficiency();
  void EstimateTOFCutEfficiency();

  void ProcessMCParticles();
  AliAODMCParticle* FindMCPrimaryMother(AliAODMCParticle* particle);
  AliAODMCParticle* FindMCPrimaryMother(Int_t primlb);
  Int_t FindMCPrimaryMotherLabel(Int_t primlb);
  Bool_t InPHOSAcc(AliAODMCParticle* particle);
  Bool_t IsPrimaryMCParticle(AliAODMCParticle* particle);
  Bool_t IsPrimaryMCParticle3D(AliAODMCParticle* particle);

  TF1* GetNonlinFunction() { return fUserNonlinFunc; };

 private:
  static const Int_t kMods = 4;     // number of PHOS modules
  static const Int_t kPIDCuts = 4;  // PHOS PID cuts
  static const Int_t kCentBins = 7; // centrality
  static const Int_t kVtxBins = 20; // z-vertex
  static const Int_t kPRBins = 6;   // Reaction plane

  THashList* fOutputContainer; //! final histogram container

  AliMCEvent* fMCEvent;
  AliAODEvent* fEvent; //! Current event

  Int_t fRunNumber;

  TClonesArray* fMCStack;

  TClonesArray* fCaloPhotonsPHOS; //! List of selected photons in current event
  TClonesArray* fCaloPhotonsMix;  //! List of selected photons in mixed events

  TList* fCaloPhotonsPHOSList;
  TList* fPHOSEvents[kVtxBins][kCentBins][kPRBins]; //! Previous events for mixing

  AliAnalysisUtils* fUtils;
  AliPHOSGeometry* fPHOSGeo;
  AliAODCaloCells* fAODCells;
  AliCaloTriggerMimicHelper* fCaloTriggerMimicHelper;
  AliCaloTriggerMimicHelper::phosTriggerType fPHOSTrigType;

  Bool_t fIsMC;
  Bool_t fIsPHOSTriggerAnalysis;
  Bool_t fIsMBTriggerAnalysis;
  Bool_t fDoClustQA;
  Bool_t fDoCellQA;

  Int_t fNMixEvents;

  Double_t fVertex[3];
  Int_t fZvtx;

  Double_t fRP; // readction plane

  Bool_t fDoCentralityStudy;
  TString fCentralityEstimator;
  AliMultSelection* fMultSelection;
  Float_t fCentrality;
  Int_t fNCenBin;
  Int_t fCentBin;
  TArrayI fCenBinEdges;
  Double_t fCentralityMin;
  Double_t fCentralityMax;

  Bool_t fDoNonlinCorr;
  TF1* fUserNonlinFunc;
  Bool_t fDoTOFEffCorr;
  TF1* fUserTOFEff;
  Double_t fPHOSModEnCorr[kMods];

  Bool_t fUseCoreEnergy;
  Double_t fEminCut;
  Double_t fDispSigma;
  Double_t fCPVSigma;
  Double_t fTOFCut;
  Bool_t fMinBiasLowEClustCut;

  TH1F* fHistInfo;
  TH1F* fHistSelectEvent;
  TH1F* fHistVertexZ;

  TH1F* fHistCentMain;
  TH2F* fHistCentV0MvsCL0;
  TH2F* fHistCentV0MvsCL1;
  TH2F* fHistCentV0MvsV0C;
  TH2F* fHistCentCL0vsCL1;
  TH2F* fHistCentZNAvsZNC;

  TH2I* fPHOSBadMap[6];

  TH2F* fHistTOFClust;
  TH2F* fHistTOFClustMod[kMods];
  TH2F* fHistCaloPhotonTOFvsE;
  TH2F* fHistCaloPhotonTOFvsEMod[kMods];
  TH1F* fHistCaloPhotonTOFCut;

  TH2F* fHistClustTOFvsDDL;
  TH2F* fHistClustTOFvsDDLEnCut;

  TH2F* fHistClustMultVsCentrality;
  TH2F* fHistClustMultVsCentralityMod[kMods];

  TH1F* fHistCaloPhotonPt[kPIDCuts];

  TH2F* fHistM02vsPt;
  TH2F* fHistM20vsPt;
  TH2F* fHistClustEvsXZMod[kMods];

  TH1F* fHistClustFullE;
  TH1F* fHistClustCoreE;
  TH1F* fHistCaloPhotonPtMod[kMods];

  TH2F* fHistNonlinTest;

  TH1F* fHistClustFullEMod[kMods];
  TH1F* fHistClustCoreEMod[kMods];

  TH2F* fHistMgg[kPIDCuts];
  TH2F* fHistMixMgg[kPIDCuts];

  TH2F* fHistMggTOFCutEffBase;
  TH2F* fHistMixMggTOFCutEffBase;
  TH2F* fHistMggTOFCutEffProbe;
  TH2F* fHistMixMggTOFCutEffProbe;

  TH2F* fHistMggCutEff[kPIDCuts];
  TH2F* fHistMixMggCutEff[kPIDCuts];

  TH2F* fHistMggCutEffMod[kMods][kPIDCuts];
  TH2F* fHistMixMggCutEffMod[kMods][kPIDCuts];

  TH2F* fHistMggMod[kMods];
  TH2F* fHistMixMggMod[kMods];
  TH2F* fHistMggPhIDCutMod[kMods];
  TH2F* fHistMixMggPhIDCutMod[kMods];

  TH2F* fHistTruePi0MggVsRecPt;
  TH2F* fHistTruePi0MggVsTruePt;
  TH2F* fHistTrueEtaMggVsRecPt;
  TH2F* fHistTrueEtaMggVsTruePt;

  TH2F* fHistMCPartIDvsPt;
  TH2F* fHistMCCaloPartIDvsPt[kPIDCuts];

  TH1F* fHistPrimPi0Pt;
  TH1F* fHistPrimPi0InAccPt;
  TH1F* fHistPrimPi0BothPhInAccPt;
  TH1F* fHistPrimPi0OnePhInAccPt;
  TH2F* fHistPrimPi0PtvsPhInAccPt;

  TH1F* fHistPrimEtaPt;
  TH1F* fHistPrimEtaInAccPt;
  TH1F* fHistPrimEtaBothPhInAccPt;
  TH1F* fHistPrimEtaOnePhInAccPt;
  TH2F* fHistPrimEtaPtvsPhInAccPt;

  TH1F* fHistPrimGammaPt;
  TH1F* fHistPrimGammaInAccPt;

  // Histograms for analysis

  ClassDef(AliAnalysisPHOSNeutralMesonsAndPhotons, 1);
};
#endif
