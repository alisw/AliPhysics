#ifndef ALIANALYSISTASKSEDVSEVENTSHAPES_H
#define ALIANALYSISTASKSEDVSEVENTSHAPES_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// Class AliAnalysisTaskSEDvsEventShapes
// AliAnalysisTaskSE for the D meson vs. Event shape analysis in different mutiplicity window
// Authors: Renu Bala, Manoj Bhanudas Jadhav
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TArrayD.h>
#include <TFile.h>
#include <TRandom.h>
#include <TProfile.h>
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliNormalizationCounter.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliVertexingHFUtils.h"
#include "AliVEvent.h"

class AliAnalysisTaskSEDvsEventShapes : public AliAnalysisTaskSE
{
public:
    
    AliAnalysisTaskSEDvsEventShapes();
    AliAnalysisTaskSEDvsEventShapes(const char *name, Int_t pdgMeson, AliRDHFCuts* cuts, Bool_t switchPPb);
    virtual ~AliAnalysisTaskSEDvsEventShapes();
    
    void SetMassLimits(Double_t lowlimit, Double_t uplimit);
    void SetMassLimits(Int_t pdg, Double_t range);
    Double_t GetUpperMassLimit() const {return fUpmasslimit;}
    Double_t GetLowerMassLimit() const {return fLowmasslimit;}
    void SetNMassBins(Int_t nbins){fNMassBins=nbins;}
    Int_t GetNMassBins() const{return fNMassBins;}
    Bool_t GetSubtractTrackletsFromDaughters() const {return fSubtractTrackletsFromDau;}
    Bool_t GetRecomputeSpherocityWithoutDau() const {return fRecomputeSpherocity;}
    Bool_t GetRemoveD0fromDstar() const {return fRemoveD0fromDstar;}
    
    void SetImpactParameterBinning(Int_t nbins, Double_t dmin, Double_t dmax){
        fNImpParBins=nbins;
        fLowerImpPar=dmin;
        fHigherImpPar=dmax;
    }
    void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
    void SetMCOption(Int_t option=0){ fMCOption = option; }
    void SetIsPPbData(Bool_t flag=kTRUE){
        fisPPbData=flag;
    }
    void SetUseBit(Bool_t use=kTRUE){fUseBit=use;}
    void SetDoImpactParameterHistos(Bool_t doImp=kTRUE){fDoImpPar=doImp;}
    
    void SetMultiplVsZProfileLHC10b(TProfile* hprof){
        if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
        fMultEstimatorAvg[0]=new TProfile(*hprof);
    }
    void SetMultiplVsZProfileLHC10c(TProfile* hprof){
        if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
        fMultEstimatorAvg[1]=new TProfile(*hprof);
    }
    void SetMultiplVsZProfileLHC10d(TProfile* hprof){
        if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
        fMultEstimatorAvg[2]=new TProfile(*hprof);
    }
    void SetMultiplVsZProfileLHC10e(TProfile* hprof){
        if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
        fMultEstimatorAvg[3]=new TProfile(*hprof);
    }
    
    void SetMultiplVsZProfileLHC13b(TProfile* hprof){
        if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
        fMultEstimatorAvg[0]=new TProfile(*hprof);
    }
    void SetMultiplVsZProfileLHC13c(TProfile* hprof){
        if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
        fMultEstimatorAvg[1]=new TProfile(*hprof);
    }
    
    void SetReferenceMultiplcity(Double_t rmu){fRefMult=rmu;}
    
    // Nch Ntrk weights on MC
    void UseMCNchWeight(Int_t flag) { fUseNchWeight = flag; }
    void SetHistoNchWeight(TH1F *h){
        if(fHistoMCNch) delete fHistoMCNch;
        fHistoMCNch = new TH1F(*h);
    }
    void SetMeasuredNchHisto(TH1F* h){
        if(fHistoMeasNch) delete fHistoMeasNch;
        fHistoMeasNch = new TH1F(*h);
    }
    
    // pT weights on MC
    void UsePtWeight(Bool_t flag) { fUsePtWeight = flag; }
    Double_t GetPtWeight(Float_t pt);
    Double_t dNdptFit(Float_t pt, Double_t* par);
    
    void SetSubtractTrackletsFromDaughters(Bool_t opt){fSubtractTrackletsFromDau=opt;}
    
    // Flag to use the zvtx correction from ( 0= none, 1= usual d2h, 2=AliESDUtils for VZERO multiplicity)
    void SetUseVZEROParameterizedVertexCorr(Int_t flag) { fDoVZER0ParamVertexCorr=flag; }
    // Flag to fill THnSparse with MultUncorr and NoPid cases ( 0 = only Mult, 1 = Mult and multUncorr, 2 = NoPid and 3 is All)
    void SetFillSoSparseForMultUncorrNoPid(Int_t flag) { fFillSoSparseChecks=flag; }
    void SetEventShapeParameters(Double_t ptMin, Double_t ptMax, Double_t etaMin, Double_t etaMax, Int_t minMult, Double_t phiStepSizeDeg, Int_t filtbit1, Int_t filtbit2) { fptMin=ptMin; fptMax=ptMax; fetaMin=etaMin; fetaMax=etaMax; fminMult=minMult; fphiStepSizeDeg=phiStepSizeDeg; ffiltbit1=filtbit1; ffiltbit2=filtbit2;}
    
    void SetCalculationsForSphericity(Bool_t CalSpheri){fCalculateSphericity=CalSpheri;}
    void SetRecomputeSpherocityWithoutDau(Bool_t RecomputeSphero){fRecomputeSpherocity=RecomputeSphero;}
    void SetRemoveD0fromDstar(Bool_t RemoveD0fromDstar){fRemoveD0fromDstar=RemoveD0fromDstar;}
    
    Int_t GetUseVZEROParameterizedVertexCorr() { return fDoVZER0ParamVertexCorr; }
    
    enum { kNtrk10=0, kNtrk10to16=1, kVZERO=2, kNtrk03=3, kNtrk05=4, kVZEROA=5, kVZEROEq=6, kVZEROAEq=7 };
    void SetMultiplicityEstimator(Int_t value){ fMultiplicityEstimator=value; }
    Int_t GetMultiplicityEstimator(){ return fMultiplicityEstimator; }
    enum { kEta10=0, kEta10to16=1, kEtaVZERO=2, kEta03=3, kEta05=5, kEtaVZEROA=5 };
    void SetMCPrimariesEstimator(Int_t value){ fMCPrimariesEstimator=value; }
    Int_t GetMCPrimariesEstimator(){ return fMCPrimariesEstimator; }
    
    void SetUseQuarkLevel(Bool_t opt){fUseQuarkTag=opt;}
    void SetEtaAccCut(Double_t etacut){fEtaAccCut=etacut;}
    void SetPtAccCut(Double_t ptcut){fPtAccCut=ptcut;}
    Bool_t CheckGenAcc(TClonesArray* arrayMC, Int_t nProng, Int_t *labDau);
    void SetKeepTrackControlHisto(Bool_t flag){fFillTrackHisto=flag;}
    
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    
private:
    
    AliAnalysisTaskSEDvsEventShapes(const AliAnalysisTaskSEDvsEventShapes &source);
    AliAnalysisTaskSEDvsEventShapes& operator=(const AliAnalysisTaskSEDvsEventShapes& source);
    
    TProfile* GetEstimatorHistogram(const AliVEvent *event);
    void CreateImpactParameterHistos();
    void CreateMeasuredNchHisto();
    Bool_t FillTrackControlHisto(AliAODEvent* aod, Int_t nSelTrkCorr, Double_t spherocity, Double_t genspherocity, Int_t nSelectedEvwithCand);
    void FillMCMassHistos(TClonesArray *arrayMC, Int_t labD, Double_t countMult, Double_t spherocity, Double_t sphericity, Double_t recSpherocity, Double_t nchWeight);
    void FillMCGenAccHistos(AliAODEvent* aod, TClonesArray *arrayMC, AliAODMCHeader *mcHeader, Double_t countMult, Double_t spherocity, Double_t sphericity, Bool_t isEvSel, Double_t nchWeight);
    
    TList  *fOutput; //! list send on output slot 1
    TList  *fListCuts; // list of cuts
    TList  *fOutputCounters; //! list send on output slot 3
    TList  *fListProfiles; // list of profile histos for z-vtx correction
    TList  *fOutputEffCorr; //! list send on output slot 5
    
    TH1F *fHistNEvents;     //! hist. for No. of events
    
    TH2F* fHistNtrVsZvtx; //!  hist of ntracklets vs Zvertex
    TH2F* fHistNtrCorrVsZvtx; //!  hist of ntracklets vs Zvertex
    TH2F* fHistNtrVsnTrackEvWithCand; //!<!  control hist of ntracklets vs nTracks passing track selection for spherocity calculation for event with atleast one D meson
    TH2F* fHistNtrVsSo; //!  hist of ntracklets vs So
    TH2F* fHistNtrCorrVsSo; //!  hist of ntracklets vs So
    TH2F* fHistNtrVsSpheri; //!  hist of ntracklets vs Spheri
    TH2F* fHistNtrCorrVsSpheri; //!  hist of ntracklets vs Spheri
    
    TH2F* fHistNtrVsNchMC; //!<!  hist of ntracklets vs Nch (Generated)
    TH2F* fHistNtrCorrVsNchMC; //!<!  hist of ntracklets vs Nch (Generated)
    TH2F* fHistNtrVsNchMCPrimary; //!<!  hist of ntracklets vs Nch (Primary)
    TH2F* fHistNtrCorrVsNchMCPrimary; //!<!  hist of ntracklets vs Nch (Primary)
    TH2F* fHistNtrVsNchMCPhysicalPrimary; //!<!  hist of ntracklets vs Nch (Physical Primary)
    TH2F* fHistNtrCorrVsNchMCPhysicalPrimary; //!<!  hist of ntracklets vs Nch (Physical Primary)
    TH1F* fHistGenPrimaryParticlesInelGt0; //!<!hist. of geenrated multiplcity
    TH3F* fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary; //!<! hist of Nch (generated) vs Nch (Primary) vs Nch (Physical Primary)
    
    TH1F* fHistNtrCorrPSSel; //! hist. of ntracklets for physics selection only selected events
    TH1F* fHistNtrCorrEvSel; //! hist. of ntracklets for selected events
    TH1F* fHistNtrCorrEvWithCand; //! hist. of ntracklets for evnts with a candidate
    TH1F* fHistNtrCorrEvWithD;//! hist. of ntracklets for evnts with a candidate in D mass peak
    
    TH3F *fHistnTrackvsEtavsPhi;  //!<! hist. of number of tracks passing track selection for spherocity calculation vs eta vs. phi
    TH3F *fHistnTrackvsEtavsPhiEvWithCand;  //!<! hist. of number of tracks passing track selection for spherocity calculation vs eta vs. phi for event with atleast one D meson
    
    TH3F *fHistTrueSovsMeasSo;  //!<! hist. of number of tracks passing track selection for spherocity calculation vs eta vs. phi
    TH3F *fHistTrueSovsMeasSoEvWithCand;  //!<! hist. of number of tracks passing track selection for spherocity calculation vs eta vs. phi for event with atleast one D meson
    
    TH3F *fHistSpheroAxisDeltaPhi;  //!<! hist. of Invariant mass, pt vs. deltaPhi of spherocity axis w.r.t. D-meson direction
    TH3F *fHistSpheroAxisDeltaGenPhi;  //!<! hist. of Invariant mass, pt vs. deltaPhi of generated spherocity axis w.r.t. D-meson direction

    THnSparseD *fSparseEvtShape;//! THnSparse histograms for Spherocity
    THnSparseD *fSparseEvtShapewithNoPid;//! THnSparse histograms for D0 vs. Spherocity
    THnSparseD *fSparseEvtShapePrompt;//! THnSparse histograms for Prompt D0 vs. Spherocity
    THnSparseD *fSparseEvtShapeFeeddown;//! THnSparse histograms for feeddown D0 vs. Spherocity
    THnSparseD *fSparseEvtShapeRecSphero;//! THnSparse histograms for Both Prompt and feeddown D0 vs. Spherocity
    THnSparseD *fMCAccGenPrompt; //! histo for StepMCGenAcc for D meson prompt
    THnSparseD *fMCAccGenFeeddown; //! histo for StepMCGenAcc for D meson feeddown
    THnSparseD *fMCRecoPrompt; //! histo for StepMCReco for D meson feeddown
    THnSparseD *fMCRecoFeeddown; //! histo for StepMCReco for D meson feeddown
    THnSparseD *fMCRecoBothPromptFD; //! histo for StepMCReco for D meson Both Prompt Feeddown
    THnSparseD *fMCAccGenPromptSpheri; //! histo for StepMCGenAcc for D meson prompt for Sphericity
    THnSparseD *fMCAccGenFeeddownSpheri; //! histo for StepMCGenAcc for D meson feeddown for Sphericity
    THnSparseD *fMCRecoPromptSpheri; //! histo for StepMCReco for D meson feeddown for Sphericity
    THnSparseD *fMCRecoFeeddownSpheri; //! histo for StepMCReco for D meson feeddown for Sphericity
    THnSparseD *fMCRecoBothPromptFDSpheri; //! histo for StepMCReco for D meson Both Prompt Feeddown for Sphericity
    
    THnSparseD *fMCAccGenPromptEvSel; //! histo for StepMCGenAcc for D meson prompt with Vertex selection (IsEvSel = kTRUE)
    THnSparseD *fMCAccGenFeeddownEvSel; //! histo for StepMCGenAcc for D meson feeddown with Vertex selection (IsEvSel = kTRUE)
    
    THnSparseF *fHistMassPtImpPar[5];//! histograms for impact paramter studies
    
    Double_t fUpmasslimit;  //upper inv mass limit for histos
    Double_t fLowmasslimit; //lower inv mass limit for histos
    Int_t   fNMassBins;    // nbins for invariant mass histos
    
    AliRDHFCuts *fRDCutsAnalysis; // Cuts for Analysis
    AliNormalizationCounter *fCounterC;           //! Counter for normalization, corrected multiplicity
    AliNormalizationCounter *fCounterU;           //! Counter for normalization, uncorrected multiplicity
    AliNormalizationCounter *fCounterCandidates;  //! Counter for normalization, corrected multiplicity for candidates
    
    Bool_t fDoImpPar;  //swicth for D impact parameter THnSparse
    Int_t  fNImpParBins;   // nunber of bins in impact parameter histos
    Double_t fLowerImpPar;  // lower limit in impact parameter (um)
    Double_t fHigherImpPar; // higher limit in impact parameter (um)
    
    Bool_t fReadMC;    //flag for access to MC
    Int_t  fMCOption;  // 0=keep all cand, 1=keep only signal, 2= keep only back
    Bool_t fisPPbData; // flag to run on pPb data (differen histogram bining)
    Bool_t fUseBit;    // flag to use bitmask
    Bool_t fSubtractTrackletsFromDau; // flag for subtracting D meson daughter contribution to N of tracklets
    Bool_t fCalculateSphericity; // flag for computing Sphericity
    Bool_t fRecomputeSpherocity; // flag for subtracting D meson daughter contribution to Spherocity calculation
    Bool_t fRemoveD0fromDstar; // flag for removal of D0 from D* meson
    Int_t fUseNchWeight; // weight on the MC on the generated multiplicity (0->no weights, 1->Nch weights, 2->Ntrk weights)
    TH1F* fHistoMCNch;    // weight histogram for the MC on the generated multiplicity
    TH1F* fHistoMeasNch;  // weight histogram on the true measured multiplicity
    Bool_t fUsePtWeight; // weight on the MC on the generated pT
    Double_t fWeight; // Total weight on the MC: nchWeight*ptWeight
    
    TProfile* fMultEstimatorAvg[4]; //TProfile with mult vs. Z per period
    Double_t fRefMult;   // refrence multiplcity (period b)
    Int_t fPdgMeson;   // pdg code of analyzed meson
    
    Int_t fMultiplicityEstimator; // Definition of the multiplicity estimator: kNtrk10=0, kNtrk10to16=1, kVZERO=2
    Int_t fMCPrimariesEstimator;  // Definition of the primaries estimator eta range: |eta|<1.0=0, -1.6<|eta|<1.0=1, VZEROrange=2
    
    Int_t fDoVZER0ParamVertexCorr; // Flag to use the zvtx correction from (0=none, 1=usual d2h, 2=AliESDUtils for VZERO multiplicity)
    
    Int_t fFillSoSparseChecks; // Flag to fill THnSparse with MultUncorr and NoPid cases ( 0 = only Mult, 1 = Mult and multUncorr, 2 = NoPid and 3 is All)
    
    Bool_t fUseQuarkTag; /// flag for quark/hadron level identification of prompt and feeddown
    Bool_t fFillTrackHisto; /// flag for filling track control histograms
    Double_t fEtaAccCut; /// eta limits for acceptance step
    Double_t fPtAccCut; /// pt limits for acceptance step
    
    Double_t fetaMin;
    Double_t fetaMax;
    Double_t fptMin;
    Double_t fptMax;
    Int_t fminMult;
    Int_t ffiltbit1;
    Int_t ffiltbit2;
    Double_t fphiStepSizeDeg;
    
    ClassDef(AliAnalysisTaskSEDvsEventShapes,14); // D vs. mult task
};

#endif
