#ifndef AliAnalysisTaskSED0BDT_H
#define AliAnalysisTaskSED0BDT_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///*************************************************************************
/// \class Class AliAnalysisTaskSED0BDT
/// \brief AliAnalysisTaskSE for D0 candidates invariant mass histogram
/// and comparison to MC truth (kinematics stored in the AOD) and cut variables
/// distributions
/// \author Authors: A.Dainese, andrea.dainese@ln.infn.it
/// \author and C.Bianchin, chiara.bianchin@pd.infn.it
///*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <THnSparse.h>
#include <TString.h>
#include <TRandom3.h>
#include <TProfile.h>
#include "AliESDUtils.h"
#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsD0toKpiBDT.h"
//~ #include "AliRDHFBDT.h"
#include "AliNormalizationCounter.h"
#include "AliEventCuts.h"

class AliAODEvent;

class AliAnalysisTaskSED0BDT : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSED0BDT();
  AliAnalysisTaskSED0BDT(const char *name,AliRDHFCutsD0toKpiBDT* cuts);
  virtual ~AliAnalysisTaskSED0BDT();


  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void CreateMCAcceptanceHistos();
  Bool_t CheckAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau);
  void FillMCAcceptanceHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader);

  void ProcessBDT(AliAODEvent *aod, AliAODRecoDecayHF2Prong *part,TClonesArray *arrMC, Float_t multi);
  void ProcessBDT2(AliAODEvent *aod, AliAODRecoDecayHF2Prong *part,TClonesArray *arrMC, Float_t multi);
  void SetArray(Int_t type=AliAnalysisTaskSED0BDT::kD0){fArray=type;}
  enum{kD0,kLS};

  void SetReadMC(Bool_t readMC=kFALSE){fReadMC=readMC;}
  void SetDoMCAcceptanceHistos(Bool_t doMCAcc=kTRUE){fStepMCAcc=doMCAcc;}
  void SetCutOnDistr(Bool_t cutondistr=kFALSE){fCutOnDistr=cutondistr;}
  void SetUsePid4Distr(Bool_t usepid=kTRUE){fUsePid4Distr=usepid;}
  void SetFillOnlyD0D0bar(Int_t flagfill){fFillOnlyD0D0bar=flagfill;}
  void SetFillVarHists(Bool_t flag) {fFillVarHists=flag;}
  void SetFillPtHistos(Bool_t flag) {fFillPtHist=flag;}
  void SetFillYHistos(Bool_t flag) {fFillYHist=flag;}
  void SetFillImpactParameterHistos(Bool_t flag) {fFillImpParHist=flag;}
  void SetFillSparses(Bool_t flag) {fFillSparses=flag;}
  void SetUseRejectionMethod(Bool_t flag=kFALSE, Float_t factor=0.01) {fUseRejectionMethod=flag; fRejectionFactor=factor;}
  void SetSystem(Int_t sys){fSys=sys; if(fSys==1) SetFillVarHists(kFALSE);}
  void SetRejectSDDClusters(Bool_t flag) { fIsRejectSDDClusters=flag; }
  void SetUseSelectionBit(Bool_t flag) { fUseSelectionBit=flag; }
  void SetWriteVariableTree(Bool_t flag) { fWriteVariableTree=flag; }
  void SetWriteProtosgnVar(Bool_t flag) { fWriteProtosgnVar=flag; }
  void SetSelectTrueD0(Bool_t flag) { fSelectTrueD0 = flag; }
  void SetUseMassWindow(Bool_t flag) { fUsedMassWindow=flag; }
  void SetDrawDetSignal(Bool_t flag) { fDrawDetSignal=flag; }
  void SetPIDCheck(Bool_t flag) { fPIDCheck=flag; }
  void SetUseQuarkLevelTag(Bool_t opt){fUseQuarkTagInKine=opt;}
  void SetAODMismatchProtection(Int_t opt=0) {fAODProtection=opt;}
  void SetPileupRejectionVZEROTPCout(Bool_t flag) {fEnablePileupRejVZEROTPCout=flag;}
  void SetPileupRejectionVZEROTPCcls(Bool_t flag, Bool_t rejpileup) {fEnablePileupRejVZEROTPCcls=flag; fRejectOutOfBunchPileUp=rejpileup;}
  void SetFillSubSampleHist(Bool_t flag) {fFillSubSampleHist=flag;}
    
    void SetSubtractTrackletsFromDaughters(Bool_t opt){fSubtractTrackletsFromDau=opt;}
    void SetMultiana(Bool_t flag) { fmultiana=flag; }
    void SetReferenceMultiplcity(Double_t rmu){fRefMult=rmu;}
    void SetMultiplicityEstimator(Int_t value){ fMultiplicityEstimator=value; }
    enum { kNtrk10=0, kNtrk10to16=1, kVZERO=2, kNtrk03=3, kNtrk05=4, kVZEROA=5, kVZEROEq=6, kVZEROAEq=7 };
    Int_t GetMultiplicityEstimator(){ return fMultiplicityEstimator; }
    enum { kEta10=0, kEta10to16=1, kEtaVZERO=2, kEta03=3, kEta05=5, kEtaVZEROA=5 };
    void SetMCPrimariesEstimator(Int_t value){ fMCPrimariesEstimator=value; }
    Int_t GetMCPrimariesEstimator(){ return fMCPrimariesEstimator; }
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
      fYearNumber = 13;
    }
    void SetMultiplVsZProfileLHC13c(TProfile* hprof){
      if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
      fMultEstimatorAvg[1]=new TProfile(*hprof);
      fYearNumber = 13;
    }
      
    void SetMultiplVsZProfileLHC16qt1stBunch(TProfile* hprof){
      if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
      fMultEstimatorAvg[0]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16qt2ndBunch(TProfile* hprof){
      if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
      fMultEstimatorAvg[1]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16qt3rdBunch(TProfile* hprof){
      if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
      fMultEstimatorAvg[2]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16qt4thBunch(TProfile* hprof){
      if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
      fMultEstimatorAvg[3]=new TProfile(*hprof);
      fYearNumber = 16;
    }

    void SetMultiplVsZProfileLHC16d(TProfile* hprof){
      if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
      fMultEstimatorAvg[0]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16e(TProfile* hprof){
      if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
      fMultEstimatorAvg[1]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16g(TProfile* hprof){
      if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
      fMultEstimatorAvg[2]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16h1(TProfile* hprof){
      if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
      fMultEstimatorAvg[3]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16h2(TProfile* hprof){
      if(fMultEstimatorAvg[4]) delete fMultEstimatorAvg[4];
      fMultEstimatorAvg[4]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16j(TProfile* hprof){
      if(fMultEstimatorAvg[5]) delete fMultEstimatorAvg[5];
      fMultEstimatorAvg[5]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16k(TProfile* hprof){
      if(fMultEstimatorAvg[6]) delete fMultEstimatorAvg[6];
      fMultEstimatorAvg[6]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16l(TProfile* hprof){
      if(fMultEstimatorAvg[7]) delete fMultEstimatorAvg[7];
      fMultEstimatorAvg[7]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16o(TProfile* hprof){
      if(fMultEstimatorAvg[8]) delete fMultEstimatorAvg[8];
      fMultEstimatorAvg[8]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC16p(TProfile* hprof){
      if(fMultEstimatorAvg[9]) delete fMultEstimatorAvg[9];
      fMultEstimatorAvg[9]=new TProfile(*hprof);
      fYearNumber = 16;
    }
    void SetMultiplVsZProfileLHC17e(TProfile* hprof){
      if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
      fMultEstimatorAvg[0]=new TProfile(*hprof);
      fYearNumber = 17;
    }
    void SetMultiplVsZProfileLHC17f(TProfile* hprof){
      if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
      fMultEstimatorAvg[1]=new TProfile(*hprof);
      fYearNumber = 17;
    }
    void SetMultiplVsZProfileLHC17h(TProfile* hprof){
      if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
      fMultEstimatorAvg[2]=new TProfile(*hprof);
      fYearNumber = 17;
    }
    void SetMultiplVsZProfileLHC17i(TProfile* hprof){
      if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
      fMultEstimatorAvg[3]=new TProfile(*hprof);
      fYearNumber = 17;
    }
    void SetMultiplVsZProfileLHC17j(TProfile* hprof){
      if(fMultEstimatorAvg[4]) delete fMultEstimatorAvg[4];
      fMultEstimatorAvg[4]=new TProfile(*hprof);
      fYearNumber = 17;
    }
    void SetMultiplVsZProfileLHC17k(TProfile* hprof){
      if(fMultEstimatorAvg[5]) delete fMultEstimatorAvg[5];
      fMultEstimatorAvg[5]=new TProfile(*hprof);
      fYearNumber = 17;
    }
    void SetMultiplVsZProfileLHC17l(TProfile* hprof){
      if(fMultEstimatorAvg[6]) delete fMultEstimatorAvg[6];
      fMultEstimatorAvg[6]=new TProfile(*hprof);
      fYearNumber = 17;
    }
    void SetMultiplVsZProfileLHC17m(TProfile* hprof){
      if(fMultEstimatorAvg[7]) delete fMultEstimatorAvg[7];
      fMultEstimatorAvg[7]=new TProfile(*hprof);
      fYearNumber = 17;
    }
    void SetMultiplVsZProfileLHC17o(TProfile* hprof){
      if(fMultEstimatorAvg[8]) delete fMultEstimatorAvg[8];
      fMultEstimatorAvg[8]=new TProfile(*hprof);
      fYearNumber = 17;
    }
    void SetMultiplVsZProfileLHC17r(TProfile* hprof){
      if(fMultEstimatorAvg[9]) delete fMultEstimatorAvg[9];
      fMultEstimatorAvg[9]=new TProfile(*hprof);
      fYearNumber = 17;
    }
    void SetMultiplVsZProfileLHC18b(TProfile* hprof){
      if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
      fMultEstimatorAvg[0]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18d(TProfile* hprof){
      if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
      fMultEstimatorAvg[1]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18e(TProfile* hprof){
      if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
      fMultEstimatorAvg[2]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18f(TProfile* hprof){
      if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
      fMultEstimatorAvg[3]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18g(TProfile* hprof){
      if(fMultEstimatorAvg[4]) delete fMultEstimatorAvg[4];
      fMultEstimatorAvg[4]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18h(TProfile* hprof){
      if(fMultEstimatorAvg[5]) delete fMultEstimatorAvg[5];
      fMultEstimatorAvg[5]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18i(TProfile* hprof){
      if(fMultEstimatorAvg[6]) delete fMultEstimatorAvg[6];
      fMultEstimatorAvg[6]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18j(TProfile* hprof){
      if(fMultEstimatorAvg[7]) delete fMultEstimatorAvg[7];
      fMultEstimatorAvg[7]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18k(TProfile* hprof){
      if(fMultEstimatorAvg[8]) delete fMultEstimatorAvg[8];
      fMultEstimatorAvg[8]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18l(TProfile* hprof){
      if(fMultEstimatorAvg[9]) delete fMultEstimatorAvg[9];
      fMultEstimatorAvg[9]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18m(TProfile* hprof){
      if(fMultEstimatorAvg[10]) delete fMultEstimatorAvg[10];
      fMultEstimatorAvg[10]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18n(TProfile* hprof){
      if(fMultEstimatorAvg[11]) delete fMultEstimatorAvg[11];
      fMultEstimatorAvg[11]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18o(TProfile* hprof){
      if(fMultEstimatorAvg[12]) delete fMultEstimatorAvg[12];
      fMultEstimatorAvg[12]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    void SetMultiplVsZProfileLHC18p(TProfile* hprof){
      if(fMultEstimatorAvg[13]) delete fMultEstimatorAvg[13];
      fMultEstimatorAvg[13]=new TProfile(*hprof);
      fYearNumber = 18;
    }
    
  
  void SetBDTPtCut(Double_t min, Double_t max) {fBDTPtCut[0]=min; fBDTPtCut[1]=max;}
  void SetBDTPtbins(AliRDHFCutsD0toKpiBDT* cut) {fCut4BDTptbin=cut;}
  void SetBDTSidebandSamplingFraction(Double_t f) {fBDTSidebandSamplingFraction=f;}
  void SetBDTSampleSideband(Bool_t sb) {fSampleSideband = sb;}
  void SetBDTFullVarString(TString str) {fBDTFullVarString = str;}
  void SetBDTClassifierVarString(TString str) {fBDTClassifierVarString = str;}
  
  void SetBDTList(TList *bdtlist) {fListRDHFBDT=bdtlist;}
  void SetBDTNamesList(TList *namelist) {fListBDTNames=namelist;}


  void SetEnableCentralityCorrCutsPbPb(Bool_t flag=kFALSE, Int_t year=2018) {
    fEnableCentralityCorrCuts=flag;
    if(year==2018){	
      fEventCuts.SetupPbPb2018();
      fEventCuts.SetManualMode();
    }else{
      fEventCuts.SetupRun2PbPb();
    }
  }


  Bool_t GetCutOnDistr() const {return fCutOnDistr;}
  Bool_t GetUsePid4Distr() const {return fUsePid4Distr;}
  Int_t  GetFillOnlyD0D0bar() const {return fFillOnlyD0D0bar;}
  Bool_t GetFillVarHists() const {return fFillVarHists;}
  Bool_t GetFillPtHistos() const {return fFillPtHist;}
  Bool_t GetFillYHistos() const {return fFillYHist;}
  Bool_t GetFillImpactParameterHistos() const {return fFillImpParHist;}
  Int_t  GetSystem() const {return fSys;}
  Bool_t GetRejectSDDClusters() const { return fIsRejectSDDClusters; }
  Bool_t GetUseSelectionBit() const { return fUseSelectionBit; }
  Bool_t GetWriteVariableTree() const {return fWriteVariableTree;}
  Bool_t GetDrawDetSignal() const {return fDrawDetSignal;}
  Bool_t GetPIDCheck() const {return fPIDCheck;}
  Bool_t GetFillSubSampleHist() const {return fFillSubSampleHist;}

 private:

  AliAnalysisTaskSED0BDT(const AliAnalysisTaskSED0BDT &source);
  AliAnalysisTaskSED0BDT& operator=(const AliAnalysisTaskSED0BDT& source);
  void	   DrawDetSignal(AliAODRecoDecayHF2Prong *part, TList *ListDetSignal);

  void     FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliAODMCHeader *mcHeader, AliRDHFCutsD0toKpiBDT *cuts, TList *listout);
  void     FillVarHists(AliAODEvent *aodev,AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpiBDT *cuts, TList *listout);
  void     FillCandVariables(AliAODEvent *aodev, AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliAODMCHeader *mcHeader, AliRDHFCutsD0toKpiBDT *cuts);
  AliAODVertex* GetPrimaryVtxSkipped(AliAODEvent *aodev);
  void CreateImpactParameterHistos();
  Int_t CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;
  Float_t GetTrueImpactParameter(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partD0) const ;
  Float_t ComputeTopomatic(AliAODEvent *aodev,AliAODRecoDecayHF2Prong *part);
    
    Bool_t fSubtractTrackletsFromDau; /// flag for subtracting D meson daughter contribution to N of tracklets
    Bool_t fmultiana;   ///
    Double_t fRefMult;   /// refrence multiplcity (period b)
    Int_t fMultiplicityEstimator; /// Definition of the multiplicity estimator: kNtrk10=0, kNtrk10to16=1, kVZERO=2
    Int_t fMCPrimariesEstimator;  /// Definition of the primaries estimator eta range: |eta|<1.0=0, -1.6<|eta|<1.0=1, VZEROrange=2
    TProfile* fMultEstimatorAvg[14]; /// TProfile with mult vs. Z per period.q
    Int_t fYearNumber; ///year number of the data taking
    TList  *fListProfiles; ///list of profile histos for z-vtx correction
    Int_t fDoVZER0ParamVertexCorr; /// Flag to use the zvtx correction from (0=none, 1=usual d2h, 2=AliESDUtils for VZERO multiplicity)
    TProfile* GetEstimatorHistogram(const AliVEvent *event);
    AliNormalizationCounter *fCounterC;           //!<!Counter for normalization, corrected multiplicity
    TH1F* fHistNtrCorrEvSel; //!<! hist. of ntracklets for selected events

  TList    *fOutputMass;          //!<! list send on output slot 1
  TList    *fOutputMassPt;        //!<! list send on output slot 6
  TList    *fOutputMassY;        //!<! list send on output slot 9
  TList    *fDistr;               //!<! list send on output slot 2
  TH1F     *fNentries;            //!<! histogram with number of events on output slot 3
  THnSparseF *fMCAccPrompt;       //!<!histo for StepMCAcc for D0 prompt (pt,y,ptB)
  THnSparseF *fMCAccBFeed;        //!<!histo for StepMCAcc for D0 FD (pt,y,ptB)
  Bool_t fStepMCAcc;              // flag to activate histos for StepMCAcc
  AliRDHFCutsD0toKpiBDT *fCuts;      //  Cuts - sent to output slot 4
  Bool_t    fEnableCentralityCorrCuts; /// flag to enable centrality correlation event cuts
  AliEventCuts  fEventCuts;       // Event cut object for centrality correlation event cuts
  THnSparseF *fHistMassPtImpParTC[5];   //!<! histograms for impact paramter studies
  Int_t     fArray;               ///  can be D0 or Like Sign candidates
  Bool_t    fReadMC;              ///  flag for MC array: kTRUE = read it, kFALSE = do not read it
  Bool_t    fCutOnDistr;          ///  flag to decide if apply cut also on distributions: 0 no cuts, 1 looser cuts, 2 tighter/ cuts
  Bool_t    fUsePid4Distr;        //  flag to use the particle identification to fill the signal histograms of distributions. It has effect only with fReadMC=kFALSE
  AliNormalizationCounter *fCounter;//!<! AliNormalizationCounter on output slot 5
  Int_t     fNPtBins;             ///  number of pt bins
  Double_t  fLsNormalization;     ///  normalization
  Int_t     fFillOnlyD0D0bar;     /// flag to fill mass histogram with D0/D0bar only (0 = fill with both, 1 = fill with D0 only, 2 = fill with D0bar only)
  TObjArray fDaughterTracks;      /// keeps the daughter tracks
  Int_t     fIsSelectedCandidate; /// selection outcome
  Bool_t    fFillVarHists;        /// flag to enable filling variable histos
  Int_t     fSys;                 /// fSys=0 -> p-p; fSys=1 ->PbPb (in this case fFillVarHists=kFALSE by default: set it to kTRUE *after* if needed)
  Bool_t    fIsRejectSDDClusters; /// flag to reject events with SDD clusters
  Bool_t    fFillPtHist;          /// flag to fill Pt and Impact Parameter Histograms
  Bool_t    fFillYHist;          /// flag to fill Y Histograms
  Bool_t    fFillImpParHist;      /// flag to fill Pt and Impact Parameter Histograms
  Bool_t    fFillSubSampleHist;    /// flag to fill SubSample histogram
  Int_t     fEventCounter; /// event counter used for sub sample test
  Bool_t    fUseSelectionBit;     /// flag to check or not the selection bit
  Int_t     fAODProtection;       /// flag to activate protection against AOD-dAOD mismatch.
                                  /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names

  Bool_t    fWriteVariableTree;       /// flag to decide whether to write the candidate variables on a tree variables
  TTree    *fVariablesTree;           //!<! tree of the candidate variables after track selection on output slot 7
  Double_t *fCandidateVariables;      //!<!  variables to be written to the tree
  Bool_t    fWriteProtosgnVar;        /// flag to decide whether to write the selected candidates variables on a tree for cut optimization
  Bool_t    fSelectTrueD0;            /// flag to decide whether to write only true D0/D0bar
  Bool_t    fUsedMassWindow;          /// flag to activate the mass window selection for output size reduction
  Bool_t	  fPIDCheck;			/// flag to decide whether to fill "PID = x" bins in fNentrie
  Bool_t    fDrawDetSignal;		/// flag to decide whether to draw the TPC dE/dx and TOF signal before/after PID
  Bool_t    fUseQuarkTagInKine;            // flag for quark/hadron level identification of prompt and feeddown
  Bool_t    fFillSparses;                  // flag to activate THnSparse
  Bool_t    fUseRejectionMethod;           // flag to activate the Rejection method
  Float_t   fRejectionFactor;              // rejection factor to be used in the rejection method
  THnSparseF *fhStudyImpParSingleTrackSign; //!<! sparse with imp par residual cuts for MC
  THnSparseF *fhStudyImpParSingleTrackCand;  //!<! sparse with imp par residual cuts for Data
  THnSparseF *fhStudyImpParSingleTrackFd;   //!<! sparse with imp par residual cuts for MC
  TList	   *fDetSignal;		//!<!Detector signal histograms (on output slot 8)
  TH2F *fhMultVZEROTPCoutTrackCorrNoCut;  //!<!
  TH2F *fhMultVZEROTPCoutTrackCorr;  //!<!
  Bool_t    fEnablePileupRejVZEROTPCout;
  TH2F *fhMultVZEROTPCclustersCorrNoCut;  //!<!
  TH2F *fhMultVZEROTPCclustersCorr;  //!<!
  Bool_t    fEnablePileupRejVZEROTPCcls;
  Bool_t    fRejectOutOfBunchPileUp;
  
  TH3F *h3Invmass[8];     //!<!
    TH3F *h3Invmass_19[8];     //!<!
    TH3F *h3Invmass_1029[8];     //!<!
    TH3F *h3Invmass_3059[8];     //!<!
    TH3F *h3Invmass_19999[8];     //!<!
    TH3F *h3Invmass_6099[8];     //!<!
  // = 																   =
  AliRDHFCutsD0toKpiBDT *fCut4BDTptbin;
  																	
  TList			*fListRDHFBDT;
  TList			*fListBDTNames;
  TList 		*fListBDTNtuple;
  TList 		*fListBDTResp;
  
  Double_t 		fBDTPtCut[2];
  Double_t		fBDTSidebandSamplingFraction;
  
  Bool_t 		fSampleSideband;
  
  TString		fBDTFullVarString;
  TString		fBDTClassifierVarString;
  
  Int_t     fIsSelectedCandidateBDT; /// selection outcome

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSED0BDT,26); /// AliAnalysisTaskSE for D0->Kpi
  /// \endcond
};

#endif
