#ifndef ALIANALYSISTASKSEDVSMULTIPLICITY_BDT_H
#define ALIANALYSISTASKSEDVSMULTIPLICITY_BDT_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
/// \class Class AliAnalysisTaskSEDvsMultiplicity_BDT
/// \brief AliAnalysisTaskSE for the D meson vs. multiplcity analysis
/// \author Authors: Renu Bala, Zaida Conesa del Valle, Francesco Prino
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
#include "AliAODMCParticle.h"
#include "AliVertexingHFUtils.h"
#include "AliVEvent.h"
#include "AliRDHFBDT.h"

class AliAnalysisTaskSEDvsMultiplicity_BDT : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEDvsMultiplicity_BDT();
  AliAnalysisTaskSEDvsMultiplicity_BDT(const char *name, Int_t pdgMeson, AliRDHFCuts* cuts, Bool_t switchPPb, Bool_t readMC, Bool_t applyBDT);
  virtual ~AliAnalysisTaskSEDvsMultiplicity_BDT();


  void SetMassLimits(Double_t lowlimit, Double_t uplimit);
  void SetMassLimits(Int_t pdg, Double_t range);
  Double_t GetUpperMassLimit() const {return fUpmasslimit;}
  Double_t GetLowerMassLimit() const {return fLowmasslimit;}
  void SetNMassBins(Int_t nbins){fNMassBins=nbins;}
  Int_t GetNMassBins() const{return fNMassBins;}
  Bool_t GetSubtractTrackletsFromDaughters() const {return fSubtractTrackletsFromDau;}

  void SetImpactParameterBinning(Int_t nbins, Double_t dmin, Double_t dmax){
    fNImpParBins=nbins;
    fLowerImpPar=dmin;
    fHigherImpPar=dmax;
  }

  void SetMCOption(Int_t option=0){ fMCOption = option; }
  void SetIsPPbData(Bool_t flag=kTRUE){ 
    fisPPbData=flag;
  }
  void SetUseBit(Bool_t use=kTRUE){fUseBit=use;}
  void SetDoImpactParameterHistos(Bool_t doImp=kTRUE){fDoImpPar=doImp;}
  void SetKeepEstimatorCorrelationPlots(Bool_t use=kTRUE){fKeepCorrPlots=use;}

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
  void SetReferenceMultiplcity(Double_t rmu){fRefMult=rmu;}

  /// Nch Ntrk weights on MC
  void UseMCNchWeight(Int_t flag) { fUseNchWeight = flag; }
  void SetHistoNchWeight(TH1F *h){
    if(fHistoMCNch) delete fHistoMCNch;
    fHistoMCNch = new TH1F(*h);
  }
  void SetMeasuredNchHisto(TH1F* h){
    if(fHistoMeasNch) delete fHistoMeasNch;
    fHistoMeasNch = new TH1F(*h);
  }

  void SetSubtractTrackletsFromDaughters(Bool_t opt){fSubtractTrackletsFromDau=opt;}
  Int_t CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;

  void SetLcToV0decay(Bool_t flag) {fLctoV0=flag;  }

  /// Flag to use the zvtx correction from ( 0= none, 1= usual d2h, 2=AliESDUtils for VZERO multiplicity)
  void SetUseVZEROParameterizedVertexCorr(Int_t flag) { fDoVZER0ParamVertexCorr=flag; }
  Int_t GetUseVZEROParameterizedVertexCorr() { return fDoVZER0ParamVertexCorr; }

  enum { kNtrk10=0, kNtrk10to16=1, kVZERO=2, kNtrk03=3, kNtrk05=4, kVZEROA=5, kVZEROEq=6, kVZEROAEq=7 };
  void SetMultiplicityEstimator(Int_t value){ fMultiplicityEstimator=value; }
  Int_t GetMultiplicityEstimator(){ return fMultiplicityEstimator; }
  enum { kEta10=0, kEta10to16=1, kEtaVZERO=2, kEta03=3, kEta05=5, kEtaVZEROA=5 };
  void SetMCPrimariesEstimator(Int_t value){ fMCPrimariesEstimator=value; }
  Int_t GetMCPrimariesEstimator(){ return fMCPrimariesEstimator; }

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 
  //BDT application
  void SetBDTPtCut(Double_t min, Double_t max) {fBDTPtCut[0]=min; fBDTPtCut[1]=max;}
  void SetBDTRespCut(Double_t cut) {fBDTRespCut=cut;}
  void SetBDTSidebandCut(Double_t lcut, Double_t rcut) {fBDTSidebandCut[0]=lcut; fBDTSidebandCut[1]=rcut;}
  void SetBDTSidebandSamplingFraction(Double_t f) {fBDTSidebandSamplingFraction=f;}
  void SetBDTSampleSideband(Bool_t sb) {fSampleSideband = sb;}
  void SetBDTGetRespTree(Bool_t rt) {fGetRespTree = rt;}
  void SetBDTFullVarString(TString str) {fBDTFullVarString = str;}
  void SetBDTClassifierVarString(TString str) {fBDTClassifierVarString = str;}
  void SetFillTree(Bool_t fillTree) {fFillTree = fillTree;}
 
  void SetBDTList(TList *bdtlist) {fListRDHFBDT=bdtlist;} 
 private:

  AliAnalysisTaskSEDvsMultiplicity_BDT(const AliAnalysisTaskSEDvsMultiplicity_BDT &source);
  AliAnalysisTaskSEDvsMultiplicity_BDT& operator=(const AliAnalysisTaskSEDvsMultiplicity_BDT& source); 

  TProfile* GetEstimatorHistogram(const AliVEvent *event);
  void CreateImpactParameterHistos();
  void CreateMeasuredNchHisto();
  void FillMCMassHistos(TClonesArray *arrayMC, Int_t labD, Int_t countMult,Double_t nchWeight);

  void ProcessBDTD0(AliAODEvent *aod, AliAODRecoDecayHF *part, AliRDHFCuts *CutsAnalysis, Double_t bfield, TClonesArray *arrayMC, Int_t passCutsValue, Int_t coutMulti);
  void ProcessBDTDs(AliAODRecoDecayHF *part, AliRDHFCuts *CutsAnalysis, Double_t bfield, TClonesArray *arrayMC, Int_t passCutsValue, Int_t coutMulti);

  TList  *fOutput; //!<! list send on output slot 1
  TList  *fListCuts; ///list of cuts
  TList  *fOutputCounters; //!<! list send on output slot 3
  TList  *fListProfiles; ///list of profile histos for z-vtx correction
  Int_t     fFillOnlyD0D0bar;     /// flag to fill mass histogram with D0/D0bar only (0 = fill with both, 1 = fill with D0 only, 2 = fill with D0bar only)
  Bool_t    fUseQuarkTagInKine;            // flag for quark/hadron level identification of prompt and feeddown
 
  TH1F *fHistNEvents;     //!<!hist. for No. of events

  TH2F* fHistNtrEta16vsNtrEta1EvSel; //!<!hist. for Ntracklets in eta<1.6 vs. eta<1.
  TH2F* fHistNtrEta05vsNtrEta1EvSel; //!<!hist. for Ntracklets in eta<0.5 vs. eta<1.
  TH2F* fHistNtrEta03vsNtrEta1EvSel; //!<!hist. for Ntracklets in eta<0.3 vs. eta<1.
  TH2F* fHistNtrEtaV0AvsNtrEta1EvSel; //!<!hist. for Ntracklets in eta-V0A vs. eta<1.
  TH2F* fHistNtrEtaV0MvsNtrEta1EvSel; //!<!hist. for Ntracklets in eta-V0M vs. eta<1.
  TH2F* fHistNtrEtaV0AvsV0AEqEvSel;   //!<!hist. for V0A raw mult vs V0A equalized multiplicity
  TH2F* fHistNtrEtaV0MvsV0MEqEvSel;   //!<!hist. for V0M raw mult vs V0M equalized multiplicity
  TH2F* fHistNtrCorrEta1vsNtrRawEta1EvSel; //!<!hist. for Ntracklets in eta<1 with and w/o corrections
  TH2F* fHistMultCorrvsMultRawEvSel;       //!<!hist. for multiplicity with and w/o corrections
  TH2F* fHistNtrEta16vsNtrEta1EvWithCand; //!<!hist. for Ntracklets in eta<1.6 vs. eta<1. for events with a candidate
  TH2F* fHistNtrEta05vsNtrEta1EvWithCand; //!<!hist. for Ntracklets in eta<0.5 vs. eta<1. for events with a candidate
  TH2F* fHistNtrEta03vsNtrEta1EvWithCand; //!<!hist. for Ntracklets in eta<0.3 vs. eta<1. for events with a candidate
  TH2F* fHistNtrEtaV0AvsNtrEta1EvWithCand; //!<!hist. for Ntracklets in eta-V0A vs. eta<1. for events with a candidate
  TH2F* fHistNtrEtaV0MvsNtrEta1EvWithCand; //!<!hist. for Ntracklets in eta-V0M vs. eta<1. for events with a candidate
  TH2F* fHistNtrEtaV0AvsV0AEqEvWithCand;     //!<!hist. for V0A raw mult vs V0A equalized multiplicity for events with a candidate
  TH2F* fHistNtrEtaV0MvsV0MEqEvWithCand;     //!<!hist. for V0M raw mult vs V0M equalized multiplicity for events with a candidate
  TH2F* fHistNtrCorrEta1vsNtrRawEta1EvWithCand; //!<!hist. for Ntracklets in eta<1 with and w/o corrections for events with a candidate
  TH2F* fHistMultCorrvsMultRawEvWithCand;       //!<!hist. for multiplicity with and w/o corrections for events with a candidate
  TH2F* fHistNtrEta16vsNtrEta1EvWithD; //!<!hist. for Ntracklets in eta<1.6 vs. eta<1. for events with a candidate in D mass peak
  TH2F* fHistNtrEta05vsNtrEta1EvWithD; //!<!hist. for Ntracklets in eta<0.5 vs. eta<1. for events with a candidate in D mass peak
  TH2F* fHistNtrEta03vsNtrEta1EvWithD; //!<!hist. for Ntracklets in eta<0.3 vs. eta<1. for events with a candidate in D mass peak
  TH2F* fHistNtrEtaV0AvsNtrEta1EvWithD; //!<!hist. for Ntracklets in eta-V0A vs. eta<1. for events with a candidate in D mass peak
  TH2F* fHistNtrEtaV0MvsNtrEta1EvWithD; //!<!hist. for Ntracklets in eta-V0M vs. eta<1. for events with a candidate in D mass peak
  TH2F* fHistNtrEtaV0AvsV0AEqEvWithD;   //!<!hist. for V0A raw mult vs V0A equalized multiplicity with a candidate in D mass peak
  TH2F* fHistNtrEtaV0MvsV0MEqEvWithD;   //!<!hist. for V0M raw mult vs V0M equalized multiplicity with a candidate in D mass peak
  TH2F* fHistNtrCorrEta1vsNtrRawEta1EvWithD; //!<!hist. for Ntracklets in eta<1 with and w/o corrections for events with a candidate in D mass peak
  TH2F* fHistMultCorrvsMultRawEvWithD;       //!<!hist. for multiplicity with and w/o corrections for events with a candidate in D mass peak

  TH2F* fHistNtrVsZvtx; //!<!  hist of ntracklets vs Zvertex
  TH2F* fHistNtrCorrVsZvtx; //!<!  hist of ntracklets vs Zvertex

  TH2F* fHistNtrVsNchMC; //!<!  hist of ntracklets vs Nch (Generated)
  TH2F* fHistNtrCorrVsNchMC; //!<!  hist of ntracklets vs Nch (Generated)
  TH2F* fHistNtrVsNchMCPrimary; //!<!  hist of ntracklets vs Nch (Primary)
  TH2F* fHistNtrCorrVsNchMCPrimary; //!<!  hist of ntracklets vs Nch (Primary)
  TH2F* fHistNtrVsNchMCPhysicalPrimary; //!<!  hist of ntracklets vs Nch (Physical Primary)
  TH2F* fHistNtrCorrVsNchMCPhysicalPrimary; //!<!  hist of ntracklets vs Nch (Physical Primary)
  TH1F* fHistGenPrimaryParticlesInelGt0; //!<!hist. of geenrated multiplcity
  TH3F* fHistNchMCVsNchMCPrimaryVsNchMCPhysicalPrimary; //!<! hist of Nch (generated) vs Nch (Primary) vs Nch (Physical Primary) 
  
  TH1F* fHistNtrUnCorrPSSel; //!<! hist. of ntracklets for physics selection only selected events
  TH1F* fHistNtrUnCorrPSTrigSel; //!<! hist. of ntracklets for physics selection + trigger name selected events
  TH1F* fHistNtrUnCorrPSTrigPileUpSel; //!<! hist. of ntracklets for physics selection + trigger name + pileup selected events
  TH1F* fHistNtrUnCorrPSTrigPileUpVtxSel; //!<! hist. of ntracklets for physics selection + trigger name + pileup + with-vertex selected events
  TH1F* fHistNtrUnCorrPSTrigPileUpVtxContSel; //!<! hist. of ntracklets for physics selection + trigger name + pileup + with-vertex-contrib selected events
  TH1F* fHistNtrUnCorrPSTrigPileUpVtxRangeSel; //!<! hist. of ntracklets for physics selection + trigger name + pileup + with-vertex-contrib-range selected events
  TH1F* fHistNtrUnCorrPSTrigPileUpVtxRangeCentrSel; //!<! hist. of ntracklets for physics selection + trigger name + pileup + with-vertex-contrib-range + centrality selected events
  TH1F* fHistNtrUnCorrEvSel; //!<! hist. of ntracklets for selected events
  TH1F* fHistNtrUnCorrEvWithCand; //!<! hist. of ntracklets for evnts with a candidate
  TH1F* fHistNtrUnCorrEvWithD;//!<! hist. of ntracklets for evnts with a candidate in D mass peak
  TH1F* fHistNtrCorrPSSel; //!<! hist. of ntracklets for physics selection only selected events
  TH1F* fHistNtrCorrEvSel; //!<! hist. of ntracklets for selected events
  TH1F* fHistNtrCorrEvWithCand; //!<! hist. of ntracklets for evnts with a candidate
  TH1F* fHistNtrCorrEvWithD;//!<! hist. of ntracklets for evnts with a candidate in D mass peak


  TH3F *fPtVsMassVsMult;  //!<! hist. of Pt vs Mult vs. mass (
  TH3F *fPtVsMassVsMultNoPid;  //!<! hist. of Pt vs Mult vs. mass (no pid)
  TH3F *fPtVsMassVsMultUncorr;  //!<! hist. of Pt vs Mult vs. mass (raw mult)
  TH3F *fPtVsMassVsMultPart;  //!<! hist. of Pt vs Mult vs. mass (particle)
  TH3F *fPtVsMassVsMultAntiPart;  //!<! hist. of Pt vs Mult vs. mass (antiparticle)
  TH3F *fPtVsMassVsMultMC;  //!<! hist. of Pt vs Mult vs. mass (MC true candidates before reconstruction)

  THnSparseF *fHistMassPtImpPar[5];//!<! histograms for impact paramter studies
 
  Double_t fUpmasslimit;  /// upper inv mass limit for histos
  Double_t fLowmasslimit; /// lower inv mass limit for histos
  Int_t   fNMassBins;    /// nbins for invariant mass histos

  AliRDHFCuts *fRDCutsAnalysis; /// Cuts for Analysis
  AliNormalizationCounter *fCounter;            //!<! Counter for normalization
  AliNormalizationCounter *fCounterC;           //!<!Counter for normalization, corrected multiplicity
  AliNormalizationCounter *fCounterU;           //!<!Counter for normalization, uncorrected multiplicity
  AliNormalizationCounter *fCounterCandidates;  //!<!Counter for normalization, corrected multiplicity for candidates


  Bool_t fDoImpPar;  /// swicth for D impact parameter THnSparse
  Int_t  fNImpParBins;   /// nunber of bins in impact parameter histos
  Double_t fLowerImpPar;  /// lower limit in impact parameter (um)
  Double_t fHigherImpPar; /// higher limit in impact parameter (um)

  Bool_t fReadMC;    /// flag for access to MC
  Int_t  fMCOption;  /// 0=keep all cand, 1=keep only signal, 2= keep only back
  Bool_t fisPPbData; /// flag to run on pPb data (differen histogram bining)
  Bool_t fUseBit;    /// flag to use bitmask
  Bool_t fSubtractTrackletsFromDau; /// flag for subtracting D meson daughter contribution to N of tracklets
  Bool_t fKeepCorrPlots; /// flag to look at the correlation of different estimators (eta ranges)
  Int_t fAODProtection;  /// flag to activate protection against AOD-dAOD mismatch.
  /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names

  Int_t fUseNchWeight; /// weight on the MC on the generated multiplicity (0->no weights, 1->Nch weights, 2->Ntrk weights)
  TH1F* fHistoMCNch;    /// weight histogram for the MC on the generated multiplicity
  TH1F* fHistoMeasNch;  /// weight histogram on the true measured multiplicity
  
  TProfile* fMultEstimatorAvg[14]; /// TProfile with mult vs. Z per period
  Double_t fRefMult;   /// refrence multiplcity (period b)
  Int_t fPdgMeson;   /// pdg code of analyzed meson
  Bool_t fLctoV0;    /// flag for Lc in K0sp decay

  Int_t fMultiplicityEstimator; /// Definition of the multiplicity estimator: kNtrk10=0, kNtrk10to16=1, kVZERO=2
  Int_t fMCPrimariesEstimator;  /// Definition of the primaries estimator eta range: |eta|<1.0=0, -1.6<|eta|<1.0=1, VZEROrange=2

  Int_t fDoVZER0ParamVertexCorr; /// Flag to use the zvtx correction from (0=none, 1=usual d2h, 2=AliESDUtils for VZERO multiplicity)
  
  Int_t fYearNumber; ///year number of the data taking

  TObjArray fDaughterTracks;      /// keeps the daughter tracks
  
  //BDT application
  TList                 *fListRDHFBDT;
  TList                 *fListBDTNtuple;
  TList                 *fListBDTResp;

  Double_t              fBDTPtCut[2];
  Double_t              fBDTRespCut;
  Double_t              fBDTSidebandSamplingFraction;
  Double_t              fBDTSidebandCut[2];
  
  Bool_t                fSampleSideband;
  Bool_t                fGetRespTree;
  Bool_t                fApplyBDT; 
  Bool_t                fFillTree;

  TString               fBDTFullVarString;
  TString               fBDTClassifierVarString;
    
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEDvsMultiplicity_BDT,21); /// charmed hadrons vs. mult task
  /// \endcond
};

#endif
