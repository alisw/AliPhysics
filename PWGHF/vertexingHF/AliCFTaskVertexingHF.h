#ifndef ALICFTASKVERTEXINGHF_H
#define ALICFTASKVERTEXINGHF_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-----------------------------------------------------------------------
/// \class Class for HF corrections as a function of many variables and step
/// \author Author : C. Zampolli, CERN
///	\author		D. Caffarri, Univ & INFN Padova caffarri@pd.infn.it
/// \brief Base class for HF Unfolding - agrelli@uu.nl
//-----------------------------------------------------------------------


#include "AliAnalysisTaskSE.h"
#include "AliCFVertexingHF2Prong.h"
#include "AliCFVertexingHF3Prong.h"
#include "AliCFVertexingHFLctoV0bachelor.h"
#include "AliCFVertexingHF.h"
#include <TH1F.h>
#include <TProfile.h>

class TH1I;
class TFile ;
class TClonesArray ;
class AliCFManager;
class AliAODRecoDecay;
class AliAODRecoDecayHF2Prong;
class AliAODMCParticle;
class THnSparse;
class TF1;
class AliRDHFCuts;
class AliCFVertexingHF2Prong;
class AliCFVertexingHF3Prong;

class AliCFTaskVertexingHF: public AliAnalysisTaskSE {
 public:

  enum {
    kStepGeneratedLimAcc = 0,
    kStepGenerated       = 1,
    kStepAcceptance      = 2,
    kStepVertex          = 3,
    kStepRefit           = 4,
    kStepReconstructed   = 5,
    kStepRecoAcceptance  = 6,
    kStepRecoITSClusters = 7,
    kStepRecoPPR         = 8,
    kStepRecoPID         = 9,
    kStepGenLimAccNoAcc  = 10
  };

  enum {
    kSnail = 0,    /// slow configuration, all variables
    kCheetah = 1,   /// fast configuration, only a subset of variables
    kFalcon = 2,   /// super fast configuration, only (pt,y,centrality)
    kESE = 3,   /// configuration with variables for ESE analysis (pt,y,centrality,q2,mult)
    kRT = 4       /// configuration with variables for RT analysis (pt,y,mult,rt)
  };

  enum {
    kAll = 0,   /// all decays (resonant + non-resonant)
    kNonResonant = 1, /// only non resonant
    kL1520 = 2,  /// Lc --> L(1520) + p
    kKstar = 3,  /// Lc --> K* + pi
    kDelta = 4   /// Lc --> Delta + K
  };

  enum { kNtrk10=0, kNtrk10to16=1, kVZERO=2 }; /// multiplicity estimators

  AliCFTaskVertexingHF();
  AliCFTaskVertexingHF(const Char_t* name, AliRDHFCuts* cuts, TF1* func = 0x0);
  AliCFTaskVertexingHF& operator= (const AliCFTaskVertexingHF& c);
  AliCFTaskVertexingHF(const AliCFTaskVertexingHF& c);
  virtual ~AliCFTaskVertexingHF();

  /// ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Init();
  void     LocalInit() {Init();}
  void     Terminate(Option_t *);

  /// UNFOLDING
  void     SetCorrelationMatrix(THnSparse* h) {fCorrelation=h;}
  void     SetAcceptanceUnf(Bool_t AcceptanceUnf) {fAcceptanceUnf = AcceptanceUnf;}
  Bool_t   GetAcceptanceUnf() const {return fAcceptanceUnf;}


  /// CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* io) {fCFManager = io;}   /// global correction manager
  AliCFManager * GetCFManager()                 {return fCFManager;} /// get corr manager

  /// Setters (and getters) for the config macro
  void    SetFillFromGenerated(Bool_t flag) {fFillFromGenerated = flag;}
  Bool_t  GetFillFromGenerated() const {return fFillFromGenerated;}
  void    SetDecayChannel (Int_t decayChannel) {fDecayChannel = decayChannel;}
  Int_t   GetDecayChannel () {return fDecayChannel;}
  void     SetUseWeight(Bool_t useWeight){fUseWeight=useWeight;}
  Bool_t   GetUseWeight() const {return fUseWeight;}
  Double_t GetWeight(Float_t pt);
  Double_t dNdptFit(Float_t pt, Double_t* par);
  Double_t GetPtWeightFromHistogram(Float_t pt);

  void SetUseFlatPtWeight(Bool_t useWeight){fUseFlatPtWeight=useWeight; fUseWeight=useWeight;}
  Bool_t GetUseFlatPtWeight() const {return fUseFlatPtWeight;}
  void SetUseZWeight(Bool_t useWeight){fUseZWeight=useWeight;}
  Bool_t GetUseZWeight() const {return fUseZWeight;}
  Double_t GetZWeight(Float_t z, Int_t runnumber);
  Double_t DodzFit(Float_t z, Double_t* par);

  void SetUseNchWeight(Bool_t useWeight){fUseNchWeight=useWeight;}
  Bool_t GetUseNchWeight() const {return fUseNchWeight;}
  void SetMCNchHisto(TH1F* h){
    if(fHistoMCNch) delete fHistoMCNch;
    fHistoMCNch=new TH1F(*h);
  }
  void CreateMeasuredNchHisto();
  void SetMeasuredNchHisto(TH1F* h){
    if(fHistoMeasNch) delete fHistoMeasNch;
    fHistoMeasNch=new TH1F(*h);
  }
  Double_t GetNchWeight(Int_t nch);
  void SetMultiplicityEstimator(Int_t value){ fMultiplicityEstimator=value; }
  Int_t GetMultiplicityEstimator(){ return fMultiplicityEstimator; }
  void SetIsPPData(Bool_t flag){ fIsPPData = flag; }
  void SetIsPPbData(Bool_t flag){ fIsPPbData = flag; }
  void SetIsPP13TeVData(Bool_t flag){ fIsPP13TeVData = flag; }

  void SetUseNchTrackletsWeight(Bool_t useWeight = kTRUE) { fUseNchWeight=useWeight; fUseTrackletsWeight=useWeight; fUseMultRatioAsWeight=useWeight; }
  Bool_t GetUseNchTrackletsWeight() const {return fUseTrackletsWeight;}
  void SetUseRatioMultiplicityDistributionsAsWeight(Bool_t flag=kTRUE){ fUseMultRatioAsWeight=flag; }
  Bool_t GetUseRatioMultiplicityDistributionsAsWeight() const {return fUseMultRatioAsWeight;}

  void SetUseZvtxCorrectedNtrkEstimator(Bool_t flag) { fZvtxCorrectedNtrkEstimator=flag; }
  Bool_t GetUseZvtxCorrectedNtrkEstimator() { return fZvtxCorrectedNtrkEstimator; }
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
  void SetMultiplVsZProfilePP13TeV(TProfile* hprof, Int_t index){
    if(fMultEstimatorAvg[index]) delete fMultEstimatorAvg[index];
    fMultEstimatorAvg[index]=new TProfile(*hprof);
  }

  TProfile* GetEstimatorHistogram(const AliVEvent* event);
  void SetReferenceMultiplcity(Double_t rmu){fRefMult=rmu;}

  void   SetDselection(UShort_t originDselection) {fOriginDselection=originDselection;}
  UShort_t GetDselection (){return fOriginDselection;}
  void SetSign(Char_t isSign) {fSign = isSign;}
  Char_t GetSign() {return fSign;}

  void SetCentralitySelection(Bool_t centSelec = kTRUE) {fCentralitySelection = centSelec;}
  Bool_t GetCentralitySelection() {return fCentralitySelection;}

  void SetFakeSelection(Int_t fakeSel = 0) {fFakeSelection=fakeSel;}
  Int_t GetFakeSelection(){return fFakeSelection;}

  void SetRejectCandidateIfNotFromQuark(Bool_t opt){fRejectIfNoQuark=opt;}
  Bool_t GetRejectCandidateIfNotFromQuark(){return fRejectIfNoQuark;}

  void SetUseMCVertex(Bool_t opt){fUseMCVertex=opt;}
  Bool_t GetUseMCVertex(){return fUseMCVertex;}


  void SetKeepDsViaPhi(){fDsOption=1;}
  void SetKeepDsViaK0star(){fDsOption=2;}
  void SetKeepAllDs(){fDsOption=3;}
  void SetCountAllDs(){fGenDsOption=AliCFVertexingHF3Prong::kCountAllDsKKpi;}
  void SetCountDsViaPhi(){fGenDsOption=AliCFVertexingHF3Prong::kCountPhipi;}
  void SetCountDsViaK0star(){fGenDsOption=AliCFVertexingHF3Prong::kCountK0stK;}
  void SetCountResonantDs(){fGenDsOption=AliCFVertexingHF3Prong::kCountResonant;}
  void SetCountNonResonantDs(){fGenDsOption=AliCFVertexingHF3Prong::kCountNonResonant;}

  Bool_t ProcessDs(Int_t returnCodeDs) const;

  void SetConfiguration(Int_t configuration) {(configuration == kSnail) ? Printf("Slow configuration chosen, all variables will be used!") : Printf("Fast configuration chosen, all variablesOnly pt, y, phi, ct, fake, z_vtx, centrality and multiplicity will be used!"); fConfiguration = configuration;}
  Int_t GetConfiguration() const {return fConfiguration;}

  void SetWeightFunction(TF1* func) {fFuncWeight = func;}
  TF1* GetWeightFunction() const {return fFuncWeight;}
  void SetWeightHistogram(TH1F* histo) {
    if(fHistoPtWeight) delete fHistoPtWeight;
    fHistoPtWeight=new TH1F(*histo);
  }
  TH1F* GetWeightHistogram() const {return (TH1F*)fHistoPtWeight;}

  void SetPtWeightsFromFONLL276overLHC12a17a();
  void SetPtWeightsFromDataPbPb276overLHC12a17a();
  void SetPtWeightsFromFONLL276overLHC12a17b();
  void SetPtWeightsFromFONLL276andBAMPSoverLHC12a17b();
  void SetPtWeightsFromFONLL276overLHC10f6a();
  void SetPtWeightsFromFONLL7overLHC10f6a();
  void SetPtWeightsFromFONLL7overLHC12a12();
  void SetPtWeightsFromFONLL7overLHC12a12bis();
  void SetPtWeightsFromFONLL7overLHC13e2fix();
  void SetPtWeightsFromFONLL5overLHC10f6a();
  void SetPtWeightsFromFONLL5overLHC13d3();
  void SetPtWeightsFromFONLL5overLHC13d3Lc();
  void SetPtWeightsFromFONLL7overLHC11b2Lc();
  void SetPtWeightsFromFONLL7overLHC10f7aLc();
  void SetPtWeightsFromFONLL8overLHC15l2a2();
  void SetPtWeightsFromFONLL13overLHC17c3a12();

  void SetPtWeightsFDFromFONLL5anddataoverLHC19c3a();
  void SetPtWeightsFDFromFONLL5anddataoverLHC19c3b();
  void SetPtWeightsFromFONLL5anddataoverLHC16i2a();
  void SetPtWeightsFromFONLL5andLBToverLHC16i2a();
  void SetPtWeightsFromFONLL5overLHC16i2abc();
  void SetPtWeightsFromFONLL5andBAMPSoverLHC16i2abc();
  void SetPtWeightsFromFONLL5andTAMUoverLHC16i2abc();
  void SetPtWeightsFromD0Cent080dataoverLHC16i2abc();
  void SetPtWeightsFromD0Cent080dataModOhoverLHC16i2abc();
  void SetPtWeightsFromD0Cent080dataModMartinezoverLHC16i2abc();
  void SetPtWeightsFromFONLL5overLHC16i6a();
  void SetPtWeightsFromFONLL5andDplusdataoverLHC16i2a();
  void SetPtWeightsFromFONLL5overLHC18a4a2();
  void SetPtWeightsFromFONLL5anddataoverLHC20g2a();
  void SetPtWeightsFromFONLL5anddataoverLHC20g2b();

  void SetResonantDecay(UInt_t resonantDecay) {fResonantDecay = resonantDecay;}
  UInt_t GetResonantDecay() const {return fResonantDecay;}

  void SetKeepLctoK0Sp() {fLctoV0bachelorOption=1;}
  void SetKeepLctoLambdaBarpi() {fLctoV0bachelorOption=2;}
  void SetKeepLctoLambdapi() {fLctoV0bachelorOption=4;}
  void SetKeepLctoV0bachelor() {fLctoV0bachelorOption=7;}

  void SetCountLctoK0Sp(){fGenLctoV0bachelorOption=AliCFVertexingHFLctoV0bachelor::kCountK0Sp;}
  void SetCountLctoLambdapi(){fGenLctoV0bachelorOption=AliCFVertexingHFLctoV0bachelor::kCountLambdapi;}

  void SetUseSelectionBit(Bool_t flag) { fUseSelectionBit=flag; }
  Bool_t GetUseSelectionBit() const { return fUseSelectionBit; }

  Bool_t ProcessLctoV0Bachelor(Int_t returnCodeDs) const;

  void SetUseAdditionalCuts(Bool_t flag) { fDecayChannel == 22 ? fUseAdditionalCuts = flag : fUseAdditionalCuts = kFALSE;}
  Bool_t GetUseAdditionalCuts() const { return fUseAdditionalCuts; }
  void SetUseCutsForTMVA(Bool_t useCutsForTMVA) { fDecayChannel == 22 ? fUseCutsForTMVA = useCutsForTMVA : fUseAdditionalCuts = kFALSE;}
  Bool_t GetUseCutsForTMVA() const {return fUseCutsForTMVA;}

  void SetUseCascadeTaskForLctoV0bachelor(Bool_t useCascadeTaskForLctoV0bachelor) {fUseCascadeTaskForLctoV0bachelor = useCascadeTaskForLctoV0bachelor;}
  Bool_t GetUseCascadeTaskForLctoV0bachelor() const {return fUseCascadeTaskForLctoV0bachelor;}

  void SetFillMinimumSteps(Bool_t FillMinimumSteps) {fFillMinimumSteps = FillMinimumSteps;}
  Bool_t GetFillMinimumSteps() const {return fFillMinimumSteps;}

  void SetCutOnMomConservation(Float_t cut) {fCutOnMomConservation = cut;}
  Float_t GetCutOnMomConservation() const {return fCutOnMomConservation;}

  Double_t ComputeTPCq2(AliAODEvent* aod, AliAODMCHeader* mcHeader, Double_t etamin, Double_t etamax, Double_t ptmin, Double_t ptmax) const;

  Double_t CalculateRTValue(AliAODEvent* esdEvent, AliAODMCHeader* mcHeader, AliCFVertexingHF* cf);
  ULong64_t  GetEventIdAsLong(AliVHeader* header);
  TObjArray* FindLeading(TObjArray* array);
  void QSortTracks(TObjArray& a, Int_t first, Int_t last);
  TObjArray* SortRegionsRT(const AliVParticle* leading, TObjArray *array);
  TObjArray* GetMinMaxRegionRT(TList *transv1, TList *transv2);

  Double_t GetMinLeadPtRT() const {return fMinLeadPtRT;}
  void SetMinLeadPtRT(Double_t opt) {fMinLeadPtRT = opt;}

  void SetAODMismatchProtection(Int_t opt=0) {fAODProtection=opt;}
  void SetRejectOOBPileupEvents() {fRejectOOBPileUpEvents=kTRUE; fKeepOnlyOOBPileupEvents=kFALSE;}
  void SetKeepOnlyOOBPileupEvents() {fRejectOOBPileUpEvents=kFALSE; fKeepOnlyOOBPileupEvents=kTRUE;}

 protected:
  AliCFManager   *fCFManager;   ///  pointer to the CF manager
  TH1I *fHistEventsProcessed;   //!<! simple histo for monitoring the number of events processed
  THnSparse* fCorrelation;      ///  response matrix for unfolding
  TList  *fListProfiles; //list of profile histos for z-vtx correction
  Int_t fCountMC;               ///  MC particle found
  Int_t fCountGenLimAcc;        ///  MC particle found in limited acceptance
  Int_t fCountGenLimAccNoAcc;   ///  MC particle found in limited acceptance that doesn't satisfy acceptance cuts
  Int_t fCountAcc;              ///  MC particle found that satisfy acceptance cuts
  Int_t fCountVertex;       ///  Reco particle found that satisfy vertex constrained
  Int_t fCountRefit;        ///  Reco particle found that satisfy kTPCrefit and kITSrefit
  Int_t fCountReco;             ///  Reco particle found that satisfy cuts
  Int_t fCountRecoAcc;          ///  Reco particle found that satisfy cuts in requested acceptance
  Int_t fCountRecoITSClusters;  ///  Reco particle found that satisfy cuts in n. of ITS clusters
  Int_t fCountRecoPPR;          ///  Reco particle found that satisfy cuts in PPR
  Int_t fCountRecoPID;          /// Reco PID step
  Int_t fEvents;                ///  n. of events
  Int_t fDecayChannel;          /// decay channel to configure the task
  Bool_t fFillFromGenerated;    ///  flag to indicate whether data container should be filled with generated values also for reconstructed particles
  UShort_t fOriginDselection;      /// flag to select D0 origins. 0 Only from charm 1 only from beauty 2 both from charm and beauty
  Bool_t fAcceptanceUnf;        ///  flag for unfolding before or after cuts.
  AliRDHFCuts* fCuts;            /// cuts
  Bool_t fUseWeight;             /// flag to decide whether to use pt-weights != 1 when filling the container or not
  Double_t fWeight;              /// weight used to fill the container
  Bool_t fUseFlatPtWeight;       /// flag to decide to use a flat pt shape
  Bool_t fUseZWeight;           /// flag to decide whether to use z-vtx weights != 1 when filling the container or not
  Bool_t fUseNchWeight;         /// flag to decide whether to use Ncharged weights != 1 when filling the container or not
  Bool_t fUseTrackletsWeight;   /// flag to decide whether to use Ncharged weights != 1 when filling the container or not
  Bool_t fUseMultRatioAsWeight; /// flag to use directly the ratio of the distributions (fHistoMCNch) instead of computing it
  Int_t fNvar;                  /// number of variables for the container
  TString fPartName;    /// D meson name
  TString fDauNames;    /// daughter in fin state
  Char_t fSign;                 /// flag to decide wheter to keep D0 only (0), D0bar only (1), or both D0 and D0bar (2)
  Bool_t fCentralitySelection;  /// flag to switch off the centrality selection
  Int_t  fFakeSelection;  /// selection flag for fakes tracks
  Bool_t fRejectIfNoQuark;  /// flag to remove events not geenrated with PYTHIA
  Bool_t fUseMCVertex;  /// flag to use MC vertex (useful when runnign in pp)
  Int_t  fDsOption;     /// Ds decay option (selection level)
  Int_t  fGenDsOption;     /// Ds decay option (generation level)
  Int_t fConfiguration; /// configuration (slow / fast) of the CF --> different variables will be allocated (all / reduced number)
  TF1* fFuncWeight;     /// user-defined function to be used to calculate weights
  TH1F* fHistoPtWeight; /// user-defined histogram to calculate the Pt weights
  TH1F* fHistoMeasNch;  /// histogram with measured Nch distribution (pp 7 TeV)
  TH1F* fHistoMCNch;  /// histogram with Nch distribution from MC production
  UInt_t fResonantDecay;  /// resonant deacy channel to be used if the CF should be run on resonant channels only
  Int_t fLctoV0bachelorOption; /// Lc->V0+bachelor decay option (selection level)
  Int_t fGenLctoV0bachelorOption; /// Lc->V0+bachelor decay option (generation level)
  Bool_t fUseSelectionBit;     /// flag to use selection bit
  UInt_t fPDGcode; /// PDG code

  Int_t fMultiplicityEstimator; /// Definition of the multiplicity estimator: kNtrk10=0, kNtrk10to16=1, kVZERO=2
  TProfile* fMultEstimatorAvg[33]; /// TProfile with mult vas. Z per period
  Double_t fRefMult;   /// refrence multiplcity (period b)
  Bool_t fZvtxCorrectedNtrkEstimator; /// flag to use the z-vtx corrected (if not use uncorrected) multiplicity estimator
  Bool_t fIsPPData; /// flag for pp data (not checking centrality)
  Bool_t fIsPPbData; /// flag for pPb data (used for multiplicity corrections)
  Bool_t fIsPP13TeVData; /// flag for pp 13 TeV data (used for multiplicity corrections)
  Bool_t fUseAdditionalCuts;  /// flag to use additional cuts needed for Lc --> K0S + p, TMVA
  Bool_t fUseCutsForTMVA;     /// flag to use additional cuts needed for Lc --> K0S + p, TMVA
  /// these are the pre-selection cuts for the TMVA
  Bool_t fUseCascadeTaskForLctoV0bachelor;   /// flag to define which task to use for Lc --> K0S+p
  Bool_t fFillMinimumSteps;   /// Skip filling the unneed steps for most of the analyses to save disk space
  Float_t fCutOnMomConservation; /// cut on momentum conservation
  Double_t fMinLeadPtRT;   /// minimum pT cut for leading particle in RT calculation
  Int_t fAODProtection;         /// flag to activate protection against AOD-dAOD mismatch.
                                /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
  Bool_t fRejectOOBPileUpEvents; /// flag to enable rejection of events with simulated pileup
  Bool_t fKeepOnlyOOBPileupEvents; /// flag to keep only events with simulated pileup

  /// \cond CLASSIMP
  ClassDef(AliCFTaskVertexingHF,30); /// class for HF corrections as a function of many variables
  /// \endcond
};

#endif
