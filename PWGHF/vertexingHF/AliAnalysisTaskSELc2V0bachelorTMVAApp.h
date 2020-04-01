#ifndef ALIANALYSISTASKSELC2V0BACHELORTMVAAPP_H
#define ALIANALYSISTASKSELC2V0BACHELORTMVAAPP_H
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

/* $Id: AliAnalysisTaskSELc2V0bachelorTMVAApp.h 61835 2013-04-05 23:07:23Z fprino $ */ 

#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliPID.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliTPCPIDResponse.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"
#include "AliAODRecoCascadeHF.h"

#include <TMVA/Tools.h>
#include <TMVA/Reader.h>
#include <TMVA/MethodCuts.h>
#include <TProfile.h>

/// \class AliAnalysisTaskSELc2V0bachelorTMVAApp

class IClassifierReader;
class ReadBDT_Default;

class TH1F;
class TH1D;

class AliAnalysisTaskSELc2V0bachelorTMVAApp : public AliAnalysisTaskSE 
{
  
 public:

  enum EBachelor {
    kBachInvalid = -1,
    kBachFake = 0,
    kBachNoProton = 1,
    kBachPrimary = 2,
    kBachNoLambdaMother = 3,
    kBachDifferentLambdaMother = 4,
    kBachCorrectLambdaMother = 5 };

  enum EK0S {
    kK0SInvalid = -1,
    kK0SFake = 0,
    kK0SNoK0S = 1,
    kK0SWithoutMother = 2,
    kK0SNotFromK0 = 3,
    kK0Primary = 4,
    kK0NoLambdaMother = 5,
    kK0DifferentLambdaMother = 6,
    kK0CorrectLambdaMother = 7 };    
  
  
  AliAnalysisTaskSELc2V0bachelorTMVAApp();
  AliAnalysisTaskSELc2V0bachelorTMVAApp(const Char_t* name, AliRDHFCutsLctoV0* cutsA,
				 Bool_t useOnTheFly=kFALSE);
  virtual ~AliAnalysisTaskSELc2V0bachelorTMVAApp();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 
  /// histos
  void FillLc2pK0Sspectrum(AliAODRecoCascadeHF *part, Int_t isLc,
			   Int_t &nSelectedAnal, AliRDHFCutsLctoV0 *cutsAnal,
			   TClonesArray *mcArray, Int_t iLctopK0s, AliAODEvent *aod);

  void MakeAnalysisForLc2prK0S(AliAODEvent *aodEvent,
			       TClonesArray *arrayLctopK0s,
			       TClonesArray *mcArray,
			       Int_t &nSelectedAnal, AliRDHFCutsLctoV0 *cutsAnal, 
			       TClonesArray *array3Prong, AliAODMCHeader *aodheader);
  
  void SetMVReader(IClassifierReader* r) {fBDTReader = r;}
  IClassifierReader* const GetMVReader() {return fBDTReader;}
  void SetTMVAlibName(const char* libName) {fTMVAlibName = libName;}
  TString GetTMVAlibName() {return fTMVAlibName;}
  void SetTMVAlibPtBin(const char* libPtBin) {fTMVAlibPtBin = libPtBin;}
  TString GetTMVAlibPtBin() {return fTMVAlibPtBin;}
  void SetNamesTMVAVariables(TString names) {fNamesTMVAVar = names;}
  TString GetNamesTMVAVariables() {return fNamesTMVAVar;}
  
  /// set MC usage
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}

  void SetK0sAnalysis(Bool_t a) {fIsK0sAnalysis=a;}
  Bool_t GetK0sAnalysis() const {return fIsK0sAnalysis;}

  void SetUseOnTheFlyV0(Bool_t a) { fUseOnTheFlyV0=a; }
  Bool_t GetUseOnTheFlyV0() { return fUseOnTheFlyV0; }

  void SetFillOnlySgn(Bool_t a) { fFillOnlySgn=a; }
  Bool_t GetFillOnlySgn() { return fFillOnlySgn; }

  void SetTopoConstraint(Bool_t a) { ftopoConstraint=a; }
  Bool_t GetTopoConstraint() { return ftopoConstraint; }

  void SetCallKFVertexing(Bool_t a) { fCallKFVertexing=a; }
  Bool_t GetCallKFVertexing() { return fCallKFVertexing; }

  void SetKeepingOnlyHIJINGBkg(Bool_t a) { fKeepingOnlyHIJINGBkg = a;}
  Bool_t GetKeepingOnlyHIJINGBkg() {return fKeepingOnlyHIJINGBkg;}

  void SetKFCutChi2NDF(Float_t a) {fCutKFChi2NDF = a;}
  Float_t GetKFCutChi2NDF() {return fCutKFChi2NDF;}

  void SetKFCutDeviationFromVtx(Float_t a) {fCutKFDeviationFromVtx = a;}
  Float_t GetKFCutDeviationFromVtx() {return fCutKFDeviationFromVtx;}

  void SetKFCutDeviationFromVtxV0(Float_t a) {fCutKFDeviationFromVtxV0 = a;}
  Float_t GetKFCutDeviationFromVtxV0() {return fCutKFDeviationFromVtxV0;}

  void SetKeepingOnlyPYTHIABkg(Bool_t a) { fKeepingOnlyPYTHIABkg = a;}
  Bool_t GetKeepingOnlyPYTHIABkg() {return fKeepingOnlyPYTHIABkg;}
  
  void SetTriggerMask(ULong64_t c) { fTriggerMask = c;}	
  
  void SetFillTree(Bool_t a) { fFillTree = a;}	
  
  void SetMCNchHisto(TH1F* h){
    if(fHistoMCNch) delete fHistoMCNch;
    fHistoMCNch = new TH1F(*h);
  }
    
  void SetDebugHistograms(Bool_t flag) {fDebugHistograms = flag;}
  Bool_t GetDebugHistograms() const {return fDebugHistograms;}

  void SetAODMismatchProtection(Int_t opt = 1) {fAODProtection = opt;}
  Int_t GetAODMismatchProtection() const {return fAODProtection;}

  void SetUsePIDresponseForNsigma(Bool_t flag) {fUsePIDresponseForNsigma = flag;}
  Bool_t GetUsePIDresponseForNsigma() const {return fUsePIDresponseForNsigma;}

  void SetNVars(Int_t n) {fNVars = n;}
  Int_t GetNVars() const {return fNVars;}

  void SetTimestampCut(UInt_t value) {fTimestampCut = value;}
  UInt_t GetTimestampCut() const {return fTimestampCut;}

  void SetTMVAReader(TMVA::Reader* r) {fReader = r;}
  TMVA::Reader* GetTMVAReader() const {return fReader;}

  void SetNVarsSpectators(Int_t n) {fNVarsSpectators = n;}
  Int_t GetNVarsSpectators() const {return fNVarsSpectators;}

  void SetNamesTMVAVariablesSpectators(TString names) {fNamesTMVAVarSpectators = names;}
  TString GetNamesTMVAVariablesSpectators() {return fNamesTMVAVarSpectators;}

  void SetUseXmlWeightsFile(Bool_t flag) {fUseXmlWeightsFile = flag;}
  Bool_t GetUseXmlWeightsFile() const {return fUseXmlWeightsFile;}

  void SetUseWeightsLibrary(Bool_t flag) {fUseWeightsLibrary = flag;}
  Bool_t GetUseWeightsLibrary() const {return fUseWeightsLibrary;}

  void SetXmlWeightsFile(TString fileName) {fXmlWeightsFile = fileName;}
  TString GetXmlWeightsFile() const {return fXmlWeightsFile;}

  void SetUseXmlFileFromCVMFS(Bool_t flag) {fUseXmlFileFromCVMFS = flag;}
  Bool_t GetUseXmlFileFromCVMFS() const {return fUseXmlFileFromCVMFS;}

  void SetXmlFileFromCVMFS(TString fileName) {fXmlFileFromCVMFS = fileName;}
  TString GetXmlFileFromCVMFS() const {return fXmlFileFromCVMFS;}

  void SetUseMultiplicityCorrection(Bool_t flag){fUseMultCorrection=flag;}

  void SetReferenceMultiplcity(Double_t rmu){fRefMult=rmu;}

  enum { kNtrk10=0, kNtrk10to16=1, kVZERO=2, kNtrk03=3, kNtrk05=4, kVZEROA=5, kVZEROEq=6, kVZEROAEq=7 };
  void SetMultiplicityEstimator(Int_t value){ fMultiplicityEstimator=value; }
  Int_t GetMultiplicityEstimator(){ return fMultiplicityEstimator; }

  /// Flag to use the zvtx correction from ( 0= none, 1= usual d2h, 2=AliESDUtils for VZERO multiplicity)
  void SetUseVZEROParameterizedVertexCorr(Int_t flag) { fDoVZER0ParamVertexCorr=flag; }
  Int_t GetUseVZEROParameterizedVertexCorr() { return fDoVZER0ParamVertexCorr; }

  void SetUseMultiplcityCut(Bool_t flag){fUseMultiplicityCut=flag;}

  void SetMinimumMultiplicity(Float_t value){fMultiplicityCutMin=value;}

  void SetMaximumMultiplicity(Float_t value){fMultiplicityCutMax=value;}

  // pp - 2010
  void SetMultiplVsZProfileLHC10b(TProfile* hprof){
    if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
    fMultEstimatorAvg[0]=new TProfile(*hprof);
    fYearNumber = 10;
  }
  void SetMultiplVsZProfileLHC10c(TProfile* hprof){
    if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
    fMultEstimatorAvg[1]=new TProfile(*hprof);
    fYearNumber = 10;
  }
  void SetMultiplVsZProfileLHC10d(TProfile* hprof){
    if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
    fMultEstimatorAvg[2]=new TProfile(*hprof);
    fYearNumber = 10;
  }
  void SetMultiplVsZProfileLHC10e(TProfile* hprof){
    if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
    fMultEstimatorAvg[3]=new TProfile(*hprof);
    fYearNumber = 10;
  }

  // pp - 2016
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

  // pp - 2017
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

  // pp - 2018
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

 private:
  
  EBachelor CheckBachelor(AliAODRecoCascadeHF *part, AliAODTrack* bachelor, TClonesArray *mcArray);
  EK0S CheckK0S(AliAODRecoCascadeHF *part, AliAODv0* v0part, TClonesArray *mcArray);
  Int_t FindV0Label(AliAODRecoDecay* v0part, TClonesArray *mcArray) const;
  Int_t FindLcLabel(AliAODRecoCascadeHF* cascade, TClonesArray *mcArray) const;
  Int_t CallKFVertexing(AliAODRecoCascadeHF *cascade, AliAODv0* v0part, AliAODTrack* bach, TClonesArray *mcArray, 
			Double_t* V0KF, Double_t* errV0KF, Double_t* LcKF, Double_t* errLcKF,
			Double_t* distances, Double_t* armPolKF);

  void FillMCHisto(TClonesArray *mcArray);

  AliAnalysisTaskSELc2V0bachelorTMVAApp(const AliAnalysisTaskSELc2V0bachelorTMVAApp &source);
  AliAnalysisTaskSELc2V0bachelorTMVAApp& operator=(const AliAnalysisTaskSELc2V0bachelorTMVAApp& source); 
  /* TProfile* GetEstimatorHistogram(const AliVEvent *event); */
  
  Bool_t fUseMCInfo;          /// Use MC info
  TList *fOutput;             //!<! User output1: list of trees

  // define the histograms
  AliPIDResponse *fPIDResponse;       //!<! PID response object
  AliPIDCombined *fPIDCombined;       //!<! combined PID response object
  Bool_t fIsK0sAnalysis;              /// switch between Lpi and K0sp
  AliNormalizationCounter *fCounter;  //!<! AliNormalizationCounter on output slot 4
  AliNormalizationCounter *fCounterC; //!<! AliNormalizationCounter on output slot 4, corrected with multiplicity dependence
  AliRDHFCutsLctoV0 *fAnalCuts;       /// Cuts - sent to output slot 5
  TList *fListCuts;                   //!<! list of cuts
  TList *fListWeight;                 /// list of weights
  TList *fListCounters;               //!<! list of counters on output slot 2
  TList *fListProfiles;               /// list of profiles for z-vtx correction of multiplicity
  Bool_t fUseOnTheFlyV0;              /// flag to analyze also on-the-fly V0 candidates
  Bool_t fIsEventSelected;            /// flag for event selected

  TTree   *fVariablesTreeSgn;         //!<! tree of the candidate variables after track selection (Signal)
  TTree   *fVariablesTreeBkg;         //!<! tree of the candidate variables after track selection (Background)
  Float_t *fCandidateVariables;       //!<! variables to be written to the tree

  TH1F* fHistoCentrality;             //!<! histogram with centrality from AliRDHFCuts
  TH1F* fHistoEvents;                 //!<! histogram with number of events analyzed
  TH1F* fHistoTracklets_1;            //!<! histogram with number of tracklets in the event in eta [-1, 1]
  TH2F* fHistoTracklets_1_cent;       //!<! histogram with number of tracklets in the event in eta [-1, 1] vs centrality
  TH1F* fHistoTracklets_All;          //!<! histogram with number of tracklets in the event in eta [-999, 999]
  TH2F* fHistoTracklets_All_cent;     //!<! histogram with number of tracklets in the event in eta [-999, 999] vs centrality
  TH1F* fHistoLc;                     //!<! histogram with number of Lc
  TH1F* fHistoLcOnTheFly;             //!<! histogram with number of Lc with on-the-fly V0
  Bool_t fFillOnlySgn;                /// flag to fill only signal (speeding up processing)
  TH1F* fHistoLcBeforeCuts;           //!<! histogram with number of Lc before any cut 
  TH1F* fHistoFiducialAcceptance;     //!<! histogram to check FiducialAcceptance cut
  TH2F* fHistoCodesSgn;               //!<! histogram with codes for bachelor and V0 for signal
  TH2F* fHistoCodesBkg;               //!<! histogram with codes for bachelor and V0 for background
  TH1F* fHistoLcpKpiBeforeCuts;       //!<! histogram number of true Lc-->pKpi (3 prong) before any cut
  AliAODVertex *fVtx1;                /// primary vertex

  TH1D* fHistoDistanceLcToPrimVtx;    //!<! KF: distance Lc vertex from primary vertex   
  TH1D* fHistoDistanceV0ToPrimVtx;    //!<! KF: distance V0 vertex from primary vertex   
  TH1D* fHistoDistanceV0ToLc;         //!<! KF: distance V0 vertex from Lc vertex    

  TH1D* fHistoDistanceLcToPrimVtxSgn; //!<! KF: distance of signal Lc vertex from primary vertex    
  TH1D* fHistoDistanceV0ToPrimVtxSgn; //!<! KF: distance for signal Lc of V0 vertex from primary vertex   
  TH1D* fHistoDistanceV0ToLcSgn;      //!<! KF: distance for signal Lc of V0 vertex from Lc vertex 
         
  TH1D* fHistoVtxLcResidualToPrimVtx; //!<! KF: residual wrt MC of distance Lc vertex from primary vertex (MC - KF)
  TH1D* fHistoVtxV0ResidualToPrimVtx; //!<! KF: residual wrt MC of distance V0 vertex from primary vertex (MC - KF)
  TH1D* fHistoVtxV0ResidualToLc;      //!<! KF: residual wrt MC of distance V0 vertex from Lc vertex (MC - KF)

  TH1D* fHistoMassV0All;              //!<! KF: mass for all V0 reconstructed with KF
  TH1D* fHistoDecayLengthV0All;       //!<! KF: decay length for all V0 reconstructed with KF
  TH1D* fHistoLifeTimeV0All;          //!<! KF: life time for all V0 reconstructed with KF

  TH1D* fHistoMassV0True;             //!<! KF: mass for true V0 reconstructed with KF
  TH1D* fHistoDecayLengthV0True;      //!<! KF: decay length for true V0 reconstructed with KF
  TH1D* fHistoLifeTimeV0True;         //!<! KF: life time for true V0 reconstructed with KF

  TH1D* fHistoMassV0TrueFromAOD;      //!<! KF: AOD mass for true V0 reconstructed with KF

  TH1D* fHistoMassV0TrueK0S;          //!<! KF: mass for true V0 which are really K0S reconstructed with KF
  TH1D* fHistoDecayLengthV0TrueK0S;   //!<! KF: decay length for true V0 which are really K0S reconstructed with KF
  TH1D* fHistoLifeTimeV0TrueK0S;      //!<! KF: life time for true V0 which are really K0S reconstructed with KF

  TH1D* fHistoMassV0TrueK0SFromAOD;   //!<! KF: AOD mass for true V0 which are really K0S reconstructed with KF

  TH1D* fHistoMassLcAll;              //!<! KF: mass for all Lc reconstructed with KF
  TH1D* fHistoDecayLengthLcAll;       //!<! KF: decay length for all Lc reconstructed with KF
  TH1D* fHistoLifeTimeLcAll;          //!<! KF: life time for all Lc reconstructed with KF

  TH1D* fHistoMassLcTrue;             //!<! KF: mass for true cascades reconstructed with KF
  TH1D* fHistoDecayLengthLcTrue;      //!<! KF: decay length for true cascades reconstructed with KF
  TH1D* fHistoLifeTimeLcTrue;         //!<! KF: life time for true cascades reconstructed with KF

  TH1D* fHistoMassLcTrueFromAOD;      //!<! KF: AOD mass for true cascades reconstructed with KF

  TH1D* fHistoMassV0fromLcAll;        //!<! KF: mass of V0 for all cascades reconstructed with KF
  TH1D* fHistoDecayLengthV0fromLcAll; //!<! KF: decay length of V0 for all cascades reconstructed with KF
  TH1D* fHistoLifeTimeV0fromLcAll;    //!<! KF: life time of V0 for all cascades reconstructed with KF

  TH1D* fHistoMassV0fromLcTrue;       //!<! KF: mass of V0 for true cascades reconstructed with KF
  TH1D* fHistoDecayLengthV0fromLcTrue;//!<! KF: decay length of V0 for true cascades reconstructed with KF
  TH1D* fHistoLifeTimeV0fromLcTrue;   //!<! KF: life time of V0 for true cascades reconstructed with KF

  TH1D* fHistoMassLcSgn;              //!<! KF: mass of signal Lc reconstructed with KF
  TH1D* fHistoMassLcSgnFromAOD;       //!<! KF: AOD mass of signal Lc reconstructed with KF
  TH1D* fHistoDecayLengthLcSgn;       //!<! KF: decay length of signal Lc reconstructed with KF
  TH1D* fHistoLifeTimeLcSgn;          //!<! KF: life time of signal Lc reconstructed with KF

  TH1D* fHistoMassV0fromLcSgn;        //!<! KF: mass of V0 for signal Lc reconstructed with KF
  TH1D* fHistoDecayLengthV0fromLcSgn; //!<! KF: decay length of V0 for signal Lc reconstructed with KF
  TH1D* fHistoLifeTimeV0fromLcSgn;    //!<! KF: life time of V0 for signal Lc reconstructed with KF

  TH2D* fHistoKF;                     //!<! KF: V0 code vs Lc code from KF (mass, decaylength, lifetime considered) 
  TH1D* fHistoKFV0;                   //!<! KF: V0 code from KF (mass, decaylength, lifetime considered) 
  TH1D* fHistoKFLc;                   //!<! KF: Lc code from KF (mass, decaylength, lifetime considered) 

  TH2D* fHistoMassKFV0;               //!<! KF: mass vs mass error for V0 from KF  
  TH2D* fHistoDecayLengthKFV0;        //!<! KF: decay length vs decay length error for V0 from KF
  TH2D* fHistoLifeTimeKFV0;           //!<! KF: life time vs life time error for V0 from KF

  TH2D* fHistoMassKFLc;               //!<! KF: mass vs mass error for Lc from KF
  TH2D* fHistoDecayLengthKFLc;        //!<! KF: decay length vs decay length error for Lc from KF
  TH2D* fHistoLifeTimeKFLc;           //!<! KF: life time vs life time error for Lc from KF

  TH2D* fHistoArmenterosPodolanskiV0KF;      //!<! KF: Armeteros-Podolanski plot for all V0 from KF
  TH2D* fHistoArmenterosPodolanskiV0KFSgn;   //!<! KF: Armeteros-Podolanski plot for V0 from signal Lc from KF
  TH2D* fHistoArmenterosPodolanskiV0AOD;     //!<! KF: AOD Armeteros-Podolanski plot for all V0 from KF
  TH2D* fHistoArmenterosPodolanskiV0AODSgn;  //!<! KF: AOD Armeteros-Podolanski plot for V0 from signal Lc from KF

  TList *fOutputKF;                   //!<! User output1: list of histograms from KF

  Int_t fmcLabelLc;                   /// label of candidate
  Bool_t ftopoConstraint;             /// flag to use topological constraints in KF
  Bool_t fCallKFVertexing;            /// flag to decide whether to call or not KF
  Bool_t fKeepingOnlyHIJINGBkg;       /// flag to fill bkg with only candidates that have daughters generated by HIJING (to be used for enriched MC)
  AliVertexingHFUtils* fUtils;        /// AliVertexingHFUtils used to check the generator of a specific candidate
  TH1F* fHistoBackground;             //!<! histo to check the number of candidates with at least one daughter for the injected signal
  Float_t fCutKFChi2NDF;              /// cut for KF on chi2/NDF
  Float_t fCutKFDeviationFromVtx;     /// cut for KF on distance to primary vtx
  Float_t fCutKFDeviationFromVtxV0;   /// cut for KF on distance to primary vtx for V0
  Int_t fCurrentEvent;                /// current event number - for debug purposes
  Double_t fBField;                   /// magnetic field of current event
  Bool_t fKeepingOnlyPYTHIABkg;       /// flag to allow to use only PYTHIA tracks for background
  TH1F* fHistoMCLcK0SpGen;            //!<! histo with MC Lc --> K0S + p
  TH1F* fHistoMCLcK0SpGenAcc;         //!<! histo with MC Lc --> K0S + p
  TH1F* fHistoMCLcK0SpGenLimAcc;      //!<! histo with MC Lc --> K0S + p

  ULong64_t fTriggerMask;	      /// mask to the trigger word returned by the physics selection

  TF1 *fFuncWeightPythia;              //!<! weight function for Pythia vs pPb prod.
  TF1 *fFuncWeightFONLL5overLHC13d3;   //!<! weight function for FONLL vs pPb prod.
  TF1 *fFuncWeightFONLL5overLHC13d3Lc; //!<! weight function for FONLL vs pPb prod.
  TH1F* fHistoMCNch;                   //!<! histogram with Nch distribution from MC production
 
  Int_t fNTracklets_1;                 /// tracklet multiplicity in event in [-1. 1]
  Int_t fNTracklets_All;               /// tracklet multiplicity in event without eta cut
  Float_t fCentrality;                 /// centrality
  
  Bool_t fFillTree;                    /// flag to decide whether to fill the sgn and bkg trees

  Bool_t fUseWeightsLibrary;           // flag to decide whether to use or not the BDT class
  IClassifierReader *fBDTReader;       //!<! BDT reader using BDT class
  TString fTMVAlibName;                /// Name of the library to load to have the TMVA weights
  TString fTMVAlibPtBin;               /// Pt bin that will be in the library to be loaded for the TMVA
  TString fNamesTMVAVar;               /// vector of the names of the input variables
  TH2D *fBDTHisto;                     //!<!
  TH2D *fBDTHistoVsMassK0S;            //!<! BDT classifier vs mass (pi+pi-) pairs
  TH2D *fBDTHistoVstImpParBach;        //!<! BDT classifier vs proton d0
  TH2D *fBDTHistoVstImpParV0;          //!<! BDT classifier vs V0 d0
  TH2D *fBDTHistoVsBachelorPt;         //!<! BDT classifier vs proton pT
  TH2D *fBDTHistoVsCombinedProtonProb; //!<! BDT classifier vs combined proton probability
  TH2D *fBDTHistoVsCtau;               //!<! BDT classifier vs V0 ctau
  TH2D *fBDTHistoVsCosPAK0S;           //!<! BDT classifier vs V0 cosine of pointing angle
  TH2D *fBDTHistoVsSignd0;             //!<! BDT classifier vs V0 proton signed d0
  TH2D *fBDTHistoVsCosThetaStar;       //!<! BDT classifier vs proton emission angle in pK0s pair rest frame
  TH2D* fBDTHistoVsnSigmaTPCpr;       //!<! BDT classifier vs nSigmaTPCpr
  TH2D* fBDTHistoVsnSigmaTOFpr;       //!<! BDT classifier vs nSigmaTOFpr
  TH2D* fBDTHistoVsnSigmaTPCpi;       //!<! BDT classifier vs nSigmaTPCpi
  TH2D* fBDTHistoVsnSigmaTPCka;       //!<! BDT classifier vs nSigmaTPCka
  TH2D* fBDTHistoVsBachelorP;       //!<! BDT classifier vs bachelor p
  TH2D* fBDTHistoVsBachelorTPCP;       //!<! BDT classifier vs bachelor p at TPC wall
  TH2D *fHistoNsigmaTPC;               //!<! 
  TH2D *fHistoNsigmaTOF;               //!<! 

  Bool_t fDebugHistograms;             /// flag to decide whether or not to have extra histograms (useful mainly for debug)

  Int_t fAODProtection;       /// flag to activate protection against AOD-dAOD mismatch.
                                  /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names

  Bool_t fUsePIDresponseForNsigma;  /// flag to decide if to take the nSigma from the PIDresponse or from AliAODPidHF

  Int_t fNVars;  /// Number of training variables

  UInt_t fTimestampCut; // cut on timestamp

  Bool_t fUseXmlWeightsFile;                   // flag to decide whether to use or not the xml file
  TMVA::Reader *fReader;                // TMVA reader using xml file
  Float_t* fVarsTMVA;                   //[fNVars] // variables to be used by TMVA
  Int_t fNVarsSpectators;               // number of spectator variables
  Float_t* fVarsTMVASpectators;         //[fNVarsSpectators] // variables to be used by TMVA
  TString fNamesTMVAVarSpectators;      // vector of the names of the spectators variables
  TString fXmlWeightsFile;              // file with TMVA weights
  TH2D *fBDTHistoTMVA;                  //!<! BDT histo file for the case in which the xml file is used
  Bool_t fUseXmlFileFromCVMFS;          // Boolean to acces Xml from CVMFS path
  TString fXmlFileFromCVMFS;            // Path in CVMFS directory
  
  // Multiplicity corrections
  TProfile* GetEstimatorHistogram(const AliVEvent *event);
  Bool_t fUseMultCorrection;          // flag to decide wether you want to correct the multiplicity
  Double_t fRefMult;                  // refrence multiplcity (period b)
  TProfile* fMultEstimatorAvg[14];    // TProfile with mult vs. Z per period
  Int_t fYearNumber;                  // year number of the data taking
  Int_t fMultiplicityEstimator;       // Definition of the multiplicity estimator: kNtrk10=0, kNtrk10to16=1, kVZERO=2
  Int_t fDoVZER0ParamVertexCorr;      // Flag to use the zvtx correction from (0=none, 1=usual d2h, 2=AliESDUtils for VZERO multiplicity)
  Bool_t fUseMultiplicityCut;
  Float_t fMultiplicityCutMin;        // Minimum multiplicity cut, minimum is included
  Float_t fMultiplicityCutMax;        // Maximum multiplicity cut, maximum is excluded

  TH1F* fHistoNtrUnCorr;             //!<! hist. with number of uncorrected tracklets
  TH1F* fHistoNtrCorr;               //!<! hist. with number of corrected tracklets
  TH2F* fHistoVzVsNtrUnCorr;         //!<! hist. Vz vs UNCORRECTED tracklets
  TH2F* fHistoVzVsNtrCorr;           //!<! hist. Vz vs corrected tracklets
  
  /// \cond CLASSIMP    
  ClassDef(AliAnalysisTaskSELc2V0bachelorTMVAApp, 12); /// class for Lc->p K0
  /// \endcond    
};

#endif

