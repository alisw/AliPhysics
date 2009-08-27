/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/********************************** 
 * flow analysis with Q-cumulants * 
 *                                * 
 * author:  Ante Bilandzic        * 
 *           (anteb@nikhef.nl)    *
 *********************************/ 

#ifndef ALIFLOWANALYSISWITHQCUMULANTS_H
#define ALIFLOWANALYSISWITHQCUMULANTS_H

#include "AliFlowCommonConstants.h" // needed as include
#include "TMatrixD.h"
#include "TH2D.h"
#include "TBits.h"

class TObjArray;
class TList;
class TFile;
class TGraph;

class TH1;
class TProfile;
class TProfile2D;

class AliFlowEventSimple;
class AliFlowVector;

class AliFlowCommonHist;
class AliFlowCommonHistResults;

//================================================================================================================

class AliFlowAnalysisWithQCumulants{
 public:
  AliFlowAnalysisWithQCumulants();
  virtual ~AliFlowAnalysisWithQCumulants(); 
  // 0.) methods called in the constructor:
  virtual void InitializeArraysForIntFlow();
  virtual void InitializeArraysForDiffFlow();
  virtual void InitializeArraysForDistributions();
  // 1.) method Init() and methods called within Init():
  virtual void Init();
    virtual void AccessConstants();
    virtual void BookAndNestAllLists();
    virtual void BookCommonHistograms();
    virtual void BookAndFillWeightsHistograms();
    virtual void BookEverythingForIntegratedFlow();
    virtual void BookEverythingForDifferentialFlow();
    virtual void BookEverythingForDistributions();  
    virtual void BookEverythingForNestedLoops();   
  // 2.) method Make() and methods called within Make():
  virtual void Make(AliFlowEventSimple *anEvent);
    // 2a.) integrated flow:
    virtual void FillAverageMultiplicities(Int_t nRP);
    virtual void CalculateCorrelationsForIntegratedFlow(); 
    virtual void CalculateCorrectionsForNonUniformAcceptanceForIntFlowCosTerms();  
    virtual void CalculateCorrectionsForNonUniformAcceptanceForIntFlowSinTerms();  
    virtual void CalculateQProductsForIntFlow();
    virtual void CalculateSumAndProductOfEventWeights();
    virtual void CalculateWeightedCorrelationsForIntegratedFlow();
    virtual void CalculateWeightedQProductsForIntFlow();
    virtual void EvaluateNestedLoopsForIntegratedFlow(AliFlowEventSimple* anEvent); 
    // 2b.) differential flow:
    virtual void CalculateCorrelationsForDifferentialFlow(TString type);
    virtual void CalculateCorrectionsForNonUniformAcceptanceForDifferentialFlowCosTerms(TString type);  
    virtual void CalculateCorrectionsForNonUniformAcceptanceForDifferentialFlowSinTerms(TString type);  
    virtual void CalculateWeightedCorrelationsForDifferentialFlow(TString type); 
    virtual void EvaluateNestedLoopsForDifferentialFlow(AliFlowEventSimple* anEvent);
  // 3.) method Finish() and methods called within Finish():
  virtual void Finish();
    // 3a.) integrated flow:
    virtual void FinalizeCorrelationsForIntFlow(Bool_t useParticleWeights, TString eventWeights);
    virtual void CalculateFinalCorrectionsForNonUniformAcceptanceForCumulantsForIntFlow(Bool_t useParticleWeights, TString eventWeights);
    virtual void CalculateCovariancesForIntFlow(Bool_t useParticleWeights, TString eventWeights);  
    virtual void CalculateCumulantsForIntFlow(Bool_t useParticleWeights, TString eventWeights); 
    virtual void ApplyCorrectionForNonUniformAcceptanceToCumulantsForIntFlow(Bool_t useParticleWeights, TString eventWeights); 
    virtual void CalculateIntFlow(Bool_t useParticleWeights, TString eventWeights, Bool_t correctedForNUA); 
    virtual void FillCommonHistResultsIntFlow(Bool_t useParticleWeights, TString eventWeights, Bool_t correctedForNUA);
    virtual void PrintQuantifyingCorrectionsForNonUniformAcceptance(Bool_t useParticleWeights, TString eventWeights);
    virtual void PrintFinalResultsForIntegratedFlow(TString type);
    virtual void CompareResultsFromNestedLoopsAndFromQVectorsForIntFlow(Bool_t useParticleWeights);
    // 3b.) differential flow:
    virtual void FinalizeCorrelationsForDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights);
    virtual void CalculateFinalCorrectionsForNonUniformAcceptanceForDifferentialFlow(Bool_t useParticleWeights, TString type);
    virtual void CalculateCumulantsForDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights); 
    virtual void CalculateDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights); 
    virtual void FillCommonHistResultsDiffFlow(TString type, Bool_t useParticleWeights, TString eventWeights, Bool_t correctedForNUA);
    virtual void CalculateFinalResultsForRPandPOIIntegratedFlow(TString type, Bool_t useParticleWeights, TString eventWeights);
    virtual void CompareResultsFromNestedLoopsAndFromQVectorsForDiffFlow(Bool_t useParticleWeights);  
  // 4.) other methods: 
  TProfile* MakePtProjection(TProfile2D *profilePtEta) const;
  TProfile* MakeEtaProjection(TProfile2D *profilePtEta) const;
  virtual void GetOutputHistograms(TList *outputListHistos); 
  virtual void WriteHistograms(TString outputFileName);
  
  // **** SETTERS and GETTERS ****
  
  // 0.) base:                                                                                              
  TList* GetHistList() const {return this->fHistList;} 
  
  // 1.) common:
  void SetCommonHists(AliFlowCommonHist* const ch) {this->fCommonHists = ch;};
  AliFlowCommonHist* GetCommonHists() const {return this->fCommonHists;};
  void SetCommonHists2nd(AliFlowCommonHist* const ch2nd) {this->fCommonHists2nd = ch2nd;};
  AliFlowCommonHist* GetCommonHists2nd() const {return this->fCommonHists2nd;};
  void SetCommonHists4th(AliFlowCommonHist* const ch4th) {this->fCommonHists4th = ch4th;};
  AliFlowCommonHist* GetCommonHists4th() const {return this->fCommonHists4th;};
  void SetCommonHists6th(AliFlowCommonHist* const ch6th) {this->fCommonHists6th = ch6th;};
  AliFlowCommonHist* GetCommonHists6th() const {return this->fCommonHists6th;};
  void SetCommonHists8th(AliFlowCommonHist* const ch8th) {this->fCommonHists8th = ch8th;};
  AliFlowCommonHist* GetCommonHists8th() const {return this->fCommonHists8th;};
  void SetCommonHistsResults2nd(AliFlowCommonHistResults* const chr2nd) {this->fCommonHistsResults2nd = chr2nd;};
  AliFlowCommonHistResults* GetCommonHistsResults2nd() const {return this->fCommonHistsResults2nd;};
  void SetCommonHistsResults4th(AliFlowCommonHistResults* const chr4th) {this->fCommonHistsResults4th = chr4th;};
  AliFlowCommonHistResults* GetCommonHistsResults4th() const {return this->fCommonHistsResults4th;};
  void SetCommonHistsResults6th(AliFlowCommonHistResults* const chr6th) {this->fCommonHistsResults6th = chr6th;};
  AliFlowCommonHistResults* GetCommonHistsResults6th() const {return this->fCommonHistsResults6th;};
  void SetCommonHistsResults8th(AliFlowCommonHistResults* const chr8th) {this->fCommonHistsResults8th = chr8th;};
  AliFlowCommonHistResults* GetCommonHistsResults8th() const {return this->fCommonHistsResults8th;};
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  Int_t GetHarmonic() const {return this->fHarmonic;};
  void SetAnalysisLabel(const char *aLabel) {this->fAnalysisLabel->Append(*aLabel);};
  TString *GetAnalysisLabel() const {return this->fAnalysisLabel;};
  
  // 2.) weights:
  void SetWeightsList(TList* wlist) {this->fWeightsList = (TList*)wlist->Clone();}
  TList* GetWeightsList() const {return this->fWeightsList;}  
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
  void SetUseParticleWeights(TProfile* const uPW) {this->fUseParticleWeights = uPW;};
  TProfile* GetUseParticleWeights() const {return this->fUseParticleWeights;};
  void SetPhiWeights(TH1F* const histPhiWeights) {this->fPhiWeights = histPhiWeights;};
  TH1F* GetPhiWeights() const {return this->fPhiWeights;};
  void SetPtWeights(TH1D* const histPtWeights) {this->fPtWeights = histPtWeights;};
  TH1D* GetPtWeights() const {return this->fPtWeights;};
  void SetEtaWeights(TH1D* const histEtaWeights) {this->fEtaWeights = histEtaWeights;};
  TH1D* GetEtaWeights() const {return this->fEtaWeights;};
  
  // 3.) integrated flow:
  // integrated flow profiles:
  void SetAvMultiplicity(TProfile* const avMultiplicity) {this->fAvMultiplicity = avMultiplicity;};
  TProfile* GetAvMultiplicity() const {return this->fAvMultiplicity;};
  void SetQCorrelations(TProfile* const qCorrel, Int_t pW, Int_t eW) {this->fQCorrelations[pW][eW] = qCorrel;};
  TProfile* GetQCorrelations(Int_t pW, Int_t eW) const {return this->fQCorrelations[pW][eW];}; 
  void SetQProducts(TProfile* const qProduct, Int_t pW, Int_t eW) {this->fQProducts[pW][eW] = qProduct;};
  TProfile* GetQProducts(Int_t pW, Int_t eW) const {return this->fQProducts[pW][eW];};
  void SetQCorrections(TProfile* const qCorrections, Int_t pW, Int_t eW, Int_t sc) {this->fQCorrections[pW][eW][sc] = qCorrections;};
  TProfile* GetQCorrections(Int_t pW, Int_t eW, Int_t sc) const {return this->fQCorrections[pW][eW][sc];};
  // integrated flow results:
  void SetCorrelations(TH1D* const correl, Int_t pW, Int_t eW) {this->fCorrelations[pW][eW] = correl;};
  TH1D* GetCorrelations(Int_t pW, Int_t eW) const {return this->fCorrelations[pW][eW];};
  void SetCorrections(TH1D* const correct, Int_t pW, Int_t eW) {this->fCorrections[pW][eW] = correct;};
  TH1D* GetCorrections(Int_t pW, Int_t eW) const {return this->fCorrections[pW][eW];};
  void SetCovariances(TH1D* const cov, Int_t pW, Int_t eW) {this->fCovariances[pW][eW] = cov;};
  TH1D* GetCovariances(Int_t pW, Int_t eW) const {return this->fCovariances[pW][eW];};
  void SetSumOfEventWeights(TH1D* const soew, Int_t pW, Int_t eW, Int_t power) {this->fSumOfEventWeights[pW][eW][power] = soew;};
  TH1D* GetSumOfEventWeights(Int_t pW, Int_t eW, Int_t power) const {return this->fSumOfEventWeights[pW][eW][power];};
  void SetProductOfEventWeights(TH1D* const poew, Int_t pW, Int_t eW) {this->fProductOfEventWeights[pW][eW] = poew;};
  TH1D* GetProductOfEventWeights(Int_t pW, Int_t eW) const {return this->fProductOfEventWeights[pW][eW];};
  void SetCumulants(TH1D* const cumulants, Int_t pW, Int_t eW, Int_t nua) {this->fCumulants[pW][eW][nua] = cumulants;};
  TH1D* GetCumulants(Int_t pW, Int_t eW, Int_t nua) const {return this->fCumulants[pW][eW][nua];};
  void SetIntFlow(TH1D* const intFlow, Int_t pW, Int_t eW, Int_t nua) {this->fIntFlow[pW][eW][nua] = intFlow;};
  TH1D* GetIntFlow(Int_t pW, Int_t eW, Int_t nua) const {return this->fIntFlow[pW][eW][nua];};
  
  // 4.) differential flow:
  // profiles:
  void SetCorrelationsPro(TProfile2D* const correlPro, Int_t i, Int_t j, Int_t k, Int_t l) {this->fCorrelationsPro[i][j][k][l] = correlPro;};
  TProfile2D* GetCorrelationsPro(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fCorrelationsPro[i][j][k][l];};
  void SetProductsOfCorrelationsPro(TProfile2D* const proOfcorrelPro, Int_t i, Int_t j, Int_t k, Int_t l) {this->fProductsOfCorrelationsPro[i][j][k][l] = proOfcorrelPro;};
  TProfile2D* GetProductsOfCorrelationsPro(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fProductsOfCorrelationsPro[i][j][k][l];};
  void SetCorrectionTermsPro(TProfile2D* const correctTermsPro, Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) {this->fCorrectionTermsPro[i][j][k][l][m] = correctTermsPro;};
  TProfile2D* GetCorrectionTermsPro(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) const {return this->fCorrectionTermsPro[i][j][k][l][m];};  
  // results:
  void SetFinalCorrelations2D(TH2D* const fCorrelations2D, Int_t i, Int_t j, Int_t k, Int_t l) {this->fFinalCorrelations2D[i][j][k][l] = fCorrelations2D;};
  TH2D* GetFinalCorrelations2D(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fFinalCorrelations2D[i][j][k][l];};
  void SetFinalCorrelations1D(TH1D* const fCorrelations1D, Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) {this->fFinalCorrelations1D[i][j][k][l][m] = fCorrelations1D;};
  TH1D* GetFinalCorrelations1D(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) const {return this->fFinalCorrelations1D[i][j][k][l][m];};
  void SetFinalCorrections2D(TH2D* const fCorrections2D, Int_t i, Int_t j, Int_t k, Int_t l) {this->fFinalCorrections2D[i][j][k][l] = fCorrections2D;};
  TH2D* GetFinalCorrections2D(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fFinalCorrections2D[i][j][k][l];};
  void SetFinalCorrections1D(TH1D* const fCorrections1D, Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) {this->fFinalCorrections1D[i][j][k][l][m] = fCorrections1D;};
  TH1D* GetFinalCorrections1D(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) const {return this->fFinalCorrections1D[i][j][k][l][m];}; 
  void SetFinalCovariances2D(TH2D* const fCovariances2D, Int_t i, Int_t j, Int_t k, Int_t l) {this->fFinalCovariances2D[i][j][k][l] = fCovariances2D;};
  TH2D* GetFinalCovariances2D(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fFinalCovariances2D[i][j][k][l];};
  void SetFinalCovariances1D(TH1D* const fCovariances1D, Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) {this->fFinalCovariances1D[i][j][k][l][m] = fCovariances1D;};
  TH1D* GetFinalCovariances1D(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) const {return this->fFinalCovariances1D[i][j][k][l][m];};  
  void SetFinalCumulants2D(TH2D* const fCumulants2D, Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) {this->fFinalCumulants2D[i][j][k][l][m] = fCumulants2D;};
  TH2D* GetFinalCumulants2D(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) const {return this->fFinalCumulants2D[i][j][k][l][m];};
  void SetFinalCumulantsPt(TH1D* const fCumulantsPt, Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) {this->fFinalCumulantsPt[i][j][k][l][m] = fCumulantsPt;};
  TH1D* GetFinalCumulantsPt(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) const {return this->fFinalCumulantsPt[i][j][k][l][m];};
  void SetFinalCumulantsEta(TH1D* const fCumulantsEta, Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) {this->fFinalCumulantsEta[i][j][k][l][m] = fCumulantsEta;};
  TH1D* GetFinalCumulantsEta(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) const {return this->fFinalCumulantsEta[i][j][k][l][m];};
  void SetFinalFlow2D(TH2D* const fFlow2D, Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) {this->fFinalFlow2D[i][j][k][l][m] = fFlow2D;};
  TH2D* GetFinalFlow2D(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) const {return this->fFinalFlow2D[i][j][k][l][m];};
  void SetFinalFlowPt(TH1D* const fFlowPt, Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) {this->fFinalFlowPt[i][j][k][l][m] = fFlowPt;};
  TH1D* GetFinalFlowPt(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) const {return this->fFinalFlowPt[i][j][k][l][m];};
  void SetFinalFlowEta(TH1D* const fFlowEta, Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) {this->fFinalFlowEta[i][j][k][l][m] = fFlowEta;};
  TH1D* GetFinalFlowEta(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) const {return this->fFinalFlowEta[i][j][k][l][m];};
  void SetNonEmptyBins2D(TH2D* const fneb2D, Int_t i) {this->fNonEmptyBins2D[i] = fneb2D;};
  TH2D* GetNonEmptyBins2D(Int_t i) const {return this->fNonEmptyBins2D[i];};
  void SetNonEmptyBins1D(TH1D* const fneb1D, Int_t i, Int_t j) {this->fNonEmptyBins1D[i][j] = fneb1D;};
  TH1D* GetNonEmptyBins1D(Int_t i, Int_t j) const {return this->fNonEmptyBins1D[i][j];};
  
   
  // x.) debugging and cross-checking:
  void SetNestedLoopsList(TList* nllist) {this->fNestedLoopsList = nllist;};
  TList* GetNestedLoopsList() const {return this->fNestedLoopsList;}; 
  void SetEvaluateNestedLoopsForIntFlow(Bool_t const enlfif) {this->fEvaluateNestedLoopsForIntFlow = enlfif;};
  Bool_t GetEvaluateNestedLoopsForIntFlow() const {return this->fEvaluateNestedLoopsForIntFlow;};
  void SetEvaluateNestedLoopsForDiffFlow(Bool_t const enlfdf) {this->fEvaluateNestedLoopsForDiffFlow = enlfdf;};
  Bool_t GetEvaluateNestedLoopsForDiffFlow() const {return this->fEvaluateNestedLoopsForDiffFlow;};
  void SetEvaluateNestedLoops(TProfile* const enl) {this->fEvaluateNestedLoops = enl;};
  TProfile* GetEvaluateNestedLoops() const {return this->fEvaluateNestedLoops;};
  void SetDirectCorrelations(TProfile* const dc) {this->fDirectCorrelations = dc;};
  TProfile* GetDirectCorrelations() const {return this->fDirectCorrelations;};
  void SetDirectCorrectionsCos(TProfile* const dcc) {this->fDirectCorrectionsCos = dcc;};
  TProfile* GetDirectCorrectionsCos() const {return this->fDirectCorrectionsCos;};
  void SetDirectCorrectionsSin(TProfile* const dcs) {this->fDirectCorrectionsSin = dcs;};
  TProfile* GetDirectCorrectionsSin() const {return this->fDirectCorrectionsSin;};
  void SetDirectCorrelationsDiffFlow(TProfile* const dcdf) {this->fDirectCorrelationsDiffFlow = dcdf;};
  TProfile* GetDirectCorrelationsDiffFlow() const {return this->fDirectCorrelationsDiffFlow;};
  void SetDirectCorrectionsDiffFlowCos(TProfile* const dcdfc) {this->fDirectCorrectionsDiffFlowCos = dcdfc;};
  TProfile* GetDirectCorrectionsDiffFlowCos() const {return this->fDirectCorrectionsDiffFlowCos;};
  void SetDirectCorrectionsDiffFlowSin(TProfile* const dcdfs) {this->fDirectCorrectionsDiffFlowSin = dcdfs;};
  TProfile* GetDirectCorrectionsDiffFlowSin() const {return this->fDirectCorrectionsDiffFlowSin;};
  void SetDirectCorrelationsW(TProfile* const dcw) {this->fDirectCorrelationsW = dcw;};
  TProfile* GetDirectCorrelationsW() const {return this->fDirectCorrelationsW;};
  void SetDirectCorrectionsCosW(TProfile* const dccw) {this->fDirectCorrectionsCosW = dccw;};
  TProfile* GetDirectCorrectionsCosW() const {return this->fDirectCorrectionsCosW;};
  void SetDirectCorrectionsSinW(TProfile* const dcsw) {this->fDirectCorrectionsSinW = dcsw;};
  TProfile* GetDirectCorrectionsSinW() const {return this->fDirectCorrectionsSinW;};
  void SetDirectCorrelationsDiffFlowW(TProfile* const dcdfw) {this->fDirectCorrelationsDiffFlowW = dcdfw;};
  TProfile* GetDirectCorrelationsDiffFlowW() const {return this->fDirectCorrelationsDiffFlowW;};
  void SetDirectCorrectionsDiffFlowCosW(TProfile* const dcdfcw) {this->fDirectCorrectionsDiffFlowCosW = dcdfcw;};
  TProfile* GetDirectCorrectionsDiffFlowCosW() const {return this->fDirectCorrectionsDiffFlowCosW;};
  void SetDirectCorrectionsDiffFlowSinW(TProfile* const dcdfsw) {this->fDirectCorrectionsDiffFlowSinW = dcdfsw;};
  TProfile* GetDirectCorrectionsDiffFlowSinW() const {return this->fDirectCorrectionsDiffFlowSinW;};
  
 private:
  
  AliFlowAnalysisWithQCumulants(const AliFlowAnalysisWithQCumulants& afawQc);
  AliFlowAnalysisWithQCumulants& operator=(const AliFlowAnalysisWithQCumulants& afawQc); 
  
  // 0.) base:
  TList* fHistList; // base list to hold all output object
  
  // 1.) common:
  AliFlowCommonHist *fCommonHists; // common control histograms (taking into account ALL events) 
  AliFlowCommonHist *fCommonHists2nd; // common control histograms (taking into account only the events with 2 and more particles) 
  AliFlowCommonHist *fCommonHists4th; // common control histograms (taking into account only the events with 4 and more particles) 
  AliFlowCommonHist *fCommonHists6th; // common control histograms (taking into account only the events with 6 and more particles) 
  AliFlowCommonHist *fCommonHists8th; // common control histograms (taking into account only the events with 8 and more particles) 
  AliFlowCommonHistResults *fCommonHistsResults2nd; // final results for 2nd order int. and diff. flow for events with 2 and more particles
  AliFlowCommonHistResults *fCommonHistsResults4th; // final results for 4th order int. and diff. flow for events with 4 and more particles 
  AliFlowCommonHistResults *fCommonHistsResults6th; // final results for 6th order int. and diff. flow for events with 6 and more particles
  AliFlowCommonHistResults *fCommonHistsResults8th; // final results for 8th order int. and diff. flow for events with 8 and more particles
  Int_t fnBinsPhi; // number of phi bins
  Double_t fPhiMin; // minimum phi   
  Double_t fPhiMax; // maximum phi 
  Double_t fPhiBinWidth; // bin width for phi histograms  
  Int_t fnBinsPt; // number of pt bins
  Double_t fPtMin; // minimum pt   
  Double_t fPtMax; // maximum pt  
  Double_t fPtBinWidth; // bin width for pt histograms  
  Int_t fnBinsEta; // number of eta bins
  Double_t fEtaMin; // minimum eta   
  Double_t fEtaMax; // maximum eta
  Double_t fEtaBinWidth; // bin width for eta histograms  
  Int_t fHarmonic; // harmonic 
  TString *fAnalysisLabel; // analysis label (all histograms and output file will have this label)
  
  // 2.) weights
  TList *fWeightsList; // list to hold all histograms with particle weights: fUseParticleWeights, fPhiWeights, fPtWeights and fEtaWeights
  Bool_t fUsePhiWeights; // use phi weights
  Bool_t fUsePtWeights; // use pt weights
  Bool_t fUseEtaWeights; // use eta weights
  TProfile *fUseParticleWeights; // profile with three bins to hold values of fUsePhiWeights, fUsePtWeights and fUseEtaWeights
  TH1F *fPhiWeights; // histogram holding phi weights
  TH1D *fPtWeights; // histogram holding phi weights
  TH1D *fEtaWeights; // histogram holding phi weights 
  
  // 3.) integrated flow       
  TList *fIntFlowList; // list to hold all histograms relevant for integrated flow 
  TList *fIntFlowProfiles; // list to hold all profiles relevant for integrated flow
  TList *fIntFlowResults; // list to hold all histograms with final results relevant for integrated flow  
  //  3a.) event-by-event quantities:
  TMatrixD *fReQ; // fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
  TMatrixD *fImQ; // fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
  TMatrixD *fSMpk; // fSM[p][k] = (sum_{i=1}^{M} w_{i}^{k})^{p}
  TH1D *fQCorrelationsEBE[2]; // [0=weights not used,1=weights used]
  TH1D *fQCorrectionsEBE[2][2]; // [0=weights not used,1=weights used][0=sin terms, 1=cos terms]
  //  3b.) profiles:
  TProfile *fAvMultiplicity; // profile to hold average multiplicities and number of events for events with nRP>=0, nRP>=1, ... , and nRP>=8
  TProfile *fQCorrelations[2][2]; // [0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights] 
  TProfile *fQProducts[2][2]; // [0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights] 
  TProfile *fQCorrections[2][2][2]; // [0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights][0=sin terms, 1=cos terms] 
  //  3c.) results:
  TH1D *fCorrelations[2][2]; // final results for average correlations: [0=pW not used,1=pW used][0=exact eW,1=non-exact eW]
  TH1D *fCorrections[2][2]; // corrections for non-uniform acceptance to integrated Q-cumulants: [0=pW not used,1=pW used][0=exact eW,1=non-exact eW]
  TH1D *fCovariances[2][2]; // covariances of multi-particle correlations: [0=pW not used,1=pW used][0=exact eW,1=non-exact eW]
  TH1D *fSumOfEventWeights[2][2][2]; // [0=pW not used,1=pW used][0=exact eW,1=non-exact eW][0=power 1,1=power 2][0=eW for <2>, 1=eW for <4>, ...]
  TH1D *fProductOfEventWeights[2][2]; // [0=pW not used,1=pW used][0=exact eW,1=non-exact eW][0=eWs for <2><4>, 1=eWs for <2><6>, ...]
  TH1D *fCumulants[2][2][2]; // integrated Q-cumulants: [0=pW not used,1=pW used][0=exact eW,1=non-exact eW][0=not corrected, 1=corrected]
  TH1D *fIntFlow[2][2][2]; // int. flow estimates from Q-cumulants: [0=pW not used,1=pW used][0=exact eW,1=non-exact eW][0=not corrected, 1=corrected]
  
  
  // 4.) differential flow
  TList *fDiffFlowList;
  // nested lists to hold profiles:
  TList *fDiffFlowProfiles; // list to hold all lists with profiles with correlations, correction terms for NUA and products of correlations
  TList *fDFPType[2]; // [0=RP,1=POI]
  TList *fDFPParticleWeights[2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used]
  TList *fDFPEventWeights[2][2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights]
  TList *fDiffFlowCorrelations[2][2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights]
  TList *fDiffFlowProductsOfCorrelations[2][2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights]
  TList *fDiffFlowCorrectionTerms[2][2][2][2]; // [0=RP,1=POI][0=pW not used,1=pW used][0=exact eWeights,1=non-exact eWeights][0=sin terms,1=cos terms]
  //  4a.) event-by-event quantities:
  TProfile2D *fReEBE[3][4][9]; // real part of r_{m*n,k}(pt,eta), p_{m*n,k}(pt,eta) and q_{m*n,k}(pt,eta)
  TProfile2D *fImEBE[3][4][9]; // imaginary part of r_{m*n,k}(pt,eta), p_{m*n,k}(pt,eta) and q_{m*n,k}(pt,eta)
  TProfile2D *fs[3][9]; // [t][k] // to be improved
  //  4b.) profiles:
  TProfile2D *fCorrelationsPro[2][2][2][4]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights][corr.'s index]
  TProfile2D *fProductsOfCorrelationsPro[2][2][2][5]; // [0=RP,1=POI][0=pW not used,1=pW used][0=exact eWeights,1=non-exact eWeights][products' index]
  TProfile2D *fCorrectionTermsPro[2][2][2][2][2]; // [0=RP,1=POI][0=pW not used,1=pW used][0=e eW,1=ne eW][0=sin terms,1=cos terms][corr. terms' index]
  // nested lists to hold histograms with results:
  TList *fDiffFlowResults; // list to hold all lists with histos with results for diff. flow, cumulants, correlations, covariances and corrections 
  TList *fDFRType[2]; // [0=RP,1=POI] 
  TList *fDFRParticleWeights[2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used]
  TList *fDFREventWeights[2][2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights]  
  TList *fDFRCorrections[2][2][2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights][0=not corr,1=corr]  
  TList *fDiffFlowFinalCorrelations[2][2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights]
  TList *fDiffFlowFinalCorrections[2][2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights] 
  TList *fDiffFlowFinalCovariances[2][2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights] 
  TList *fDiffFlowFinalCumulants[2][2][2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eW,1=non-exact eW][0=not corr,1=corr] 
  TList *fDiffFlowFinalFlow[2][2][2][2]; // [0=RP,1=POI][0=pWeights not used,1=pWeights used][0=exact eW,1=non-exact eW][0=not corr,1=corr]
  TH2D *fFinalCorrelations2D[2][2][2][4]; // [0=RP,1=POI][0=pW not used,1=pW used][0=exact eW,1=non-exact eW][index of correlations] 
  TH1D *fFinalCorrelations1D[2][2][2][2][4]; // [0=RP,1=POI][0=pW not used,1=pW used][0=exact eW,1=non-exact eW][0=pt,1=eta][index of correlations] 
  TH2D *fFinalCorrections2D[2][2][2][4]; // [0=RP,1=POI][0=pW not used,1=pW used][0=exact eW,1=non-exact eW][corr.to cumulant of order] 
  TH1D *fFinalCorrections1D[2][2][2][2][4]; // [0=RP,1=POI][0=pW not used,1=pW used][0=exact eW,1=non-exact eW][0=pt,1=eta][corr. to cumulant of order] 
  TH2D *fFinalCovariances2D[2][2][2][4]; // [0=RP,1=POI][0=pW not used,1=pW used][0=exact eW,1=non-exact eW][index of covariances] 
  TH1D *fFinalCovariances1D[2][2][2][2][4]; // [0=RP,1=POI][0=pW not used,1=pW used][0=exact eW,1=non-exact eW][0=pt,1=eta][index of covariances] 
  TH2D *fFinalCumulants2D[2][2][2][2][4]; // [0=RP,1=POI][0=pW not used,1=pW used][0=e eW,1=ne eW][0=not corr,1=corr][cumulant's order] 
  TH1D *fFinalCumulantsPt[2][2][2][2][4]; // [0=RP,1=POI][0=pW nu,1=pW u][0=e eW,1=ne eW][0=nc,1=corr][cumulant's order] 
  TH1D *fFinalCumulantsEta[2][2][2][2][4]; // [0=RP,1=POI][0=pW nu,1=pW u][0=e eW,1=ne eW][0=nc,1=corr][cumulant's order]
  TH2D *fFinalFlow2D[2][2][2][2][4]; // [0=RP,1=POI][0=pW not used,1=pW used][0=e eW,1=ne eW][0=not corr,1=corr][order of flow estimate] 
  TH1D *fFinalFlowPt[2][2][2][2][4]; // [0=RP,1=POI][0=pW not used,1=pW used][0=e eW,1=ne eW][0=not corr,1=corr][order of flow estimate] 
  TH1D *fFinalFlowEta[2][2][2][2][4]; // [0=RP,1=POI][0=pW not used,1=pW used][0=e eW,1=ne eW][0=not corr,1=corr][order of flow estimate] 
  TH2D *fNonEmptyBins2D[2]; // [0=RP,1=POI]
  TH1D *fNonEmptyBins1D[2][2]; // [0=RP,1=POI][0=pt,1=eta]
      
  // 5.) distributions:
  TList *fDistributionsList; // list to hold all distributions
  TH1D *fDistributions[2][2][4]; // [0=pWeights not used,1=pWeights used][0=exact eWeights,1=non-exact eWeights][0=<2>,1=<4>,2=<6>,3=<8>]
    
  // x.) debugging and cross-checking:
  // Remark: for weighted correlations cross-checking is performed only with phi weights (this is sufficient) 
  TList *fNestedLoopsList; // list to hold all profiles filled with nested loops
  Bool_t fEvaluateNestedLoopsForIntFlow; // evaluate nested loops relevant for integrated flow
  Bool_t fEvaluateNestedLoopsForDiffFlow; // evaluate nested loops relevant for differential flow
  TProfile *fEvaluateNestedLoops; // profile with two bins to hold values of fEvaluateNestedLoopsForIntFlow and fEvaluateNestedLoopsForDiffFlow
  TProfile *fDirectCorrelations; // reduced multi-particle correlations calculated with nested loop relevant for int. flow 
  TProfile *fDirectCorrectionsCos; // corrections for non-uniform acceptance (cos terms) calculated with nested loops (int. flow)
  TProfile *fDirectCorrectionsSin; // corrections for non-uniform acceptance (sin terms) calculated with nested loops (int. flow)
  TProfile *fDirectCorrelationsDiffFlow; // multi-particle correlations calculated with nested loop relevant for diff. flow
  TProfile *fDirectCorrectionsDiffFlowCos; // corrections for non-uniform acceptance (cos terms) calculated with nested loops (diff. flow)
  TProfile *fDirectCorrectionsDiffFlowSin; // corrections for non-uniform acceptance (sin terms) calculated with nested loops (diff. flow)
  TProfile *fDirectCorrelationsW; // weighted multi-particle correlations calculated with nested loop relevant for int. flow 
  TProfile *fDirectCorrectionsCosW; // weighted corrections for non-uniform acceptance (cos terms) calculated with nested loops (int. flow)
  TProfile *fDirectCorrectionsSinW; // weighted corrections for non-uniform acceptance (sin terms) calculated with nested loops (int. flow)
  TProfile *fDirectCorrelationsDiffFlowW; // weighted reduced multi-particle correlations calculated with nested loop relevant for diff. flow
  TProfile *fDirectCorrectionsDiffFlowCosW; // weighted corrections for non-uniform acceptance (cos terms) calculated with nested loops (diff. flow)
  TProfile *fDirectCorrectionsDiffFlowSinW; // weighted corrections for non-uniform acceptance (sin terms) calculated with nested loops (diff. flow)
                  
  ClassDef(AliFlowAnalysisWithQCumulants, 0);
};

//================================================================================================================

#endif





