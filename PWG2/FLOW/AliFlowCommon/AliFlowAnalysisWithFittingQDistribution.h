/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/******************************** 
 * estimating reference flow by *
 *   fitting q-distribution     * 
 *                              *
 * author: Ante Bilandzic       * 
 *       (abilandzic@gmail.com) *
 *                              *  
 *  based on the macro written  *
 *     by Sergei Voloshin       *
 *******************************/ 
 
#ifndef ALIFLOWANALYSISWITHFITTINGQDISTRIBUTION_H
#define ALIFLOWANALYSISWITHFITTINGQDISTRIBUTION_H

#include "AliFlowCommonConstants.h"

class TObjArray;
class TList;
class TFile;
class TDirectoryFile;

class TH1F;
class TH1D;
class TProfile;
class TF1;

class AliFlowEventSimple;
class AliFlowTrackSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;
class AliFlowVector;

//================================================================================================================

class AliFlowAnalysisWithFittingQDistribution{
 public:
  AliFlowAnalysisWithFittingQDistribution();
  virtual ~AliFlowAnalysisWithFittingQDistribution(); 
  // 0.) methods called in the constructor:
  virtual void InitializeArrays();
  // 1.) method Init() and methods called within Init():
  virtual void Init();
   virtual void AccessConstants();
   virtual void BookCommonHistograms();
   virtual void BookAndFillWeightsHistograms();
   virtual void BookEverythingForDistributions();
   virtual void StoreFittingParameters();
   virtual void AccessFittingParameters();
  // 2.) method Make() and methods called within Make(): 
  virtual void Make(AliFlowEventSimple* anEvent);
   virtual void CheckPointersUsedInMake();
  // 3.) method Finish() and methods called within Finish(): 
  virtual void Finish(Bool_t doFit = kTRUE);
   virtual void CheckPointersUsedInFinish();
   virtual void DoFit(Bool_t sigma2Fitted);
   virtual void FillCommonHistResults(Bool_t sigma2Fitted);
   virtual void PrintOnTheScreen();
  // 4.) other methods:
  virtual void GetOutputHistograms(TList *outputListHistos); 
  virtual void WriteHistograms(TString *outputFileName);
  virtual void WriteHistograms(TString outputFileName);
  virtual void WriteHistograms(TDirectoryFile *outputFileName);
    
  // **** SETTERS and GETTERS ****
  
  // 0.) base:                                                                                              
  TList* GetHistList() const {return this->fHistList;} 
  // 1.) common:
  void SetCommonHists(AliFlowCommonHist* const ch) {this->fCommonHists = ch;};
  AliFlowCommonHist* GetCommonHists() const {return this->fCommonHists;};
  void SetCommonHistsResults(AliFlowCommonHistResults* const chr) {this->fCommonHistsResults = chr;};
  AliFlowCommonHistResults* GetCommonHistsResults() const {return this->fCommonHistsResults;};
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  Int_t GetHarmonic() const {return this->fHarmonic;};
  void SetAnalysisLabel(const char *aLabel) {this->fAnalysisLabel->Append(*aLabel);};
  TString *GetAnalysisLabel() const {return this->fAnalysisLabel;};
  // 2.) weights:
  void SetWeightsList(TList* wlist) {this->fWeightsList = (TList*)wlist->Clone();};
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
  // 3.) distributions:
  void SetSumOfParticleWeights(TH1D* const sopW) {this->fSumOfParticleWeights = sopW;};
  TH1D* GetSumOfParticleWeights() const {return this->fSumOfParticleWeights;};
  void SetqDistribution(TH1D* const qd) {this->fqDistribution = qd;};
  TH1D* GetqDistribution() const {return this->fqDistribution;};
  void SetqMin(Double_t const qmin) {this->fqMin = qmin;};
  Double_t GetqMin() const {return this->fqMin;};
  void SetqMax(Double_t const qmax) {this->fqMax = qmax;};
  Double_t GetqMax() const {return this->fqMax;};
  void SetqNbins(Int_t const qNbins) {this->fqNbins = qNbins;};
  Int_t GetqNbins() const {return this->fqNbins;};  
  // 4.) final results of fitting:
  void SetIntFlow(TH1D* const intFlow, Int_t sigmaFitted) {this->fIntFlow[sigmaFitted] = intFlow;};
  TH1D* GetIntFlow(Int_t sigmaFitted) const {return this->fIntFlow[sigmaFitted];};
  void SetSigma2(TH1D* const sigma2, Int_t sigmaFitted) {this->fSigma2[sigmaFitted] = sigma2;};
  TH1D* GetSigma2(Int_t sigmaFitted) const {return this->fSigma2[sigmaFitted];};
  void SetChi2(TH1D* const chi2, Int_t sigmaFitted) {this->fChi2[sigmaFitted] = chi2;};
  TH1D* GetChi2(Int_t sigmaFitted) const {return this->fChi2[sigmaFitted];};
  void SetFittingFunction(TF1* const ff, Int_t sigmaFitted) {this->fFittingFunction[sigmaFitted] = ff;};
  TF1* GetFittingFunction(Int_t sigmaFitted) const {return this->fFittingFunction[sigmaFitted];};
  // 5.) fitting parameters:
  void SetFittingParameters(TProfile* const fp) {this->fFittingParameters = fp;};
  TProfile* GetFittingParameters() const {return this->fFittingParameters;};
  void SetTreshold(Double_t const treshold) {this->fTreshold = treshold;};
  Double_t GetTreshold() const {return this->fTreshold;};
  void SetvStart(Double_t const vStart) {this->fvStart = vStart;};
  Double_t GetvStart() const {return this->fvStart;};
  void SetvMin(Double_t const vMin) {this->fvMin = vMin;};
  Double_t GetvMin() const {return this->fvMin;};
  void SetvMax(Double_t const vMax) {this->fvMax = vMax;};
  Double_t GetvMax() const {return this->fvMax;};
  void SetSigma2Start(Double_t const Sigma2Start) {this->fSigma2Start = Sigma2Start;};
  Double_t GetSigma2Start() const {return this->fSigma2Start;};
  void SetSigma2Min(Double_t const Sigma2Min) {this->fSigma2Min = Sigma2Min;};
  Double_t GetSigma2Min() const {return this->fSigma2Min;};
  void SetSigma2Max(Double_t const Sigma2Max) {this->fSigma2Max = Sigma2Max;};
  Double_t GetSigma2Max() const {return this->fSigma2Max;};
  void SetFinalResultIsFromSigma2Fitted(Bool_t frifs2f) {this->fFinalResultIsFromSigma2Fitted = frifs2f;};
  Bool_t GetFinalResultIsFromSigma2Fitted() const {return this->fFinalResultIsFromSigma2Fitted;};  
  void SetPrintOnTheScreen(Bool_t pots) {this->fPrintOnTheScreen = pots;};
  Bool_t GetPrintOnTheScreen() const {return this->fPrintOnTheScreen;};  
  
 private:
  AliFlowAnalysisWithFittingQDistribution(const AliFlowAnalysisWithFittingQDistribution &afawfqd);
  AliFlowAnalysisWithFittingQDistribution& operator=(const AliFlowAnalysisWithFittingQDistribution &afawfqd);           
  // 0.) base:
  TList *fHistList; // base list to hold all output object
  // 1.) common:
  AliFlowCommonHist *fCommonHists; // common control histograms 
  AliFlowCommonHistResults *fCommonHistsResults; // final results in common histograms
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
  // 2.) particle weights (abbreviated to 'pWeights' or even to 'pW' throughout the code):
  TList *fWeightsList; // list to hold all histograms with particle weights: fUseParticleWeights, fPhiWeights, fPtWeights and fEtaWeights
  Bool_t fUsePhiWeights; // use phi weights
  Bool_t fUsePtWeights; // use pt weights
  Bool_t fUseEtaWeights; // use eta weights
  TProfile *fUseParticleWeights; // profile with three bins to hold values of fUsePhiWeights, fUsePtWeights and fUseEtaWeights
  TH1F *fPhiWeights; // histogram holding phi weights
  TH1D *fPtWeights; // histogram holding pt weights
  TH1D *fEtaWeights; // histogram holding eta weights 
  // 3.) distributions:
  TH1D *fSumOfParticleWeights; // distribution of sum of particle weights (for unit weights this equals to multiplicity)
  TH1D *fqDistribution; // distribution of Q/sqrt{sum of phi weights}
  Double_t fqMin; // lower boundary of TH1D *fqDistribution
  Double_t fqMax; // upper boundary of TH1D *fqDistribution
  Int_t fqNbins; // number of bins of TH1D *fqDistribution
  // 4.) final results of fitting:
  TH1D *fIntFlow[2]; // final result for integrated flow [0=sigma^2 not fitted, 1=sigma^2 fitted]  
  TH1D *fSigma2[2]; // final results for sigma^2 [0=sigma^2 not fitted, 1=sigma^2 fitted]
  TH1D *fChi2[2]; // final results for chi^2 from Minuit [0=sigma^2 not fitted, 1=sigma^2 fitted]
  TF1 *fFittingFunction[2]; // resulting fitting function of q-distribution [0=sigma^2 not fitted, 1=sigma^2 fitted]
  // 5.) fitting parameters:
  TProfile *fFittingParameters; // profile to hold all fitting parameters
  Double_t fTreshold; // the first bin taken for the fitting is the first bin with nEntries >= fTreshold (analogously for the last bin)
  Double_t fvStart; // fitting of v will start from this point
  Double_t fvMin; // v range, lower boundary
  Double_t fvMax; // v range, upper boundary
  Double_t fSigma2Start; // fitting of sigma2 will start from this point
  Double_t fSigma2Min; // sigma2 range, lower boundary (this should be kept above 0.5 according to theorists...)
  Double_t fSigma2Max; // sigma2 range, upper boundary
  Bool_t fFinalResultIsFromSigma2Fitted; // the result obtained with sigma^2 fitted or sigma^2 fixed is being stored
  Bool_t fPrintOnTheScreen; // print or not the final results on the screen
  
  ClassDef(AliFlowAnalysisWithFittingQDistribution, 0);
};

//================================================================================================================

#endif





