/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/******************************** 
 * integrated flow estimate by  *
 *   fitting q-distribution     * 
 *                              *
 * author: Ante Bilandzic       * 
 *          (anteb@nikhef.nl)   *
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

class TH1;
class TProfile;

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
  // 2.) method Make() and methods called within Make(): 
  virtual void Make(AliFlowEventSimple* anEvent);
  // 3.) method Finish() and methods called within Finish(): 
  virtual void Finish(Bool_t doFit = kTRUE);
   virtual void DoFit(Bool_t useParticleWeights, Bool_t sigma2NotFixed);
   virtual void FillCommonHistResultsIntFlow(Bool_t useParticleWeights, Bool_t sigma2NotFixed);
   virtual void PrintFinalResultsForIntegratedFlow();
  // 4.) other methods:
  virtual void GetOutputHistograms(TList *outputListHistos); 
  virtual void WriteHistograms(TString *outputFileName);
  virtual void WriteHistograms(TString outputFileName);
    
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
  void SetSumOfParticleWeights(TH1D* const sopW, Int_t pW) {this->fSumOfParticleWeights[pW] = sopW;};
  TH1D* GetSumOfParticleWeights(Int_t pW) const {return this->fSumOfParticleWeights[pW];};
  void SetqDistribution(TH1D* const qd, Int_t pW) {this->fqDistribution[pW] = qd;};
  TH1D* GetqDistribution(Int_t pW) const {return this->fqDistribution[pW];};
  // 4.) final results of fitting:
  void SetIntFlow(TH1D* const intFlow, Int_t pW, Int_t sigmaFixed) {this->fIntFlow[pW][sigmaFixed] = intFlow;};
  TH1D* GetIntFlow(Int_t pW, Int_t sigmaFixed) const {return this->fIntFlow[pW][sigmaFixed];};
  void SetSigma2(TH1D* const sigma2, Int_t pW, Int_t sigmaFixed) {this->fSigma2[pW][sigmaFixed] = sigma2;};
  TH1D* GetSigma2(Int_t pW, Int_t sigmaFixed) const {return this->fSigma2[pW][sigmaFixed];};
  
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
  TH1D *fPtWeights; // histogram holding phi weights
  TH1D *fEtaWeights; // histogram holding phi weights 
  // 3.) distributions:
  TH1D *fSumOfParticleWeights[2]; // [0=particle weights are unit (not used), 1=particle weights are used]
  TH1D *fqDistribution[2]; // distribution of Q/sqrt{sum of particle weights} [0=particle weights are unit (not used), 1=particle weights are used]
  // 4.) final results of fitting:
  TH1D *fIntFlow[2][2]; // final result for integrated flow [0=pWeights are unit (not used), 1=pWeights are used][0=sigma^2 not fixed, 1=sigma^2 fixed to 1/2]  
  TH1D *fSigma2[2][2]; // final results for sigma^2 [0=pWeights are unit (not used), 1=pWeights are used][0=sigma^2 not fixed, 1=sigma^2 fixed to 1/2]
    
  ClassDef(AliFlowAnalysisWithFittingQDistribution, 0);
};

//================================================================================================================

#endif





