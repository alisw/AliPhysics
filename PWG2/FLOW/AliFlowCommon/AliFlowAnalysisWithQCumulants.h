/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/********************************** 
 * flow analysis with Q-cumulants * 
 *                                * 
 * author: Ante Bilandzic         * 
 *        (abilandzic@gmail.com)  *
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
class TDirectoryFile;

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
  virtual void InitializeArraysForVarious();
  virtual void InitializeArraysForNestedLoops();
  // 1.) method Init() and methods called within Init():
  virtual void Init();
    virtual void CrossCheckSettings();
    virtual void CommonConstants(TString method);
    virtual void BookAndNestAllLists();
      virtual void BookAndNestListsForDifferentialFlow();
    virtual void BookCommonHistograms();
    virtual void BookAndFillWeightsHistograms();
    virtual void BookEverythingForIntegratedFlow();
    virtual void BookEverythingForDifferentialFlow();
    virtual void BookEverythingFor2DDifferentialFlow();
    virtual void BookEverythingForDistributions(); 
    virtual void BookEverythingForVarious();
    virtual void BookEverythingForNestedLoops();   
    virtual void StoreIntFlowFlags();
    virtual void StoreDiffFlowFlags();
    virtual void StoreFlagsForDistributions();   
    virtual void StoreHarmonic();
  // 2.) method Make() and methods called within Make():
  virtual void Make(AliFlowEventSimple *anEvent);
    // 2a.) Common:
    virtual void CheckPointersUsedInMake();     
    virtual void FillAverageMultiplicities(Int_t nRP);
    virtual void FillCommonControlHistograms(AliFlowEventSimple *anEvent);
    virtual void ResetEventByEventQuantities();
    // 2b.) Reference flow:
    virtual void CalculateIntFlowCorrelations(); 
    virtual void CalculateIntFlowCorrelationsUsingParticleWeights();
    virtual void CalculateIntFlowProductOfCorrelations();
    virtual void CalculateIntFlowSumOfEventWeights();
    virtual void CalculateIntFlowSumOfProductOfEventWeights();
    virtual void CalculateIntFlowCorrectionsForNUACosTerms();
    virtual void CalculateIntFlowCorrectionsForNUACosTermsUsingParticleWeights();
    virtual void CalculateIntFlowCorrectionsForNUASinTerms();  
    virtual void CalculateIntFlowCorrectionsForNUASinTermsUsingParticleWeights();    
    virtual void CalculateIntFlowProductOfCorrectionTermsForNUA();
    virtual void CalculateIntFlowSumOfEventWeightsNUA();
    virtual void CalculateIntFlowSumOfProductOfEventWeightsNUA();
    // 2c.) Cross-checking reference flow correlations with nested loops: 
    virtual void EvaluateIntFlowNestedLoops(AliFlowEventSimple* const anEvent);
    virtual void EvaluateIntFlowCorrelationsWithNestedLoops(AliFlowEventSimple* const anEvent); 
    virtual void EvaluateIntFlowCorrelationsWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* const anEvent); 
    virtual void EvaluateIntFlowCorrectionsForNUAWithNestedLoops(AliFlowEventSimple* const anEvent); 
    virtual void EvaluateIntFlowCorrectionsForNUAWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* const anEvent);
    // 2d.) Differential flow:
    virtual void CalculateDiffFlowCorrelations(TString type, TString ptOrEta); // type = RP or POI
    virtual void CalculateDiffFlowCorrelationsUsingParticleWeights(TString type, TString ptOrEta); // type = RP or POI 
    virtual void CalculateDiffFlowProductOfCorrelations(TString type, TString ptOrEta); // type = RP or POI
    virtual void CalculateDiffFlowSumOfEventWeights(TString type, TString ptOrEta); // type = RP or POI
    virtual void CalculateDiffFlowSumOfProductOfEventWeights(TString type, TString ptOrEta); // type = RP or POI
    virtual void CalculateDiffFlowCorrectionsForNUACosTerms(TString type, TString ptOrEta);
    virtual void CalculateDiffFlowCorrectionsForNUACosTermsUsingParticleWeights(TString type, TString ptOrEta);
    virtual void CalculateDiffFlowCorrectionsForNUASinTerms(TString type, TString ptOrEta);  
    virtual void CalculateDiffFlowCorrectionsForNUASinTermsUsingParticleWeights(TString type, TString ptOrEta);  
    // 2e.) 2D differential flow:
    virtual void Calculate2DDiffFlowCorrelations(TString type); // type = RP or POI
    // 2f.) Other differential correlators (i.e. Teaney-Yan correlator):    
    virtual void CalculateOtherDiffCorrelators(TString type, TString ptOrEta); // type = RP or POI    
    // 2g.) Distributions of reference flow correlations:
    virtual void StoreDistributionsOfCorrelations();
    // 2h.) Store phi distibution for one event to vizualize flow:
    virtual void StorePhiDistributionForOneEvent(AliFlowEventSimple* const anEvent);    
    // 2i.) Cross-checking differential flow correlations with nested loops:
    virtual void EvaluateDiffFlowNestedLoops(AliFlowEventSimple* const anEvent);
    virtual void EvaluateDiffFlowCorrelationsWithNestedLoops(AliFlowEventSimple* const anEvent, TString type, TString ptOrEta);
    virtual void EvaluateDiffFlowCorrelationsWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* const anEvent, TString type, TString ptOrEta); 
    virtual void EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoops(AliFlowEventSimple* const anEvent, TString type, TString ptOrEta);
    virtual void EvaluateDiffFlowCorrectionTermsForNUAWithNestedLoopsUsingParticleWeights(AliFlowEventSimple* const anEvent, TString type, TString ptOrEta);
    virtual void EvaluateOtherDiffCorrelatorsWithNestedLoops(AliFlowEventSimple* const anEvent, TString type, TString ptOrEta);
  // 3.) method Finish() and methods called within Finish():
  virtual void Finish();
    virtual void CheckPointersUsedInFinish();     
    // 3a.) integrated flow:
    virtual void FinalizeCorrelationsIntFlow();
    virtual void FinalizeCorrectionTermsForNUAIntFlow(); 
    virtual void CalculateCovariancesIntFlow();  
    virtual void CalculateCovariancesNUAIntFlow();  
    virtual void CalculateCumulantsIntFlow(); 
    virtual void CalculateReferenceFlow(); 
    virtual void FillCommonHistResultsIntFlow();
    //  nua:   
    virtual void CalculateQcumulantsCorrectedForNUAIntFlow(); 
    virtual void PrintFinalResultsForIntegratedFlow(TString type);
    virtual void CrossCheckIntFlowCorrelations();
    virtual void CrossCheckIntFlowExtraCorrelations(); // extra correlations which appear only when particle weights are used
    virtual void CrossCheckIntFlowCorrectionTermsForNUA(); 
    // 3b.) differential flow:
    virtual void FinalizeReducedCorrelations(TString type, TString ptOrEta);
    virtual void CalculateDiffFlowCovariances(TString type, TString ptOrEta); 
    virtual void CalculateDiffFlowCumulants(TString type, TString ptOrEta); 
    virtual void CalculateDiffFlow(TString type, TString ptOrEta); 
    virtual void FinalizeCorrectionTermsForNUADiffFlow(TString type, TString ptOrEta); 
    virtual void CalculateDiffFlowCumulantsCorrectedForNUA(TString type, TString ptOrEta);   
    virtual void CalculateDiffFlowCorrectedForNUA(TString type, TString ptOrEta); 
    virtual void CalculateFinalResultsForRPandPOIIntegratedFlow(TString type); // to be improved (add also possibility to integrate over eta yield)
    virtual void FillCommonHistResultsDiffFlow(TString type);   
    virtual void CrossCheckDiffFlowCorrelations(TString type, TString ptOrEta);
    virtual void PrintNumberOfParticlesInSelectedBin();     
    virtual void CrossCheckDiffFlowCorrectionTermsForNUA(TString type, TString ptOrEta); 
    // 3c.) 2D:
    virtual void Calculate2DDiffFlowCumulants(TString type);    
    virtual void Calculate2DDiffFlow(TString type);    
    // 3d.) Other differential correlators:
    virtual void CrossCheckOtherDiffCorrelators(TString type, TString ptOrEta);
    
  // 4.)  method GetOutputHistograms() and methods called within GetOutputHistograms(): 
  virtual void GetOutputHistograms(TList *outputListHistos);
    virtual void GetPointersForCommonHistograms(); 
    virtual void GetPointersForParticleWeightsHistograms();
    virtual void GetPointersForIntFlowHistograms(); 
    virtual void GetPointersForDiffFlowHistograms(); 
    virtual void GetPointersFor2DDiffFlowHistograms(); 
    virtual void GetPointersForOtherDiffCorrelators(); 
    virtual void GetPointersForNestedLoopsHistograms(); 
    
  // 5.) other methods:   
  TProfile* MakePtProjection(TProfile2D *profilePtEta) const;
  TProfile* MakeEtaProjection(TProfile2D *profilePtEta) const;
  virtual void WriteHistograms(TString outputFileName);
  virtual void WriteHistograms(TDirectoryFile *outputFileName);
  
  // **** SETTERS and GETTERS ****
  
  // 0.) base:
  void SetHistList(TList* const hlist) {this->fHistList = hlist;} 
  TList* GetHistList() const {return this->fHistList;} 
  
  // 1.) common:
  void SetBookOnlyBasicCCH(Bool_t const bobcch) {this->fBookOnlyBasicCCH = bobcch;};
  Bool_t GetBookOnlyBasicCCH() const {return this->fBookOnlyBasicCCH;};  
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
  void SetCommonConstants(TProfile* const cc) {this->fCommonConstants = cc;};
  TProfile* GetCommonConstants() const {return this->fCommonConstants;};  
  void SetFillMultipleControlHistograms(Bool_t const fmch) {this->fFillMultipleControlHistograms = fmch;};
  Bool_t GetFillMultipleControlHistograms() const {return this->fFillMultipleControlHistograms;};  
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  Int_t GetHarmonic() const {return this->fHarmonic;};
  void SetAnalysisLabel(const char *aLabel) {this->fAnalysisLabel->Append(*aLabel);}; // to be improved (Append(*aLabel) changed into Append(aLabel)) 
  TString *GetAnalysisLabel() const {return this->fAnalysisLabel;};
  void SetPrintFinalResults(Bool_t const printOrNot, Int_t const i) {this->fPrintFinalResults[i] = printOrNot;};
  Bool_t GetPrintFinalResults(Int_t i) const {return this->fPrintFinalResults[i];};  
   
  // 2a.) particle weights:
  void SetWeightsList(TList* const wlist) {this->fWeightsList = (TList*)wlist->Clone();}
  TList* GetWeightsList() const {return this->fWeightsList;}  
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
  void SetUseTrackWeights(Bool_t const uTrackW) {this->fUseTrackWeights = uTrackW;};
  Bool_t GetUseTrackWeights() const {return this->fUseTrackWeights;};
  void SetUseParticleWeights(TProfile* const uPW) {this->fUseParticleWeights = uPW;};
  TProfile* GetUseParticleWeights() const {return this->fUseParticleWeights;};
  void SetPhiWeights(TH1F* const histPhiWeights) {this->fPhiWeights = histPhiWeights;};
  TH1F* GetPhiWeights() const {return this->fPhiWeights;};
  void SetPtWeights(TH1D* const histPtWeights) {this->fPtWeights = histPtWeights;};
  TH1D* GetPtWeights() const {return this->fPtWeights;};
  void SetEtaWeights(TH1D* const histEtaWeights) {this->fEtaWeights = histEtaWeights;};
  TH1D* GetEtaWeights() const {return this->fEtaWeights;};
  
  // 2b.) event weights:
  void SetMultiplicityWeight(const char *multiplicityWeight) {*this->fMultiplicityWeight = multiplicityWeight;};
  
  // 3.) Reference flow:
  // Flags:
  void SetIntFlowFlags(TProfile* const intFlowFlags) {this->fIntFlowFlags = intFlowFlags;};
  TProfile* GetIntFlowFlags() const {return this->fIntFlowFlags;};
  void SetApplyCorrectionForNUA(Bool_t const applyCorrectionForNUA) {this->fApplyCorrectionForNUA = applyCorrectionForNUA;};
  Bool_t GetApplyCorrectionForNUA() const {return this->fApplyCorrectionForNUA;};
  void SetApplyCorrectionForNUAVsM(Bool_t const applyCorrectionForNUAVsM) {this->fApplyCorrectionForNUAVsM = applyCorrectionForNUAVsM;};
  Bool_t GetApplyCorrectionForNUAVsM() const {return this->fApplyCorrectionForNUAVsM;};  
  void SetnBinsMult(Int_t const nbm) {this->fnBinsMult = nbm;};
  Int_t GetnBinsMult() const {return this->fnBinsMult;};  
  void SetMinMult(Double_t const minm) {this->fMinMult = minm;};
  Double_t GetMinMult() const {return this->fMinMult;};
  void SetMaxMult(Double_t const maxm) {this->fMaxMult = maxm;};
  Double_t GetMaxMult() const {return this->fMaxMult;};
  void SetPropagateErrorAlsoFromNIT(Bool_t const peafNIT) {this->fPropagateErrorAlsoFromNIT = peafNIT;};
  Bool_t GetPropagateErrorAlsoFromNIT() const {return this->fPropagateErrorAlsoFromNIT;};  
  void SetCalculateCumulantsVsM(Bool_t const ccvm) {this->fCalculateCumulantsVsM = ccvm;};
  Bool_t GetCalculateCumulantsVsM() const {return this->fCalculateCumulantsVsM;};   
  void SetCalculateAllCorrelationsVsM(Bool_t const cacvm) {this->fCalculateAllCorrelationsVsM = cacvm;};
  Bool_t GetCalculateAllCorrelationsVsM() const {return this->fCalculateAllCorrelationsVsM;};   
  void SetMinimumBiasReferenceFlow(Bool_t const mmrf) {this->fMinimumBiasReferenceFlow = mmrf;};
  Bool_t GetMinimumBiasReferenceFlow() const {return this->fMinimumBiasReferenceFlow;};  
  void SetForgetAboutCovariances(Bool_t const fac) {this->fForgetAboutCovariances = fac;};
  Bool_t GetForgetAboutCovariances() const {return this->fForgetAboutCovariances;};
  void SetStorePhiDistributionForOneEvent(Bool_t const spdfoe) {this->fStorePhiDistributionForOneEvent = spdfoe;};
  Bool_t GetStorePhiDistributionForOneEvent() const {return this->fStorePhiDistributionForOneEvent;};
  void SetPhiDistributionForOneEventSettings(Double_t const pdfoes, Int_t const i) {this->fPhiDistributionForOneEventSettings[i] = pdfoes;};
  Double_t GetPhiDistributionForOneEventSettings(Int_t const i) const {return this->fPhiDistributionForOneEventSettings[i];};

  // Reference flow profiles:
  void SetAvMultiplicity(TProfile* const avMultiplicity) {this->fAvMultiplicity = avMultiplicity;};
  TProfile* GetAvMultiplicity() const {return this->fAvMultiplicity;};
  void SetIntFlowCorrelationsPro(TProfile* const intFlowCorrelationsPro) {this->fIntFlowCorrelationsPro = intFlowCorrelationsPro;};
  TProfile* GetIntFlowCorrelationsPro() const {return this->fIntFlowCorrelationsPro;};
  void SetIntFlowSquaredCorrelationsPro(TProfile* const ifscp) {this->fIntFlowSquaredCorrelationsPro = ifscp;};
  TProfile* GetIntFlowSquaredCorrelationsPro() const {return this->fIntFlowSquaredCorrelationsPro;};
  void SetIntFlowCorrelationsVsMPro(TProfile* const ifcvp, Int_t const ci) {this->fIntFlowCorrelationsVsMPro[ci] = ifcvp;};
  TProfile* GetIntFlowCorrelationsVsMPro(Int_t const ci) const {return this->fIntFlowCorrelationsVsMPro[ci];};    
  void SetIntFlowSquaredCorrelationsVsMPro(TProfile* const ifscvp, Int_t const ci) {this->fIntFlowSquaredCorrelationsVsMPro[ci] = ifscvp;};
  TProfile* GetIntFlowSquaredCorrelationsVsMPro(Int_t const ci) const {return this->fIntFlowSquaredCorrelationsVsMPro[ci];};   
  void SetIntFlowCorrelationsAllPro(TProfile* const intFlowCorrelationsAllPro) {this->fIntFlowCorrelationsAllPro = intFlowCorrelationsAllPro;};
  TProfile* GetIntFlowCorrelationsAllPro() const {return this->fIntFlowCorrelationsAllPro;};  
  void SetIntFlowExtraCorrelationsPro(TProfile* const intFlowExtraCorrelationsPro) {this->fIntFlowExtraCorrelationsPro = intFlowExtraCorrelationsPro;};
  TProfile* GetIntFlowExtraCorrelationsPro() const {return this->fIntFlowExtraCorrelationsPro;};  
  void SetIntFlowProductOfCorrelationsPro(TProfile* const intFlowProductOfCorrelationsPro) {this->fIntFlowProductOfCorrelationsPro = intFlowProductOfCorrelationsPro;};
  TProfile* GetIntFlowProductOfCorrelationsPro() const {return this->fIntFlowProductOfCorrelationsPro;};      
  void SetIntFlowProductOfCorrelationsVsMPro(TProfile* const ifpocvm, Int_t const pi) {this->fIntFlowProductOfCorrelationsVsMPro[pi] = ifpocvm;};
  TProfile* GetIntFlowProductOfCorrelationsVsMPro(Int_t const pi) const {return this->fIntFlowProductOfCorrelationsVsMPro[pi];};    
  void SetIntFlowProductOfCorrectionTermsForNUAPro(TProfile* const ifpoctfNUA) {this->fIntFlowProductOfCorrectionTermsForNUAPro = ifpoctfNUA;};
  TProfile* GetIntFlowProductOfCorrectionTermsForNUAPro() const {return this->fIntFlowProductOfCorrectionTermsForNUAPro;};  
  void SetIntFlowCorrectionTermsForNUAPro(TProfile* const ifctfnp, Int_t const sc) {this->fIntFlowCorrectionTermsForNUAPro[sc] = ifctfnp;};
  TProfile* GetIntFlowCorrectionTermsForNUAPro(Int_t sc) const {return this->fIntFlowCorrectionTermsForNUAPro[sc];};    
  void SetIntFlowCorrectionTermsForNUAVsMPro(TProfile* const ifctfnpvm, Int_t const sc, Int_t const ci) {this->fIntFlowCorrectionTermsForNUAVsMPro[sc][ci] = ifctfnpvm;};
  TProfile* GetIntFlowCorrectionTermsForNUAVsMPro(Int_t sc, Int_t ci) const {return this->fIntFlowCorrectionTermsForNUAVsMPro[sc][ci];};    
  // integrated flow histograms holding all results:
  void SetIntFlowCorrelationsHist(TH1D* const intFlowCorrelationsHist) {this->fIntFlowCorrelationsHist = intFlowCorrelationsHist;};
  TH1D* GetIntFlowCorrelationsHist() const {return this->fIntFlowCorrelationsHist;};
  void SetIntFlowCorrelationsVsMHist(TH1D* const ifcvmh, Int_t const ci) {this->fIntFlowCorrelationsVsMHist[ci] = ifcvmh;};
  TH1D* GetIntFlowCorrelationsVsMHist(Int_t const ci) const {return this->fIntFlowCorrelationsVsMHist[ci];};    
  void SetIntFlowCorrelationsAllHist(TH1D* const intFlowCorrelationsAllHist) {this->fIntFlowCorrelationsAllHist = intFlowCorrelationsAllHist;};
  TH1D* GetIntFlowCorrelationsAllHist() const {return this->fIntFlowCorrelationsAllHist;};  
  void SetIntFlowCorrectionTermsForNUAHist(TH1D* const ifctfnh, Int_t const sc) {this->fIntFlowCorrectionTermsForNUAHist[sc] = ifctfnh;};
  TH1D* GetIntFlowCorrectionTermsForNUAHist(Int_t sc) const {return this->fIntFlowCorrectionTermsForNUAHist[sc];};  
  void SetIntFlowCovariances(TH1D* const intFlowCovariances) {this->fIntFlowCovariances = intFlowCovariances;};
  TH1D* GetIntFlowCovariances() const {return this->fIntFlowCovariances;};
  void SetIntFlowSumOfEventWeights(TH1D* const intFlowSumOfEventWeights, Int_t const power) {this->fIntFlowSumOfEventWeights[power] = intFlowSumOfEventWeights;};
  TH1D* GetIntFlowSumOfEventWeights(Int_t power) const {return this->fIntFlowSumOfEventWeights[power];};
  void SetIntFlowSumOfProductOfEventWeights(TH1D* const intFlowSumOfProductOfEventWeights) {this->fIntFlowSumOfProductOfEventWeights = intFlowSumOfProductOfEventWeights;};
  TH1D* GetIntFlowSumOfProductOfEventWeights() const {return this->fIntFlowSumOfProductOfEventWeights;}; 
  void SetIntFlowCovariancesVsM(TH1D* const ifcvm, Int_t ci) {this->fIntFlowCovariancesVsM[ci] = ifcvm;};
  TH1D* GetIntFlowCovariancesVsM(Int_t ci) const {return this->fIntFlowCovariancesVsM[ci];};    
  void SetIntFlowSumOfEventWeightsVsM(TH1D* const ifsoewvm, Int_t si, Int_t lc) {this->fIntFlowSumOfEventWeightsVsM[si][lc] = ifsoewvm;};
  TH1D* GetIntFlowSumOfEventWeightsVsM(Int_t si, Int_t lc) const {return this->fIntFlowSumOfEventWeightsVsM[si][lc];};    
  void SetIntFlowSumOfProductOfEventWeightsVsM(TH1D* const ifsopoevm, Int_t si) {this->fIntFlowSumOfProductOfEventWeightsVsM[si] = ifsopoevm;};
  TH1D* GetIntFlowSumOfProductOfEventWeightsVsM(Int_t si) const {return this->fIntFlowSumOfProductOfEventWeightsVsM[si];};      
  void SetIntFlowCovariancesNUA(TH1D* const intFlowCovariancesNUA) {this->fIntFlowCovariancesNUA = intFlowCovariancesNUA;};
  TH1D* GetIntFlowCovariancesNUA() const {return this->fIntFlowCovariancesNUA;};
  void SetIntFlowSumOfEventWeightsNUA(TH1D* const ifsoewNUA, Int_t const sc, Int_t const power) {this->fIntFlowSumOfEventWeightsNUA[sc][power] = ifsoewNUA;};
  TH1D* GetIntFlowSumOfEventWeightsNUA(Int_t sc, Int_t power) const {return this->fIntFlowSumOfEventWeightsNUA[sc][power];};
  void SetIntFlowSumOfProductOfEventWeightsNUA(TH1D* const ifsopoewNUA) {this->fIntFlowSumOfProductOfEventWeightsNUA = ifsopoewNUA;};
  TH1D* GetIntFlowSumOfProductOfEventWeightsNUA() const {return this->fIntFlowSumOfProductOfEventWeightsNUA;}; 
  void SetIntFlowQcumulants(TH1D* const intFlowQcumulants) {this->fIntFlowQcumulants = intFlowQcumulants;};
  TH1D* GetIntFlowQcumulants() const {return this->fIntFlowQcumulants;}; 
  void SetIntFlowQcumulantsVsM(TH1D* const intFlowQcumulantsVsM, Int_t co) {this->fIntFlowQcumulantsVsM[co] = intFlowQcumulantsVsM;};
  TH1D* GetIntFlowQcumulantsVsM(Int_t co) const {return this->fIntFlowQcumulantsVsM[co];};  
  void SetIntFlowQcumulantsRebinnedInM(TH1D* const ifqcrim) {this->fIntFlowQcumulantsRebinnedInM = ifqcrim;};
  TH1D* GetIntFlowQcumulantsRebinnedInM() const {return this->fIntFlowQcumulantsRebinnedInM;};    
  void SetIntFlowQcumulantsErrorSquaredRatio(TH1D* const ifqcesr) {this->fIntFlowQcumulantsErrorSquaredRatio = ifqcesr;};
  TH1D* GetIntFlowQcumulantsErrorSquaredRatio() const {return this->fIntFlowQcumulantsErrorSquaredRatio;}; 
  void SetIntFlow(TH1D* const intFlow) {this->fIntFlow = intFlow;};
  TH1D* GetIntFlow() const {return this->fIntFlow;};
  void SetIntFlowVsM(TH1D* const intFlowVsM, Int_t co) {this->fIntFlowVsM[co] = intFlowVsM;};
  TH1D* GetIntFlowVsM(Int_t co) const {return this->fIntFlowVsM[co];};     
  void SetIntFlowRebinnedInM(TH1D* const ifrim) {this->fIntFlowRebinnedInM = ifrim;};
  TH1D* GetIntFlowRebinnedInM() const {return this->fIntFlowRebinnedInM;};
  void SetIntFlowDetectorBias(TH1D* const ifdb) {this->fIntFlowDetectorBias = ifdb;};
  TH1D* GetIntFlowDetectorBias() const {return this->fIntFlowDetectorBias;};  
  void SetIntFlowDetectorBiasVsM(TH1D* const ifdbvm, Int_t ci) {this->fIntFlowDetectorBiasVsM[ci] = ifdbvm;};
  TH1D* GetIntFlowDetectorBiasVsM(Int_t ci) const {return this->fIntFlowDetectorBiasVsM[ci];};  
  // 4.) Differential flow:
  //  Flags:
  void SetDiffFlowFlags(TProfile* const diffFlowFlags) {this->fDiffFlowFlags = diffFlowFlags;};
  TProfile* GetDiffFlowFlags() const {return this->fDiffFlowFlags;};
  void SetCalculateDiffFlow(Bool_t const cdf) {this->fCalculateDiffFlow = cdf;};
  Bool_t GetCalculateDiffFlow() const {return this->fCalculateDiffFlow;};
  void SetCalculate2DDiffFlow(Bool_t const c2ddf) {this->fCalculate2DDiffFlow = c2ddf;};
  Bool_t GetCalculate2DDiffFlow() const {return this->fCalculate2DDiffFlow;};
  void SetCalculateDiffFlowVsEta(Bool_t const cdfve) {this->fCalculateDiffFlowVsEta = cdfve;};
  Bool_t GetCalculateDiffFlowVsEta() const {return this->fCalculateDiffFlowVsEta;};
  //  Profiles:
  //   1D:
  void SetDiffFlowCorrelationsPro(TProfile* const diffFlowCorrelationsPro, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowCorrelationsPro[i][j][k] = diffFlowCorrelationsPro;};
  TProfile* GetDiffFlowCorrelationsPro(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowCorrelationsPro[i][j][k];};
  void SetDiffFlowSquaredCorrelationsPro(TProfile* const diffFlowSquaredCorrelationsPro, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowSquaredCorrelationsPro[i][j][k] = diffFlowSquaredCorrelationsPro;};
  TProfile* GetDiffFlowSquaredCorrelationsPro(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowSquaredCorrelationsPro[i][j][k];}; 
  void SetDiffFlowProductOfCorrelationsPro(TProfile* const dfpocp, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowProductOfCorrelationsPro[i][j][k][l] = dfpocp;};
  TProfile* GetDiffFlowProductOfCorrelationsPro(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowProductOfCorrelationsPro[i][j][k][l];};
  void SetDiffFlowCorrectionTermsForNUAPro(TProfile* const dfctfnp, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowCorrectionTermsForNUAPro[i][j][k][l] = dfctfnp;};
  TProfile* GetDiffFlowCorrectionTermsForNUAPro(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowCorrectionTermsForNUAPro[i][j][k][l];};  
  //   2D:
  void Set2DDiffFlowCorrelationsPro(TProfile2D* const p2ddfcp, Int_t const i, Int_t const k) {this->f2DDiffFlowCorrelationsPro[i][k] = p2ddfcp;};
  TProfile2D* Get2DDiffFlowCorrelationsPro(Int_t i, Int_t k) const {return this->f2DDiffFlowCorrelationsPro[i][k];};
  //   Other differential correlators:
  void SetOtherDiffCorrelators(TProfile* const odc,Int_t const i,Int_t const j,Int_t const k,Int_t const l) {this->fOtherDiffCorrelators[i][j][k][l] = odc;};
  TProfile* GetOtherDiffCorrelators(Int_t i,Int_t j,Int_t k,Int_t l) const {return this->fOtherDiffCorrelators[i][j][k][l];};   
  // histograms:
  void SetDiffFlowCorrelationsHist(TH1D* const diffFlowCorrelationsHist, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowCorrelationsHist[i][j][k] = diffFlowCorrelationsHist;};
  TH1D* GetDiffFlowCorrelationsHist(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowCorrelationsHist[i][j][k];};
  void SetDiffFlowCovariances(TH1D* const diffFlowCovariances, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowCovariances[i][j][k] = diffFlowCovariances;};
  TH1D* GetDiffFlowCovariances(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowCovariances[i][j][k];};  
  void SetDiffFlowCumulants(TH1D* const diffFlowCumulants, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowCumulants[i][j][k] = diffFlowCumulants;};
  TH1D* GetDiffFlowCumulants(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowCumulants[i][j][k];};
  void SetDiffFlowDetectorBias(TH1D* const dfdb, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlowDetectorBias[i][j][k] = dfdb;};
  TH1D* GetDiffFlowDetectorBias(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowDetectorBias[i][j][k];};
  void SetDiffFlow(TH1D* const diffFlow, Int_t const i, Int_t const j, Int_t const k) {this->fDiffFlow[i][j][k] = diffFlow;};
  TH1D* GetDiffFlow(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlow[i][j][k];};
  void SetDiffFlowSumOfEventWeights(TH1D* const dfsoew, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowSumOfEventWeights[i][j][k][l] = dfsoew;};
  TH1D* GetDiffFlowSumOfEventWeights(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowSumOfEventWeights[i][j][k][l];};
  void SetDiffFlowSumOfProductOfEventWeights(TH1D* const dfsopoew, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowSumOfProductOfEventWeights[i][j][k][l] = dfsopoew;};
  TH1D* GetDiffFlowSumOfProductOfEventWeights(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowSumOfProductOfEventWeights[i][j][k][l];};
  void SetDiffFlowCorrectionTermsForNUAHist(TH1D* const dfctfnh, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowCorrectionTermsForNUAHist[i][j][k][l] = dfctfnh;};
  TH1D* GetDiffFlowCorrectionTermsForNUAHist(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowCorrectionTermsForNUAHist[i][j][k][l];};  
  //  2D:
  void Set2DDiffFlowCumulants(TH2D* const h2ddfc, Int_t const i, Int_t const j) {this->f2DDiffFlowCumulants[i][j] = h2ddfc;};
  TH2D* Get2DDiffFlowCumulants(Int_t i, Int_t j) const {return this->f2DDiffFlowCumulants[i][j];};  
  void Set2DDiffFlow(TH2D* const h2ddf, Int_t const i, Int_t const j) {this->f2DDiffFlow[i][j] = h2ddf;};
  TH2D* Get2DDiffFlow(Int_t i, Int_t j) const {return this->f2DDiffFlow[i][j];};  
  // 5.) distributions of correlations:
  // flags:
  void SetStoreDistributions(Bool_t const storeDistributions) {this->fStoreDistributions = storeDistributions;};
  Bool_t GetStoreDistributions() const {return this->fStoreDistributions;};
  // profile:
  void SetDistributionsFlags(TProfile* const distributionsFlags) {this->fDistributionsFlags = distributionsFlags;};
  TProfile* GetDistributionsFlags() const {return this->fDistributionsFlags;};  
  // histograms:
  void SetDistributions(TH1D* const distributions, Int_t const i) {this->fDistributions[i] = distributions;};
  TH1D* GetDistributions(Int_t i) const {return this->fDistributions[i];};  
  // min and max values of correlations (ci is correlations index [0=<2>,1=<4>,2=<6>,3=<8>]):
  void SetMinValueOfCorrelation(Int_t const ci, Double_t const minValue) {this->fMinValueOfCorrelation[ci] = minValue;};
  Double_t GetMinValueOfCorrelation(Int_t ci) const {return this->fMinValueOfCorrelation[ci];};
  void SetMaxValueOfCorrelation(Int_t const ci, Double_t const maxValue) {this->fMaxValueOfCorrelation[ci] = maxValue;};
  Double_t GetMaxValueOfCorrelation(Int_t ci) const {return this->fMaxValueOfCorrelation[ci];};
    
  // x.) debugging and cross-checking:
  void SetNestedLoopsList(TList* const nllist) {this->fNestedLoopsList = nllist;};
  TList* GetNestedLoopsList() const {return this->fNestedLoopsList;}; 
  void SetEvaluateIntFlowNestedLoops(Bool_t const eifnl) {this->fEvaluateIntFlowNestedLoops = eifnl;};
  Bool_t GetEvaluateIntFlowNestedLoops() const {return this->fEvaluateIntFlowNestedLoops;};
  void SetEvaluateDiffFlowNestedLoops(Bool_t const edfnl) {this->fEvaluateDiffFlowNestedLoops = edfnl;};
  Bool_t GetEvaluateDiffFlowNestedLoops() const {return this->fEvaluateDiffFlowNestedLoops;};  
  void SetMaxAllowedMultiplicity(Int_t const maxAllowedMultiplicity) {this->fMaxAllowedMultiplicity = maxAllowedMultiplicity;};
  Int_t GetMaxAllowedMultiplicity() const {return this->fMaxAllowedMultiplicity;};
  void SetEvaluateNestedLoops(TProfile* const enl) {this->fEvaluateNestedLoops = enl;};
  TProfile* GetEvaluateNestedLoops() const {return this->fEvaluateNestedLoops;}; 
  void SetIntFlowDirectCorrelations(TProfile* const ifdc) {this->fIntFlowDirectCorrelations = ifdc;};
  TProfile* GetIntFlowDirectCorrelations() const {return this->fIntFlowDirectCorrelations;};
  void SetIntFlowExtraDirectCorrelations(TProfile* const ifedc) {this->fIntFlowExtraDirectCorrelations = ifedc;};
  TProfile* GetIntFlowExtraDirectCorrelations() const {return this->fIntFlowExtraDirectCorrelations;};
  void SetIntFlowDirectCorrectionTermsForNUA(TProfile* const ifdctfn, Int_t const sc) {this->fIntFlowDirectCorrectionTermsForNUA[sc] = ifdctfn;};
  TProfile* GetIntFlowDirectCorrectionTermsForNUA(Int_t sc) const {return this->fIntFlowDirectCorrectionTermsForNUA[sc];};  
  void SetCrossCheckInPtBinNo(Int_t const crossCheckInPtBinNo) {this->fCrossCheckInPtBinNo = crossCheckInPtBinNo;};
  Int_t GetCrossCheckInPtBinNo() const {return this->fCrossCheckInPtBinNo;};
  void SetCrossCheckInEtaBinNo(Int_t const crossCheckInEtaBinNo) {this->fCrossCheckInEtaBinNo = crossCheckInEtaBinNo;};
  Int_t GetCrossCheckInEtaBinNo() const {return this->fCrossCheckInEtaBinNo;};
  void SetNoOfParticlesInBin(TH1D* const noOfParticlesInBin) {this->fNoOfParticlesInBin = noOfParticlesInBin;};
  TH1D* GetNoOfParticlesInBin() const {return this->fNoOfParticlesInBin;};  
  void SetDiffFlowDirectCorrelations(TProfile* const diffFlowDirectCorrelations,Int_t const i,Int_t const j,Int_t const k){this->fDiffFlowDirectCorrelations[i][j][k]=diffFlowDirectCorrelations;};
  TProfile* GetDiffFlowDirectCorrelations(Int_t i, Int_t j, Int_t k) const {return this->fDiffFlowDirectCorrelations[i][j][k];};
  void SetDiffFlowDirectCorrectionTermsForNUA(TProfile* const dfdctfn, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fDiffFlowDirectCorrectionTermsForNUA[i][j][k][l] = dfdctfn;};
  TProfile* GetDiffFlowDirectCorrectionTermsForNUA(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fDiffFlowDirectCorrectionTermsForNUA[i][j][k][l];};          
  void SetOtherDirectDiffCorrelators(TProfile* const oddc, Int_t const i, Int_t const j, Int_t const k, Int_t const l) {this->fOtherDirectDiffCorrelators[i][j][k][l] = oddc;};
  TProfile* GetOtherDirectDiffCorrelators(Int_t i, Int_t j, Int_t k, Int_t l) const {return this->fOtherDirectDiffCorrelators[i][j][k][l];};  

 private:
  
  AliFlowAnalysisWithQCumulants(const AliFlowAnalysisWithQCumulants& afawQc);
  AliFlowAnalysisWithQCumulants& operator=(const AliFlowAnalysisWithQCumulants& afawQc); 
  
  // 0.) base:
  TList* fHistList; // base list to hold all output object
  
  // 1.) common:
  Bool_t fBookOnlyBasicCCH; // book only basis common control histrograms (by default book them all)
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
  TProfile *fCommonConstants; // profile to hold common constants
  Bool_t fFillMultipleControlHistograms; // fill separately control histos for events with >= 2, 4, 6 and 8 particles 
  Int_t fHarmonic; // harmonic 
  TString *fAnalysisLabel; // analysis label (all histograms and output file will have this label)
  Bool_t fPrintFinalResults[4]; // print on the screen the final results (0=RF, 1=RP, 2=POI, 3=RF rebinned in M)
  
  // 2a.) particle weights:
  TList *fWeightsList; // list to hold all histograms with particle weights: fUseParticleWeights, fPhiWeights, fPtWeights and fEtaWeights
  Bool_t fUsePhiWeights; // use phi weights
  Bool_t fUsePtWeights; // use pt weights
  Bool_t fUseEtaWeights; // use eta weights
  Bool_t fUseTrackWeights; // use track weights (e.g. VZERO sector weights)
  TProfile *fUseParticleWeights; // profile with three bins to hold values of fUsePhiWeights, fUsePtWeights and fUseEtaWeights
  TH1F *fPhiWeights; // histogram holding phi weights
  TH1D *fPtWeights; // histogram holding phi weights
  TH1D *fEtaWeights; // histogram holding phi weights 
  
  // 2b.) event weights:
  TString *fMultiplicityWeight; // event-by-event weights for multiparticle correlations
  
  // 3.) integrated flow       
  //  3a.) lists:
  TList *fIntFlowList; // list to hold all histograms and profiles relevant for integrated flow 
  TList *fIntFlowProfiles; // list to hold all profiles relevant for integrated flow
  TList *fIntFlowResults; // list to hold all histograms with final results relevant for integrated flow  
  TList *fIntFlowAllCorrelationsVsM; // list to hold all profiles with correlations vs M
  //  3b.) flags:
  TProfile *fIntFlowFlags; // profile to hold all flags for integrated flow
  Bool_t fApplyCorrectionForNUA; // apply correction for non-uniform acceptance 
  Bool_t fApplyCorrectionForNUAVsM; // apply correction for non-uniform acceptance versus M  
  Int_t fnBinsMult; // number of multiplicity bins for flow analysis versus multiplicity  
  Double_t fMinMult; // minimal multiplicity for flow analysis versus multiplicity  
  Double_t fMaxMult; // maximal multiplicity for flow analysis versus multiplicity  
  Bool_t fPropagateErrorAlsoFromNIT; // propagate error by taking into account also non-isotropic terms (not sure if resulting error then is correct - to be improved)
  Bool_t fCalculateCumulantsVsM; // calculate cumulants versus multiplicity  
  Bool_t fCalculateAllCorrelationsVsM; // calculate all correlations versus multiplicity   
  Bool_t fMinimumBiasReferenceFlow; // store as reference flow in AliFlowCommonHistResults the minimum bias result (kFALSE by default)   
  Bool_t fForgetAboutCovariances; // when propagating error forget about the covariances  
  Bool_t fStorePhiDistributionForOneEvent; // store phi distribution for one event to illustrate flow
  Double_t fPhiDistributionForOneEventSettings[4]; // [v_min,v_max,refMult_min,refMult_max]
  //  3c.) event-by-event quantities:
  TMatrixD *fReQ; // fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
  TMatrixD *fImQ; // fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
  TMatrixD *fSpk; // fSM[p][k] = (sum_{i=1}^{M} w_{i}^{k})^{p+1}
  TH1D *fIntFlowCorrelationsEBE; // 1st bin: <2>, 2nd bin: <4>, 3rd bin: <6>, 4th bin: <8>
  TH1D *fIntFlowEventWeightsForCorrelationsEBE; // 1st bin: eW_<2>, 2nd bin: eW_<4>, 3rd bin: eW_<6>, 4th bin: eW_<8>
  TH1D *fIntFlowCorrelationsAllEBE; // to be improved (add comment)
  TH1D *fIntFlowCorrectionTermsForNUAEBE[2]; // [0=sin terms,1=cos terms], NUA = non-uniform acceptance
  TH1D *fIntFlowEventWeightForCorrectionTermsForNUAEBE[2]; // [0=sin terms,1=cos terms], NUA = non-uniform acceptance 
  Double_t fReferenceMultiplicityEBE; // reference multiplicity 
  //  3d.) profiles:
  TProfile *fAvMultiplicity; // profile to hold average multiplicities and number of events for events with nRP>=0, nRP>=1, ... , and nRP>=8
  TProfile *fIntFlowCorrelationsPro; // average correlations <<2>>, <<4>>, <<6>> and <<8>> (with wrong errors!) 
  TProfile *fIntFlowSquaredCorrelationsPro; // average correlations squared <<2>^2>, <<4>^2>, <<6>^2> and <<8>^2>  
  TProfile *fIntFlowCorrelationsVsMPro[4]; // average correlations <<2>>, <<4>>, <<6>> and <<8>> versus multiplicity (error is wrong here!)
  TProfile *fIntFlowSquaredCorrelationsVsMPro[4]; // average correlations <<2>^2>, <<4>^2>, <<6>^2> and <<8>^2> versus multiplicity  
  TProfile *fIntFlowCorrelationsAllPro; // average all correlations for integrated flow (with wrong errors!)
  TProfile *fIntFlowCorrelationsAllVsMPro[64]; // average all correlations vs M (errors via Sumw2 - to me improved)
  TProfile *fIntFlowExtraCorrelationsPro; // when particle weights are used some extra correlations appear 
  TProfile *fIntFlowProductOfCorrelationsPro; // average product of correlations <2>, <4>, <6> and <8>  
  TProfile *fIntFlowProductOfCorrelationsVsMPro[6]; // average product of correlations <2>, <4>, <6> and <8>  
                                                    // [0=<<2><4>>,1=<<2><6>>,2=<<2><8>>,3=<<4><6>>,4=<<4><8>>,5=<<6><8>>]  
  TProfile *fIntFlowProductOfCorrectionTermsForNUAPro; // average product of correction terms for NUA  
  TProfile *fIntFlowCorrectionTermsForNUAPro[2]; // average correction terms for non-uniform acceptance (with wrong errors!) [0=sin terms,1=cos terms] 
  TProfile *fIntFlowCorrectionTermsForNUAVsMPro[2][4]; // average correction terms for non-uniform acceptance (with wrong errors!) [0=sin terms,1=cos terms][correction term index] vs multiplicity   
  //  3e.) histograms with final results:
  TH1D *fIntFlowCorrelationsHist; // final results for average correlations <<2>>, <<4>>, <<6>> and <<8>> (with correct errors!) 
  TH1D *fIntFlowCorrelationsVsMHist[4]; // average correlations <<2>>, <<4>>, <<6>> and <<8>> versus multiplicity (error is correct here!)
  TH1D *fIntFlowCorrelationsAllHist; // final results for all average correlations (with correct errors!) 
  TH1D *fIntFlowCorrectionTermsForNUAHist[2];// final results for correction terms for non-uniform acceptance (with correct errors!) [0=sin terms,1=cos terms]
  TH1D *fIntFlowCovariances; // final result for covariances of correlations (multiplied with weight dependent prefactor)
  TH1D *fIntFlowSumOfEventWeights[2]; // sum of linear and quadratic event weights for <2>, <4>, <6> and <8>: [0=linear 1,1=quadratic]
  TH1D *fIntFlowSumOfProductOfEventWeights; // sum of products of event weights for correlations <2>, <4>, <6> and <8>  
  TH1D *fIntFlowCovariancesVsM[6]; // final result for covariances of correlations (multiplied with weight dependent prefactor) versus M
                                   // [0=Cov(2,4),1=Cov(2,6),2=Cov(2,8),3=Cov(4,6),4=Cov(4,8),5=Cov(6,8)]
  TH1D *fIntFlowSumOfEventWeightsVsM[4][2]; // sum of linear and quadratic event weights for <2>, <4>, <6> and <8> versum multiplicity
                                            // [0=sum{w_{<2>}},1=sum{w_{<4>}},2=sum{w_{<6>}},3=sum{w_{<8>}}][0=linear 1,1=quadratic]
  TH1D *fIntFlowSumOfProductOfEventWeightsVsM[6]; // sum of products of event weights for correlations <2>, <4>, <6> and <8> vs M
                                                  // [0=sum{w_{<2>}w_{<4>}},1=sum{w_{<2>}w_{<6>}},2=sum{w_{<2>}w_{<8>}},
                                                  //  3=sum{w_{<4>}w_{<6>}},4=sum{w_{<4>}w_{<8>}},5=sum{w_{<6>}w_{<8>}}]  
  TH1D *fIntFlowCovariancesNUA; // final result for covariances of all terms needed for NUA (multiplied with weight dependent prefactor)
  TH1D *fIntFlowSumOfEventWeightsNUA[2][2]; // sum of linear and quadratic event weights for NUA terms: [0=sin,1=cos][0=linear 1,1=quadratic]
  TH1D *fIntFlowSumOfProductOfEventWeightsNUA; // sum of products of event weights for NUA terms
  TH1D *fIntFlowQcumulants; // final results for integrated Q-cumulants QC{2}, QC{4}, QC{6} and QC{8}
  TH1D *fIntFlowQcumulantsVsM[4]; // final results for integrated Q-cumulants QC{2}, QC{4}, QC{6} and QC{8} versus multiplicity
  TH1D *fIntFlowQcumulantsRebinnedInM; // final results for reference Q-cumulants QC{2}, QC{4}, QC{6} and QC{8} rebinned in M
  TH1D *fIntFlowQcumulantsErrorSquaredRatio; // ratio between error squared: with/without non-isotropic terms
  TH1D *fIntFlow; // final results for integrated flow estimates v_n{2,QC}, v_n{4,QC}, v_n{6,QC} and v_n{8,QC}
  TH1D *fIntFlowVsM[4]; // final results for integrated flow estimates v_n{2,QC}, v_n{4,QC}, v_n{6,QC} and v_n{8,QC} versus multiplicity 
  TH1D *fIntFlowRebinnedInM; // final results for ref. flow estimates v_n{2,QC}, v_n{4,QC}, v_n{6,QC} and v_n{8,QC} rebinned in M
  TH1D *fIntFlowDetectorBias; // bias coming from detector inefficiencies to <<2>>, <<4>>, <<6>> and <<8>> (corrected/measured)  
  TH1D *fIntFlowDetectorBiasVsM[4]; // bias coming from detector inefficiencies to <<2>>, <<4>>, <<6>> and <<8>> vs M (corrected/measured)  
  // 4.) differential flow
  //  4a.) lists:
  TList *fDiffFlowList; // list to hold list with all histograms (fDiffFlowResults) and list with profiles (fDiffFlowProfiles) relevant for differential flow 
  TList *fDiffFlowProfiles; // list to hold all profiles relevant for differential flow
  TList *fDiffFlowResults; // list to hold all histograms with final results relevant for differential flow  
  TList *fDiffFlow2D; // list to hold all objects relevant for 2D differential flow  
  //    4aa.) nested list in list fDiffFlowProfiles: 
  TList *fDiffFlowCorrelationsProList[2][2]; // list to hold profiles with all correlations for differential flow [0=RP,1=POI][0=pt,1=eta] 
  TList *fDiffFlowProductOfCorrelationsProList[2][2]; // list to hold profiles with products of all correlations for differential flow [0=RP,1=POI][0=pt,1=eta] 
  TList *fDiffFlowCorrectionsProList[2][2]; // list to hold profiles with correction term for NUA for differential flow [0=RP,1=POI][0=pt,1=eta] 
  TList *f2DDiffFlowCorrelationsProList[2]; // list to hold profiles with all correlations for 2D differential flow [0=RP,1=POI]  
  //    4ab.) nested list in list fDiffFlowResults: 
  TList *fDiffFlowCorrelationsHistList[2][2]; // list to hold histograms with all correlations for differential flow [0=RP,1=POI][0=pt,1=eta] 
  TList *fDiffFlowSumOfEventWeightsHistList[2][2][2]; // list to hold histograms with sum of linear/quadratic event weights [0=RP,1=POI][0=pt,1=eta][0=linear 1,1=quadratic]
  TList *fDiffFlowSumOfProductOfEventWeightsHistList[2][2]; // list to hold histograms with sum of products of event weights [0=RP,1=POI][0=pt,1=eta]
  TList *fDiffFlowCorrectionsHistList[2][2]; // list to hold histograms with correction term for NUA for differential flow [0=RP,1=POI][0=pt,1=eta] 
  TList *fDiffFlowCovariancesHistList[2][2]; // list to hold histograms with all covariances for differential flow [0=RP,1=POI][0=pt,1=eta] 
  TList *fDiffFlowCumulantsHistList[2][2]; // list to hold histograms with all cumulants for differential flow [0=RP,1=POI][0=pt,1=eta] 
  TList *fDiffFlowDetectorBiasHistList[2][2]; // list to hold histograms which quantify detector bias to differential cumulants [0=RP,1=POI][0=pt,1=eta] 
  TList *fDiffFlowHistList[2][2]; // list to hold histograms with final results for differential flow [0=RP,1=POI][0=pt,1=eta]
  //  4b.) flags:  
  TProfile *fDiffFlowFlags; // profile to hold all flags for differential flow
  Bool_t fCalculateDiffFlow; // if you set kFALSE only reference flow will be calculated
  Bool_t fCalculate2DDiffFlow; // calculate 2D differential flow vs (pt,eta) (Remark: this is expensive in terms of CPU time)
  Bool_t fCalculateDiffFlowVsEta; // if you set kFALSE only differential flow vs pt is calculated
  //  4c.) event-by-event quantities:
  //   1D:
  TProfile *fReRPQ1dEBE[3][2][4][9]; // real part [0=r,1=p,2=q][0=pt,1=eta][m][k]
  TProfile *fImRPQ1dEBE[3][2][4][9]; // imaginary part [0=r,1=p,2=q][0=pt,1=eta][m][k]
  TProfile *fs1dEBE[3][2][9]; // [0=r,1=p,2=q][0=pt,1=eta][k] // to be improved
  TH1D *fDiffFlowCorrelationsEBE[2][2][4]; // [0=RP,1=POI][0=pt,1=eta][reduced correlation index]
  TH1D *fDiffFlowEventWeightsForCorrelationsEBE[2][2][4]; // [0=RP,1=POI][0=pt,1=eta][event weights for reduced correlation index]
  TH1D *fDiffFlowCorrectionTermsForNUAEBE[2][2][2][10]; // [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correction term index]
  //   2D:
  TProfile2D *fReRPQ2dEBE[3][4][9]; // real part of r_{m*n,k}(pt,eta), p_{m*n,k}(pt,eta) and q_{m*n,k}(pt,eta)
  TProfile2D *fImRPQ2dEBE[3][4][9]; // imaginary part of r_{m*n,k}(pt,eta), p_{m*n,k}(pt,eta) and q_{m*n,k}(pt,eta)
  TProfile2D *fs2dEBE[3][9]; // [t][k] // to be improved
  //  4d.) profiles:
  //   1D:
  TProfile *fDiffFlowCorrelationsPro[2][2][4]; // [0=RP,1=POI][0=pt,1=eta][correlation index]
  TProfile *fDiffFlowSquaredCorrelationsPro[2][2][4]; // [0=RP,1=POI][0=pt,1=eta][correlation index]
  TProfile *fDiffFlowProductOfCorrelationsPro[2][2][8][8]; // [0=RP,1=POI][0=pt,1=eta] [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>] x 
                                                           //                          [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>]
  TProfile *fDiffFlowCorrectionTermsForNUAPro[2][2][2][10]; // [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correction term index]
  //   2D:                                                            
  TProfile2D *f2DDiffFlowCorrelationsPro[2][4]; // [0=RP,1=POI][correlation index]
  //   Other differential correlators:
  TList *fOtherDiffCorrelatorsList; // list to hold profiles with other differential correlators
  TProfile *fOtherDiffCorrelators[2][2][2][1]; // // [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correlator index] 
  //  4e.) histograms holding final results:
  //   1D:
  TH1D *fDiffFlowCorrelationsHist[2][2][4]; // [0=RP,1=POI][0=pt,1=eta][correlation index]
  TH1D *fDiffFlowCovariances[2][2][5]; // [0=RP,1=POI][0=pW not used,1=pW used][0=exact eW,1=non-exact eW][0=pt,1=eta][index of covariances] 
  TH1D *fDiffFlowCumulants[2][2][4]; // [0=RP,1=POI][0=pt,1=eta][0=QC{2'},1=QC{4'},2=QC{6'},3=QC{8'}]
  TH1D *fDiffFlowDetectorBias[2][2][4]; // [0=RP,1=POI][0=pt,1=eta][0=gQC{2'}/QC{2'},1=gQC{4'}/QC{4'},2=gQC{6'}/QC{6'},3=gQC{8'}/QC{8'}]
  TH1D *fDiffFlow[2][2][4]; // [0=RP,1=POI][0=pt,1=eta][0=v'{2},1=v'{4},2=v'{6},3=v'{8}]
  TH1D *fDiffFlowSumOfEventWeights[2][2][2][4]; // [0=RP,1=POI][0=pt,1=eta][0=linear 1,1=quadratic][0=<2'>,1=<4'>,2=<6'>,3=<8'>]
  TH1D *fDiffFlowSumOfProductOfEventWeights[2][2][8][8]; // [0=RP,1=POI][0=pt,1=eta]  [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>] x 
                                                         //                           [0=<2>,1=<2'>,2=<4>,3=<4'>,4=<6>,5=<6'>,6=<8>,7=<8'>]
  TH1D *fDiffFlowCorrectionTermsForNUAHist[2][2][2][10]; // [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correction term index]        
  //   2D:                                                            
  TH2D *f2DDiffFlowCumulants[2][4]; // 2D differential cumulants [0=RP,1=POI][cumulant order]
  TH2D *f2DDiffFlow[2][4]; // 2D differential flow [0=RP,1=POI][cumulants order]
  // 5.) distributions:
  TList *fDistributionsList; // list to hold all distributions of correlations
  TProfile *fDistributionsFlags; // profile to hold all flags for distributions of correlations
  Bool_t fStoreDistributions; // store or not distributions of correlations
  TH1D *fDistributions[4]; // [0=distribution of <2>,1=distribution of <4>,2=distribution of <6>,3=distribution of <8>]
  Double_t fMinValueOfCorrelation[4]; // min values of <2>, <4>, <6> and <8>
  Double_t fMaxValueOfCorrelation[4]; // max values of <2>, <4>, <6> and <8>
  
  // 6.) various:
  TList *fVariousList; // list to hold various unclassified objects
  TH1D *fPhiDistributionForOneEvent; // store phi distribution for one event to illustrate flow
    
  // x.) debugging and cross-checking:
  TList *fNestedLoopsList; // list to hold all profiles filled with nested loops
  Bool_t fEvaluateIntFlowNestedLoops; // evaluate nested loops relevant for integrated flow
  Bool_t fEvaluateDiffFlowNestedLoops; // evaluate nested loops relevant for differential flow
  Int_t fMaxAllowedMultiplicity; // nested loops will be evaluated only for events with multiplicity <= fMaxAllowedMultiplicity
  TProfile *fEvaluateNestedLoops; // profile with four bins: fEvaluateIntFlowNestedLoops, fEvaluateDiffFlowNestedLoops, fCrossCheckInPtBinNo and fCrossCheckInEtaBinNo 
  // integrated flow:
  TProfile *fIntFlowDirectCorrelations; // multiparticle correlations relevant for int. flow calculated with nested loops  
  TProfile *fIntFlowExtraDirectCorrelations; // when particle weights are used some extra correlations appear   
  TProfile *fIntFlowDirectCorrectionTermsForNUA[2]; // average correction terms for non-uniform acceptance evaluated with nested loops [0=sin terms,1=cos terms] 
  // differential flow:
  Int_t fCrossCheckInPtBinNo; // cross-check results for reduced correlations and corrections in this pt bin
  Int_t fCrossCheckInEtaBinNo; // cross-check results for reduced correlations and corrections in this eta bin
  TH1D *fNoOfParticlesInBin; // bin: 1 = # of RPs in pt bin, 2 = # of RPs in eta bin, 3 = # of POIs in pt bin, 4 = # of POIs in eta bin 
  TProfile *fDiffFlowDirectCorrelations[2][2][4]; // [0=RP,1=POI][0=pt,1=eta][correlation index]
  TProfile *fDiffFlowDirectCorrectionTermsForNUA[2][2][2][10]; // [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correction term index]
  // other differential correlators: 
  TProfile *fOtherDirectDiffCorrelators[2][2][2][1]; // [0=RP,1=POI][0=pt,1=eta][0=sin terms,1=cos terms][correlator index]
                  
  ClassDef(AliFlowAnalysisWithQCumulants, 0);
};

//================================================================================================================

#endif





