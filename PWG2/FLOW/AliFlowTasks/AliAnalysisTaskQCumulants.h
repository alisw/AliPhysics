/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/**************************************
 * analysis task for Q-cumulants      * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/

#ifndef ALIANALYSISTASKQCUMULANTS_H
#define ALIANALYSISTASKQCUMULANTS_H

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class TList;
class AliFlowEventSimple;
class AliFlowAnalysisWithQCumulants;

//================================================================================================================

class AliAnalysisTaskQCumulants : public AliAnalysisTaskSE{
 public:
  AliAnalysisTaskQCumulants();
  AliAnalysisTaskQCumulants(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskQCumulants(){}; 
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  // Common:
  void SetBookOnlyBasicCCH(Bool_t const bobcch) {this->fBookOnlyBasicCCH = bobcch;};
  Bool_t GetBookOnlyBasicCCH() const {return this->fBookOnlyBasicCCH;};  
  void SetFillMultipleControlHistograms(Bool_t const fmch) {this->fFillMultipleControlHistograms = fmch;};
  Bool_t GetFillMultipleControlHistograms() const {return this->fFillMultipleControlHistograms;}; 
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  Int_t GetHarmonic() const {return this->fHarmonic;};
  void SetApplyCorrectionForNUA(Bool_t const applyCorrectionForNUA) {this->fApplyCorrectionForNUA = applyCorrectionForNUA;};
  Bool_t GetApplyCorrectionForNUA() const {return this->fApplyCorrectionForNUA;};
  void SetApplyCorrectionForNUAVsM(Bool_t const applyCorrectionForNUAVsM) {this->fApplyCorrectionForNUAVsM = applyCorrectionForNUAVsM;};
  Bool_t GetApplyCorrectionForNUAVsM() const {return this->fApplyCorrectionForNUAVsM;};      
  void SetPropagateErrorAlsoFromNIT(Bool_t const peafNIT) {this->fPropagateErrorAlsoFromNIT = peafNIT;};
  Bool_t GetPropagateErrorAlsoFromNIT() const {return this->fPropagateErrorAlsoFromNIT;};  
  void SetCalculateDiffFlow(Bool_t const calculateDiffFlow) {this->fCalculateDiffFlow = calculateDiffFlow;};
  Bool_t GetCalculateDiffFlow() const {return this->fCalculateDiffFlow;};
  void SetCalculate2DDiffFlow(Bool_t const calculate2DDiffFlow) {this->fCalculate2DDiffFlow = calculate2DDiffFlow;};
  Bool_t GetCalculate2DDiffFlow() const {return this->fCalculate2DDiffFlow;};
  void SetCalculateDiffFlowVsEta(Bool_t const cdfve) {this->fCalculateDiffFlowVsEta = cdfve;};
  Bool_t GetCalculateDiffFlowVsEta() const {return this->fCalculateDiffFlowVsEta;};  
  void SetStoreDistributions(Bool_t const storeDistributions) {this->fStoreDistributions = storeDistributions;};
  Bool_t GetStoreDistributions() const {return this->fStoreDistributions;};
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
 
  // Multiparticle correlations vs multiplicity:
  void SetnBinsMult(Int_t const nbm) {this->fnBinsMult = nbm;};
  Int_t GetnBinsMult() const {return this->fnBinsMult;};  
  void SetMinMult(Double_t const minm) {this->fMinMult = minm;};
  Double_t GetMinMult() const {return this->fMinMult;};
  void SetMaxMult(Double_t const maxm) {this->fMaxMult = maxm;};
  Double_t GetMaxMult() const {return this->fMaxMult;};
  // Particle weights:
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
  void SetUseTrackWeights(Bool_t const uTrackW) {this->fUseTrackWeights = uTrackW;};
  Bool_t GetUseTrackWeights() const {return this->fUseTrackWeights;};
  // Event weights:
  void SetMultiplicityWeight(const char *multiplicityWeight) {*this->fMultiplicityWeight = multiplicityWeight;};
 
 private:
  AliAnalysisTaskQCumulants(const AliAnalysisTaskQCumulants& aatqc);
  AliAnalysisTaskQCumulants& operator=(const AliAnalysisTaskQCumulants& aatqc);
  
  AliFlowEventSimple *fEvent;         // the input event
  AliFlowAnalysisWithQCumulants *fQC; // Q-cumulant object
  TList *fListHistos;                 // collection of output 
  // Common:
  Bool_t fBookOnlyBasicCCH;              // book only basis common control histrograms (by default book them all) 
  Bool_t fFillMultipleControlHistograms; // fill separately control histos for events with >= 2, 4, 6 and 8 particles
  Int_t fHarmonic;                       // harmonic  
  Bool_t fApplyCorrectionForNUA;         // apply correction for non-uniform acceptance 
  Bool_t fApplyCorrectionForNUAVsM;      // apply correction for non-uniform acceptance versus M    
  Bool_t fPropagateErrorAlsoFromNIT;     // propagate error by taking into account also non-isotrpic terms  
  Bool_t fCalculateDiffFlow;             // calculate differential flow in pt or eta
  Bool_t fCalculate2DDiffFlow;           // calculate differential flow in (pt,eta) (Remark: this is very expensive in terms of CPU time)
  Bool_t fCalculateDiffFlowVsEta;        // if you set kFALSE only differential flow vs pt is calculated  
  Bool_t fStoreDistributions;            // store or not distributions of correlations
  Bool_t fCalculateCumulantsVsM;         // calculate cumulants versus multiplicity  
  Bool_t fCalculateAllCorrelationsVsM;   // calculate all correlations versus multiplicity     
  Bool_t fMinimumBiasReferenceFlow;      // store as reference flow in AliFlowCommonHistResults the minimum bias result (kFALSE by default)     
  Bool_t fForgetAboutCovariances;        // when propagating error forget about the covariances  
  Bool_t fStorePhiDistributionForOneEvent; // store phi distribution for one event to illustrate flow
  Double_t fPhiDistributionForOneEventSettings[4]; // [v_min,v_max,refMult_min,refMult_max]        
  // Multiparticle correlations vs multiplicity:
  Int_t fnBinsMult;                   // number of multiplicity bins for flow analysis versus multiplicity  
  Double_t fMinMult;                  // minimal multiplicity for flow analysis versus multiplicity  
  Double_t fMaxMult;                  // maximal multiplicity for flow analysis versus multiplicity    
  // Particle weights:
  Bool_t fUseParticleWeights;         // use any particle weights
  Bool_t fUsePhiWeights;              // use phi weights
  Bool_t fUsePtWeights;               // use pt weights
  Bool_t fUseEtaWeights;              // use eta weights  
  Bool_t fUseTrackWeights;            // use track weights (e.g. VZERO sector weights)  
  TList *fWeightsList;                // list with weights
  // Event weights:
  TString *fMultiplicityWeight;       // event-by-event weights for multiparticle correlations ("combinations","unit" or "multiplicity")  
  
  ClassDef(AliAnalysisTaskQCumulants, 1); 
};

//================================================================================================================

#endif











