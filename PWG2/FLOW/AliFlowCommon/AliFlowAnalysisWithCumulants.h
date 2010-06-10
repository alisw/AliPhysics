/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************************* 
 * Flow analysis with cumulants. In this class   *
 * cumulants are calculated by making use of the *
 * formalism of generating functions proposed by *
 * Ollitrault et al.                             *
 *                                               * 
 *      Author: Ante Bilandzic                   * 
 *              (abilandzic@gmail.com)           *
 *************************************************/ 

#ifndef ALIFLOWANALYSISWITHCUMULANTS_H
#define ALIFLOWANALYSISWITHCUMULANTS_H

#include "AliFlowCommonConstants.h" 
#include "TMatrixD.h"

class TList;
class TFile;

class TH1D;
class TH1F;
class TProfile;
class TProfile2D;
class TProfile3D;
class TDirectoryFile;

class AliFlowEventSimple;
class AliFlowTrackSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;
class AliFlowVector;

//================================================================================================================

class AliFlowAnalysisWithCumulants{
 public:
  AliFlowAnalysisWithCumulants();
  virtual ~AliFlowAnalysisWithCumulants(); 
  // 0.) Methods called in the constructor:
    virtual void InitializeArrays();
  // 1.) Method Init() and methods called within Init():
  virtual void Init();
    virtual void CrossCheckSettings();
    virtual void AccessConstants();
    virtual void BookAndNestAllLists();
    virtual void BookProfileHoldingSettings();
    virtual void BookCommonHistograms();
    virtual void BookAndFillWeightsHistograms();
    virtual void BookEverythingForReferenceFlow();
    virtual void BookEverythingForDiffFlow();
    virtual void StoreReferenceFlowFlags();
    virtual void StoreDiffFlowFlags();  
    virtual void BookEverythingForTuning();  
  // 2.) Method Make() and methods called within Make():
  virtual void Make(AliFlowEventSimple* anEvent);
     virtual void CheckPointersUsedInMake(); 
     virtual void FillGeneratingFunctionForReferenceFlow(AliFlowEventSimple *anEvent);
     virtual void FillQvectorComponents(AliFlowEventSimple *anEvent);    
     virtual void FillGeneratingFunctionForDiffFlow(AliFlowEventSimple *anEvent); 
     virtual void FillGeneratingFunctionsForDifferentTuningParameters(AliFlowEventSimple *anEvent);     
  // 3.) Method Finish() and methods called within Finish():
  virtual void Finish();  
    virtual void CheckPointersUsedInFinish(); 
    virtual void AccessSettings();
    virtual void GetAvMultAndNoOfEvts();
    virtual void CalculateCumulantsForReferenceFlow();
    virtual void CalculateReferenceFlow();
    virtual void CalculateReferenceFlowError();
    virtual void FillCommonHistResultsForReferenceFlow();
    virtual void CalculateCumulantsForDiffFlow(TString rpPoi,TString ptEta);
    virtual void CalculateDifferentialFlow(TString rpPoi,TString ptEta);
    virtual void CalculateDifferentialFlowErrors(TString rpPoi,TString ptEta);
    virtual void FillCommonHistResultsForDifferentialFlow(TString rpPoi);
    virtual void CalculateIntegratedFlow(TString rpPoi); // to be improved (add also possibility to integrate over eta yield)        
    virtual void PrintFinalResults(TString rpPoi);
    virtual void FinalizeTuning();
  // 4.) Method GetOutputHistograms() and method called within it:
  virtual void GetOutputHistograms(TList *outputListHistos);
    virtual void GetPointersForBaseHistograms();
    virtual void GetPointersForCommonControlHistograms();
    virtual void GetPointersForCommonResultsHistograms();
    virtual void GetPointersForParticleWeightsHistograms();
    virtual void GetPointersForReferenceFlowObjects();
    virtual void GetPointersForDiffFlowObjects();
  virtual void WriteHistograms(TString *outputFileName);
  virtual void WriteHistograms(TString outputFileName);
  virtual void WriteHistograms(TDirectoryFile *outputFileName);
  // 5.) Other methods:   
  //     ...  
  // 6.) Setters and getters:
  //  6.0.) base:
  void SetHistList(TList* const hl) {this->fHistList = hl;}
  TList* GetHistList() const {return this->fHistList;}  
  void SetHistListName(const char *hln) {this->fHistListName->Append(*hln);}; 
  TString *GetHistListName() const {return this->fHistListName;};
  void SetAnalysisSettings(TProfile* const as) {this->fAnalysisSettings = as;};
  TProfile* GetAnalysisSettings() const {return this->fAnalysisSettings;};
  //  6.1.) common:
  void SetCommonHists(AliFlowCommonHist* const ch) {this->fCommonHists = ch;};
  AliFlowCommonHist* GetCommonHists() const {return this->fCommonHists;};
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
  void SetMultiple(Int_t const multiple) {this->fMultiple = multiple;};
  Int_t GetMultiple() const {return this->fMultiple;};    
  void SetR0(Double_t const r0) {this->fR0 = r0;};
  Double_t GetR0() const {return this->fR0;}; 
  void SetPrintFinalResults(Bool_t const printOrNot, Int_t const i) {this->fPrintFinalResults[i] = printOrNot;};
  Bool_t GetPrintFinalResults(Int_t i) const {return this->fPrintFinalResults[i];};   
  //  6.2.0.) particle weights:
  void SetWeightsList(TList* const wlist) {this->fWeightsList = (TList*)wlist->Clone();}
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
  //  6.2.1.) event weights:
  void SetMultiplicityWeight(const char *multiplicityWeight) {*this->fMultiplicityWeight = multiplicityWeight;};  
  //  6.3.) reference flow:
  void SetReferenceFlowFlags(TProfile* const rff) {this->fReferenceFlowFlags = rff;};
  TProfile* GetReferenceFlowFlags() const {return this->fReferenceFlowFlags;};
  void SetnBinsMult(Int_t const nbm) {this->fnBinsMult = nbm;};
  Int_t GetnBinsMult() const {return this->fnBinsMult;};  
  void SetMinMult(Double_t const minm) {this->fMinMult = minm;};
  Double_t GetMinMult() const {return this->fMinMult;};
  void SetMaxMult(Double_t const maxm) {this->fMaxMult = maxm;};
  Double_t GetMaxMult() const {return this->fMaxMult;};
  //  6.3.0.) profiles:
  void SetReferenceFlowGenFun(TProfile2D* const rfgf) {this->fReferenceFlowGenFun = rfgf;};
  TProfile2D* GetReferenceFlowGenFun() const {return this->fReferenceFlowGenFun;};  
  void SetQvectorComponents(TProfile* const qvc) {this->fQvectorComponents = qvc;};
  TProfile* GetQvectorComponents() const {return this->fQvectorComponents;};
  void SetAverageOfSquaredWeight(TProfile* const aosw) {this->fAverageOfSquaredWeight = aosw;};
  TProfile* GetSumOfSquaredWeight() const {return this->fAverageOfSquaredWeight;}; 
  //  6.3.1.) results: 
  void SetReferenceFlowCumulants(TH1D* const rfc) {this->fReferenceFlowCumulants = rfc;};
  TH1D* GetReferenceFlowCumulants() const {return this->fReferenceFlowCumulants;}; 
  void SetReferenceFlow(TH1D* const rf) {this->fReferenceFlow = rf;};
  TH1D* GetReferenceFlow() const {return this->fReferenceFlow;}; 
  void SetChi(TH1D* const c) {this->fChi = c;};
  TH1D* GetChi() const {return this->fChi;}; 
  // 6.4.) differential flow:
  void SetDiffFlowFlags(TProfile* const dff) {this->fDiffFlowFlags = dff;};
  TProfile* GetDiffFlowFlags() const {return this->fDiffFlowFlags;};
  //  6.4.0.) profiles: 
  void SetDiffFlowGenFun(TProfile3D* const dfgf, Int_t const ri, Int_t const rp, Int_t const pe) {this->fDiffFlowGenFun[ri][rp][pe] = dfgf;};
  TProfile3D* GetDiffFlowGenFun(Int_t const ri, Int_t const rp, Int_t const pe) const {return this->fDiffFlowGenFun[ri][rp][pe];};
  void SetNoOfParticlesInBin(TProfile* const nopib, Int_t const rp, Int_t const pe) {this->fNoOfParticlesInBin[rp][pe] = nopib;};
  TProfile* GetNoOfParticlesInBin(Int_t const rp, Int_t const pe) const {return this->fNoOfParticlesInBin[rp][pe];};
  //  6.4.1.) results:
  void SetDiffFlowCumulants(TH1D* const dfc, Int_t const rp, Int_t const pe, Int_t const co) {this->fDiffFlowCumulants[rp][pe][co] = dfc;};
  TH1D* GetDiffFlowCumulants(Int_t const rp, Int_t const pe, Int_t const co) const {return this->fDiffFlowCumulants[rp][pe][co];};
  void SetDiffFlow(TH1D* const df, Int_t const rp, Int_t const pe, Int_t const co) {this->fDiffFlow[rp][pe][co] = df;};
  TH1D* GetDiffFlow(Int_t const rp, Int_t const pe, Int_t const co) const {return this->fDiffFlow[rp][pe][co];};
  // 6.x.) Tuning the interpolating parameter r0 and using cutoff at different order in series:
  void SetTuningFlags(TProfile* const tf) {this->fTuningFlags = tf;};
  TProfile* GetTuningFlags() const {return this->fTuningFlags;};
  void SetTuneParameters(Bool_t const tp) {this->fTuneParameters = tp;};
  Bool_t GetTuneParameters() const {return this->fTuneParameters;};  
  void SetTuningR0(Double_t const tr0, Int_t const r) {this->fTuningR0[r] = tr0;};
  Double_t GetTuningR0(Int_t const r) const {return this->fTuningR0[r];};
  //  6.x.0.) profiles:
  void SetTuningGenFun(TProfile2D* const tgf, Int_t const r, Int_t const pq) {this->fTuningGenFun[r][pq] = tgf;};
  TProfile2D* GetTuningGenFun(Int_t const r, Int_t const pq) const {return this->fTuningGenFun[r][pq];};
  //  6.x.1.) results:  
  void SetTuningCumulants(TH1D* const tc, Int_t const r, Int_t const pq) {this->fTuningCumulants[r][pq] = tc;};
  TH1D* GetTuningCumulants(Int_t const r, Int_t const pq) const {return this->fTuningCumulants[r][pq];};
  void SetTuningFlow(TH1D* const tf, Int_t const r, Int_t const pq) {this->fTuningFlow[r][pq] = tf;};
  TH1D* GetTuningFlow(Int_t const r, Int_t const pq) const {return this->fTuningFlow[r][pq];};
  
 private:
  AliFlowAnalysisWithCumulants(const AliFlowAnalysisWithCumulants& afawc);
  AliFlowAnalysisWithCumulants& operator=(const AliFlowAnalysisWithCumulants& afawc); 
  // 0.) Base:
  TList *fHistList; // base list to hold all output objects
  TString *fHistListName; // name of base list  
  TProfile *fAnalysisSettings; // profile to hold analysis settings  
  // 1.) Common:
  AliFlowCommonHist *fCommonHists; // common control histograms (filled only with events with 3 or more tracks for 3-p correlators) 
  AliFlowCommonHistResults *fCommonHistsResults2nd; // common result histograms for 2nd order cumulant
  AliFlowCommonHistResults *fCommonHistsResults4th; // common result histograms for 4th order cumulant 
  AliFlowCommonHistResults *fCommonHistsResults6th; // common result histograms for 6th order cumulant
  AliFlowCommonHistResults *fCommonHistsResults8th; // common result histograms for 8th order cumulant
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
  Int_t fMultiple; // the multiple m in p=m*n, where n is harmonic (relevant for differential flow) 
  Double_t fR0; // r_{0} parameter
  Bool_t fPrintFinalResults[3]; // print on the screen the final results [0=RF,1=RP,2=POI]
  // 2a.) Particle weights:
  TList *fWeightsList; // list to hold all histograms with particle weights: fUseParticleWeights, fPhiWeights, fPtWeights and fEtaWeights
  Bool_t fUsePhiWeights; // use phi weights
  Bool_t fUsePtWeights; // use pt weights
  Bool_t fUseEtaWeights; // use eta weights
  TProfile *fUseParticleWeights; // profile with three bins to hold values of fUsePhiWeights, fUsePtWeights and fUseEtaWeights
  TH1F *fPhiWeights; // histogram holding phi weights
  TH1D *fPtWeights; // histogram holding phi weights
  TH1D *fEtaWeights; // histogram holding phi weights 
  // 2b.) Event weights:
  TString *fMultiplicityWeight; // event-by-event weight for reference flow generating function (can be "unit" or "multiplicity")
  // 3.) Reference flow:       
  //  3a.) lists:
  TList *fReferenceFlowList; // list to hold all histograms and profiles relevant for reference flow 
  TList *fReferenceFlowProfiles; // list to hold all profiles relevant for reference flow
  TList *fReferenceFlowResults; // list to hold all histograms with final results relevant for reference flow  
  //  3b.) flags:
  TProfile *fReferenceFlowFlags; // profile to hold all flags for reference flow
  Int_t fnBinsMult; // number of multiplicity bins for flow analysis versus multiplicity  
  Double_t fMinMult; // minimal multiplicity for flow analysis versus multiplicity  
  Double_t fMaxMult; // maximal multiplicity for flow analysis versus multiplicity  
  //  3c.) event-by-event quantities:
  TMatrixD *fGEBE; // reference flow generating function only for current event   
  //  3d.) profiles:
  TProfile2D *fReferenceFlowGenFun; // all-event average of the generating function used to calculate reference flow 
  TProfile *fQvectorComponents; // averages of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>)
  TProfile *fAverageOfSquaredWeight; // <<w^2>>, where w = wPhi*wPt*wEta       
  //  3e.) results:
  Double_t fAvM; // average multiplicity
  Int_t fnEvts; // number of events
  TH1D *fReferenceFlowCumulants; // final results for isotropic cumulants for reference flow   
  TH1D *fReferenceFlow; // final results for reference flow  
  TH1D *fChi; // final results for resolution 
  // 4.) Differential flow:       
  //  4a.) lists:
  TList *fDiffFlowList; // list to hold all histograms and profiles relevant for differential flow 
  TList *fDiffFlowProfiles; // list to hold all profiles relevant for differential flow
  TList *fDiffFlowResults; // list to hold all histograms with final results relevant for differential flow  
  //  4b.) flags:
  TProfile *fDiffFlowFlags; // profile to hold all flags for reference flow
  //  4c.) profiles:
  TProfile3D *fDiffFlowGenFun[2][2][2]; // all-event avarage of generating function used for differential flow [0=Re,1=Im][0=RP,1=POI][0=pt,1=eta]  
  TProfile *fNoOfParticlesInBin[2][2]; // number of particles in pt/eta bin for RPs/POIs [0=RP,1=POI][0=pt,1=eta]
  //  4d.) results:
  TH1D *fDiffFlowCumulants[2][2][4]; // isotropic differential flow cumulants [0=RP,1=POI][0=pt,1=eta][cumulant order]
  TH1D *fDiffFlow[2][2][4]; // differential flow [0=RP,1=POI][0=pt,1=eta][cumulant order]
  // x.) Tuning the interpolating parameter r0 and using cutoff at different order in series:
  //  xa.) lists:
  TList *fTuningList; // list to hold all histograms and profiles relevant for tuning 
  TList *fTuningProfiles; // list to hold all profiles relevant for tuning
  TList *fTuningResults; // list to hold all histograms with final results relevant for tuning  
  //  xb.) flags:
  TProfile *fTuningFlags; // profile to hold all flags for tuning
  Bool_t fTuneParameters; // tune r0 and cut series at different order
  Double_t fTuningR0[10]; // different r0 values (at maximum 10 different values allowed)
  //  xc.) profiles:
  TProfile2D *fTuningGenFun[10][5]; // generating function G evaluated for 10 different r0s and 5 different sets of (pmax,qmax)
  //  xd.) results:  
  TH1D *fTuningCumulants[10][5]; // isotropic cumulants for reference flow for 10 different r0s and 5 different sets of (pmax,qmax)
  TH1D *fTuningFlow[10][5]; // reference flow for 10 different r0s and 5 different sets of (pmax,qmax) 
    
  ClassDef(AliFlowAnalysisWithCumulants, 0);

};

//================================================================================================================

#endif





