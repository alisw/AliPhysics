/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/********************************************************** 
 * In this class azimuthal correlators in mixed harmonics *
 * are implemented in terms of Q-vectors. This approach   *
 * doesn't require evaluation of nested loops. This class *
 * can be used to:                                        *
 *                                                        *  
 *  a) Extract subdominant harmonics (like v1 and v4);    *
 *  b) Study flow of two-particle resonances;             *
 *  c) Study strong parity violation.                     * 
 *                                                        * 
 * Author: Ante Bilandzic (abilandzic@gmail.com)          *
 *********************************************************/ 

#ifndef ALIFLOWANALYSISWITHMIXEDHARMONICS_H
#define ALIFLOWANALYSISWITHMIXEDHARMONICS_H

#include "AliFlowCommonConstants.h" // needed as include
#include "TMatrixD.h"

class TDirectoryFile;
class TList;
class TFile;
class TH1F;
class TH1D;
class TH2;
class TH2D;
class TProfile;
class TProfile2D;

class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

//================================================================================================================

class AliFlowAnalysisWithMixedHarmonics
{
 public:
  AliFlowAnalysisWithMixedHarmonics();
  virtual ~AliFlowAnalysisWithMixedHarmonics(); 
  // 0.) Methods called in the constructor:
    virtual void InitializeArrays();
  // 1.) Method Init() and methods called within Init():
  virtual void Init();
    virtual void CrossCheckSettings();
    virtual void AccessConstants();
    virtual void BookAndNestAllLists();
    virtual void BookProfileHoldingSettings();
    virtual void BookCommonHistograms();
    virtual void BookAllEventByEventQuantities();
    virtual void BookAllAllEventQuantities();
    virtual void BookAndFillWeightsHistograms();
    // 2.) Method Make() and methods called within Make():
  virtual void Make(AliFlowEventSimple *anEvent);
    virtual void CheckPointersUsedInMake();
    virtual void Calculate3pCorrelator();
    virtual void CalculateNonIsotropicTerms();
    virtual void CalculateDifferential3pCorrelator();
    virtual void ResetEventByEventQuantities();
  // 3.) Method Finish() and methods called within Finish():
  virtual void Finish();  
    virtual void CheckPointersUsedInFinish(); 
    virtual void AccessSettings();       
    virtual void CorrectForDetectorEffects();
    virtual void PrintOnTheScreen();  
  // 4.) Method GetOutputHistograms and method called within it:
  virtual void GetOutputHistograms(TList *outputListHistos);
    virtual void GetPointersForBaseHistograms();
    virtual void GetPointersForCommonHistograms();
    virtual void GetPointersForAllEventProfiles();
    virtual void GetPointersForResultsHistograms();
  // 5.) Other methods:   
  virtual void WriteHistograms(TString outputFileName);
  virtual void WriteHistograms(TDirectoryFile *outputFileName);  
  // 6.) Setters and getters:
  void SetHistList(TList* const hl) {this->fHistList = hl;}
  TList* GetHistList() const {return this->fHistList;}  
  void SetHistListName(const char *hln) {this->fHistListName->Append(*hln);}; 
  TString *GetHistListName() const {return this->fHistListName;};
  void SetAnalysisLabel(const char *al) {this->fAnalysisLabel->Append(*al);}; 
  TString *GetAnalysisLabel() const {return this->fAnalysisLabel;};
  void SetAnalysisSettings(TProfile* const as) {this->fAnalysisSettings = as;};
  TProfile* GetAnalysisSettings() const {return this->fAnalysisSettings;};
  void SetCorrelatorInteger(Int_t const ci) {this->fCorrelatorInteger = ci;};
  Int_t GetCorrelatorInteger() const {return this->fCorrelatorInteger;}; 
  void SetNoOfMultipicityBins(Int_t const nomb) {this->fNoOfMultipicityBins = nomb;};
  Int_t GetNoOfMultipicityBins() const {return this->fNoOfMultipicityBins;};   
  void SetMultipicityBinWidth(Double_t const mbw) {this->fMultipicityBinWidth = mbw;};
  Double_t GetMultipicityBinWidth() const {return this->fMultipicityBinWidth;};   
  void SetMinMultiplicity(Double_t const mm) {this->fMinMultiplicity = mm;};
  Double_t GetMinMultiplicity() const {return this->fMinMultiplicity;}; 
  void SetCorrectForDetectorEffects(Bool_t const cfde) {this->fCorrectForDetectorEffects = cfde;};
  Bool_t GetCorrectForDetectorEffects() const {return this->fCorrectForDetectorEffects;}; 
  void SetPrintOnTheScreen(Bool_t const pots) {this->fPrintOnTheScreen = pots;};
  Bool_t GetPrintOnTheScreen() const {return this->fPrintOnTheScreen;};   
  void SetCommonHists(AliFlowCommonHist* const ch) {this->fCommonHists = ch;};
  AliFlowCommonHist* GetCommonHists() const {return this->fCommonHists;};
  void SetWeightsList(TList* const wl) {this->fWeightsList = (TList*)wl->Clone();}
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
  void SetProfileList(TList* const plist) {this->fProfileList = plist;}
  TList* GetProfileList() const {return this->fProfileList;}  
  void Set3pCorrelatorPro(TProfile* const s3pPro) {this->f3pCorrelatorPro = s3pPro;};
  TProfile* Get3pCorrelatorPro() const {return this->f3pCorrelatorPro;};
  void SetNonIsotropicTermsPro(TProfile* const nitPro) {this->fNonIsotropicTermsPro = nitPro;};
  TProfile* GetNonIsotropicTermsPro() const {return this->fNonIsotropicTermsPro;};
  void Set3pCorrelatorVsMPro(TProfile* const s3pVsMPro) {this->f3pCorrelatorVsMPro = s3pVsMPro;};
  TProfile* Get3pCorrelatorVsMPro() const {return this->f3pCorrelatorVsMPro;};
  void SetNonIsotropicTermsVsMPro(TProfile2D* const nitVsMPro) {this->fNonIsotropicTermsVsMPro = nitVsMPro;};
  TProfile2D* GetNonIsotropicTermsVsMPro() const {return this->fNonIsotropicTermsVsMPro;};
  void SetResultsList(TList* const rlist) {this->fResultsList = rlist;}
  TList* GetResultsList() const {return this->fResultsList;}    
  void Set3pCorrelatorHist(TH1D* const s3pHist) {this->f3pCorrelatorHist = s3pHist;};
  TH1D* Get3pCorrelatorHist() const {return this->f3pCorrelatorHist;};  
  void SetDetectorBiasHist(TH1D* const dbHist) {this->fDetectorBiasHist = dbHist;};
  TH1D* GetDetectorBiasHist() const {return this->fDetectorBiasHist;};  
  void SetDetectorBiasVsMHist(TH1D* const dbVsMHist) {this->fDetectorBiasVsMHist = dbVsMHist;};
  TH1D* GetDetectorBiasVsMHist() const {return this->fDetectorBiasVsMHist;};  
  void Set3pCorrelatorVsPtSumDiffPro(TProfile* const s3pcvpsd, Int_t const sd) {this->f3pCorrelatorVsPtSumDiffPro[sd] = s3pcvpsd;};
  TProfile* Get3pCorrelatorVsPtSumDiffPro(Int_t sd) const {return this->f3pCorrelatorVsPtSumDiffPro[sd];};
    
 private:
  AliFlowAnalysisWithMixedHarmonics(const AliFlowAnalysisWithMixedHarmonics& afawQc);
  AliFlowAnalysisWithMixedHarmonics& operator=(const AliFlowAnalysisWithMixedHarmonics& afawQc); 
  // 0.) Base:
  TList *fHistList; // base list to hold all output objects
  TString *fHistListName; // name of base list
  TString *fAnalysisLabel; // analysis label 
  TProfile *fAnalysisSettings; // profile to hold analysis settings
  Int_t fCorrelatorInteger; // integer n in cos[n(2phi1-phi2-phi3)]
  Int_t fNoOfMultipicityBins; // number of multiplicity bins
  Double_t fMultipicityBinWidth; // width of multiplicity bin
  Double_t fMinMultiplicity; // minimal multiplicity
  Bool_t fCorrectForDetectorEffects; // correct 3-p correlator for detector effects
  Bool_t fPrintOnTheScreen; // print or not the final results on the screen
  // 1.) Common:
  AliFlowCommonHist *fCommonHists; // common control histograms (filled only with events with 3 or more tracks for 3-p correlators) 
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
  // 2a.) Particle weights:
  TList *fWeightsList; // list to hold all histograms with particle weights: fUseParticleWeights, fPhiWeights, fPtWeights and fEtaWeights
  Bool_t fUsePhiWeights; // use phi weights
  Bool_t fUsePtWeights; // use pt weights
  Bool_t fUseEtaWeights; // use eta weights
  TProfile *fUseParticleWeights; // profile with three bins to hold values of fUsePhiWeights, fUsePtWeights and fUseEtaWeights
  TH1F *fPhiWeights; // histogram holding phi weights
  TH1D *fPtWeights; // histogram holding phi weights
  TH1D *fEtaWeights; // histogram holding phi weights 
  // 3.) Event-by-event quantities:
  TMatrixD *fReQnk; // fReQ[n][k] = Re[Q_{n,k}] = sum_{i=1}^{M} w_{i}^{k} cos(n*phi_{i})
  TMatrixD *fImQnk; // fImQ[n][k] = Im[Q_{n,k}] = sum_{i=1}^{M} w_{i}^{k} sin(n*phi_{i})
  TMatrixD *fSpk; // fS[p][k] = S_{p,k} = (sum_{i=1}^{M} w_{i}^{k})^{p+1} // note p+1 in the power to use 0th index in p in non-trivial way
  TH1D *f3pCorrelatorEBE; // 3-p correlator <cos[n(2phi1-phi2-phi3)]> for single event
  TH1D *fNonIsotropicTermsEBE; // correction terms to 3-p correlator <cos[n(2phi1-phi2-phi3)]> for single event
  TProfile *fRePEBE[2]; // real part of p_n vs [(p1+p2)/2,|p1-p2|]
  TProfile *fImPEBE[2]; // imaginary part of p_n vs [(p1+p2)/2,|p1-p2|]
  // 4.) Profiles:
  TList *fProfileList; // list holding all all-event profiles 
  TProfile *f3pCorrelatorPro; // 3-p correlator <<cos[n(2phi1-phi2-phi3)]>> not corrected for detector effects
  TProfile *fNonIsotropicTermsPro; // non-isotropic terms in the decomposition of 3-p correlator <<cos[n(2phi1-phi2-phi3)]>>
  TProfile *f3pCorrelatorVsMPro; // 3-p correlator <<cos[n(2phi1-phi2-phi3)]>> vs multiplicity
  TProfile2D *fNonIsotropicTermsVsMPro; // non-isotropic terms in the decomposition of <cos[n(2phi1-phi2-phi3)]> vs multiplicity
  TProfile *f3pCorrelatorVsPtSumDiffPro[2]; // differential 3-p correlator <<cos[n(2phi1-psi2-psi3)]>> vs [(p1+p2)/2,|p1-p2|]
  // 5.) Final results:
  TList *fResultsList; // list holding objects with final results 
  TH1D *f3pCorrelatorHist; // 3-p correlator <<cos[n(2phi1-phi2-phi3)]>> corrected for detector effects
  TH1D *fDetectorBiasHist; // bias comming from detector inefficiencies to 3-p correlator <<cos[n(2phi1-phi2-phi3)]>> (in %)
  TH1D *fDetectorBiasVsMHist; // bias comming from detector inefficiencies to 3-p correlator <<cos[n(2phi1-phi2-phi3)]>> (in %) versus multiplicity

  ClassDef(AliFlowAnalysisWithMixedHarmonics, 0);

};

//================================================================================================================

#endif





