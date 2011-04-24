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
 *     Author: Ante Bilandzic (abilandzic@gmail.com)      *
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
  virtual void BookDefault(); // book histos and profiles without any binning in multiplicity, pt or eta
  virtual void BookVsM();
  virtual void BookDifferential();
  virtual void BookAndFillWeightsHistograms();
  virtual void StoreHarmonic();    
  
  // 2.) Method Make() and methods called within Make():
  virtual void Make(AliFlowEventSimple *anEvent);
  virtual void CheckPointersUsedInMake();
  virtual void Calculate3pCorrelator();
  virtual void Calculate5pCorrelator();
  virtual void CalculateNonIsotropicTerms();
  virtual void CalculateDifferential3pCorrelator(Double_t &gIntegratedValue);
						 
  virtual void ResetEventByEventQuantities();
  
  // 3.) Method Finish() and methods called within Finish():
  virtual void Finish();  
  virtual void AccessSettings();       
  virtual void CheckPointersUsedInFinish(); 
  virtual void CorrectForDetectorEffects();
  virtual void CorrectForDetectorEffectsVsM();
  virtual void PrintOnTheScreen();  
  virtual void GetCorrelatorAndError(TProfile *g3pCorrelatorVsPt, 
				     Double_t &g3pCorrelatorValue, 
				     Double_t &g3pCorrelatorError);

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
  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
  Int_t GetHarmonic() const {return this->fHarmonic;};  
  void SetAnalysisLabel(const char *al) {this->fAnalysisLabel->Append(*al);}; 
  TString *GetAnalysisLabel() const {return this->fAnalysisLabel;};
  void SetAnalysisSettings(TProfile* const as) {this->fAnalysisSettings = as;};
  TProfile* GetAnalysisSettings() const {return this->fAnalysisSettings;};
  void SetNoOfMultipicityBins(Int_t const nomb) {this->fNoOfMultipicityBins = nomb;};
  Int_t GetNoOfMultipicityBins() const {return this->fNoOfMultipicityBins;};   
  void SetMultipicityBinWidth(Double_t const mbw) {this->fMultipicityBinWidth = mbw;};
  Double_t GetMultipicityBinWidth() const {return this->fMultipicityBinWidth;};   
  void SetMinMultiplicity(Double_t const mm) {this->fMinMultiplicity = mm;};
  Double_t GetMinMultiplicity() const {return this->fMinMultiplicity;}; 
  void SetOppositeChargesPOI(Bool_t const ocp) {this->fOppositeChargesPOI = ocp;};
  Bool_t GetOppositeChargesPOI() const {return this->fOppositeChargesPOI;};   
  void SetEvaluateDifferential3pCorrelator(Bool_t const ed3pc) {this->fEvaluateDifferential3pCorrelator = ed3pc;};
  Bool_t GetEvaluateDifferential3pCorrelator() const {return this->fEvaluateDifferential3pCorrelator;}; 
  void SetCorrectForDetectorEffects(Bool_t const cfde) {this->fCorrectForDetectorEffects = cfde;};
  Bool_t GetCorrectForDetectorEffects() const {return this->fCorrectForDetectorEffects;}; 
  void SetPrintOnTheScreen(Bool_t const pots) {this->fPrintOnTheScreen = pots;};
  Bool_t GetPrintOnTheScreen() const {return this->fPrintOnTheScreen;};  
  void SetCalculateVsM(Bool_t const cvm) {this->fCalculateVsM = cvm;};
  Bool_t GetCalculateVsM() const {return this->fCalculateVsM;};  
  void SetShowBinLabelsVsM(Bool_t const sblvm) {this->fShowBinLabelsVsM = sblvm;};
  Bool_t GetShowBinLabelsVsM() const {return this->fShowBinLabelsVsM;};  
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
  void Set5pCorrelatorPro(TProfile* const s5pPro) {this->f5pCorrelatorPro = s5pPro;};
  TProfile* Get5pCorrelatorPro() const {return this->f5pCorrelatorPro;};
  void SetNonIsotropicTermsPro(TProfile* const nitPro) {this->fNonIsotropicTermsPro = nitPro;};
  TProfile* GetNonIsotropicTermsPro() const {return this->fNonIsotropicTermsPro;};
  void Set3pCorrelatorVsMPro(TProfile* const s3pVsMPro) {this->f3pCorrelatorVsMPro = s3pVsMPro;};
  TProfile* Get3pCorrelatorVsMPro() const {return this->f3pCorrelatorVsMPro;};
  void Set3pPOICorrelatorVsM(TProfile* const s3pPOIVsM) {this->f3pPOICorrelatorVsM = s3pPOIVsM;};
  TProfile* Get3pPOICorrelatorVsM() const {return this->f3pPOICorrelatorVsM;};
  void SetNonIsotropicTermsVsMPro(TProfile2D* const nitVsMPro) {this->fNonIsotropicTermsVsMPro = nitVsMPro;};
  TProfile2D* GetNonIsotropicTermsVsMPro() const {return this->fNonIsotropicTermsVsMPro;};

  //2p correlators - pt diff
  void Set2pCorrelatorCosPsiDiffPtDiff(TProfile* const g2pCorrelatorCosPsiDiffPtDiff) {this->f2pCorrelatorCosPsiDiffPtDiff = g2pCorrelatorCosPsiDiffPtDiff;};
  TProfile* Get2pCorrelatorCosPsiDiffPtDiff() const {return this->f2pCorrelatorCosPsiDiffPtDiff;};
  void Set2pCorrelatorCosPsiSumPtDiff(TProfile* const g2pCorrelatorCosPsiSumPtDiff) {this->f2pCorrelatorCosPsiSumPtDiff = g2pCorrelatorCosPsiSumPtDiff;};
  TProfile* Get2pCorrelatorCosPsiSumPtDiff() const {return this->f2pCorrelatorCosPsiSumPtDiff;};
  void Set2pCorrelatorSinPsiDiffPtDiff(TProfile* const g2pCorrelatorSinPsiDiffPtDiff) {this->f2pCorrelatorSinPsiDiffPtDiff = g2pCorrelatorSinPsiDiffPtDiff;};
  TProfile* Get2pCorrelatorSinPsiDiffPtDiff() const {return this->f2pCorrelatorSinPsiDiffPtDiff;};
  void Set2pCorrelatorSinPsiSumPtDiff(TProfile* const g2pCorrelatorSinPsiSumPtDiff) {this->f2pCorrelatorSinPsiSumPtDiff = g2pCorrelatorSinPsiSumPtDiff;};
  TProfile* Get2pCorrelatorSinPsiSumPtDiff() const {return this->f2pCorrelatorSinPsiSumPtDiff;};

  //2p correlators - pt sum
  void Set2pCorrelatorCosPsiDiffPtSum(TProfile* const g2pCorrelatorCosPsiDiffPtSum) {this->f2pCorrelatorCosPsiDiffPtSum = g2pCorrelatorCosPsiDiffPtSum;};
  TProfile* Get2pCorrelatorCosPsiDiffPtSum() const {return this->f2pCorrelatorCosPsiDiffPtSum;};
  void Set2pCorrelatorCosPsiSumPtSum(TProfile* const g2pCorrelatorCosPsiSumPtSum) {this->f2pCorrelatorCosPsiSumPtSum = g2pCorrelatorCosPsiSumPtSum;};
  TProfile* Get2pCorrelatorCosPsiSumPtSum() const {return this->f2pCorrelatorCosPsiSumPtSum;};
  void Set2pCorrelatorSinPsiDiffPtSum(TProfile* const g2pCorrelatorSinPsiDiffPtSum) {this->f2pCorrelatorSinPsiDiffPtSum = g2pCorrelatorSinPsiDiffPtSum;};
  TProfile* Get2pCorrelatorSinPsiDiffPtSum() const {return this->f2pCorrelatorSinPsiDiffPtSum;};
  void Set2pCorrelatorSinPsiSumPtSum(TProfile* const g2pCorrelatorSinPsiSumPtSum) {this->f2pCorrelatorSinPsiSumPtSum = g2pCorrelatorSinPsiSumPtSum;};
  TProfile* Get2pCorrelatorSinPsiSumPtSum() const {return this->f2pCorrelatorSinPsiSumPtSum;};

  //2p correlators - eta diff
  void Set2pCorrelatorCosPsiDiffEtaDiff(TProfile* const g2pCorrelatorCosPsiDiffEtaDiff) {this->f2pCorrelatorCosPsiDiffEtaDiff = g2pCorrelatorCosPsiDiffEtaDiff;};
  TProfile* Get2pCorrelatorCosPsiDiffEtaDiff() const {return this->f2pCorrelatorCosPsiDiffEtaDiff;};
  void Set2pCorrelatorCosPsiSumEtaDiff(TProfile* const g2pCorrelatorCosPsiSumEtaDiff) {this->f2pCorrelatorCosPsiSumEtaDiff = g2pCorrelatorCosPsiSumEtaDiff;};
  TProfile* Get2pCorrelatorCosPsiSumEtaDiff() const {return this->f2pCorrelatorCosPsiSumEtaDiff;};
  void Set2pCorrelatorSinPsiDiffEtaDiff(TProfile* const g2pCorrelatorSinPsiDiffEtaDiff) {this->f2pCorrelatorSinPsiDiffEtaDiff = g2pCorrelatorSinPsiDiffEtaDiff;};
  TProfile* Get2pCorrelatorSinPsiDiffEtaDiff() const {return this->f2pCorrelatorSinPsiDiffEtaDiff;};
  void Set2pCorrelatorSinPsiSumEtaDiff(TProfile* const g2pCorrelatorSinPsiSumEtaDiff) {this->f2pCorrelatorSinPsiSumEtaDiff = g2pCorrelatorSinPsiSumEtaDiff;};
  TProfile* Get2pCorrelatorSinPsiSumEtaDiff() const {return this->f2pCorrelatorSinPsiSumEtaDiff;};

  //2p correlators - eta sum
  void Set2pCorrelatorCosPsiDiffEtaSum(TProfile* const g2pCorrelatorCosPsiDiffEtaSum) {this->f2pCorrelatorCosPsiDiffEtaSum = g2pCorrelatorCosPsiDiffEtaSum;};
  TProfile* Get2pCorrelatorCosPsiDiffEtaSum() const {return this->f2pCorrelatorCosPsiDiffEtaSum;};
  void Set2pCorrelatorCosPsiSumEtaSum(TProfile* const g2pCorrelatorCosPsiSumEtaSum) {this->f2pCorrelatorCosPsiSumEtaSum = g2pCorrelatorCosPsiSumEtaSum;};
  TProfile* Get2pCorrelatorCosPsiSumEtaSum() const {return this->f2pCorrelatorCosPsiSumEtaSum;};
  void Set2pCorrelatorSinPsiDiffEtaSum(TProfile* const g2pCorrelatorSinPsiDiffEtaSum) {this->f2pCorrelatorSinPsiDiffEtaSum = g2pCorrelatorSinPsiDiffEtaSum;};
  TProfile* Get2pCorrelatorSinPsiDiffEtaSum() const {return this->f2pCorrelatorSinPsiDiffEtaSum;};
  void Set2pCorrelatorSinPsiSumEtaSum(TProfile* const g2pCorrelatorSinPsiSumEtaSum) {this->f2pCorrelatorSinPsiSumEtaSum = g2pCorrelatorSinPsiSumEtaSum;};
  TProfile* Get2pCorrelatorSinPsiSumEtaSum() const {return this->f2pCorrelatorSinPsiSumEtaSum;};

  void SetResultsList(TList* const rlist) {this->fResultsList = rlist;}
  TList* GetResultsList() const {return this->fResultsList;}    
  void Set3pCorrelatorHist(TH1D* const s3pHist) {this->f3pCorrelatorHist = s3pHist;};
  TH1D* Get3pCorrelatorHist() const {return this->f3pCorrelatorHist;};    
  void Set3pCorrelatorVsMHist(TH1D* const s3pVsMHist) {this->f3pCorrelatorVsMHist = s3pVsMHist;};
  TH1D* Get3pCorrelatorVsMHist() const {return this->f3pCorrelatorVsMHist;};
  void SetDetectorBiasHist(TH1D* const dbHist) {this->fDetectorBiasHist = dbHist;};
  TH1D* GetDetectorBiasHist() const {return this->fDetectorBiasHist;};  
  void SetDetectorBiasVsMHist(TH1D* const dbVsMHist) {this->fDetectorBiasVsMHist = dbVsMHist;};
  TH1D* GetDetectorBiasVsMHist() const {return this->fDetectorBiasVsMHist;};  
  void Set3pCorrelatorVsPtSumDiffPro(TProfile* const s3pcvpsd, Int_t const sd) {this->f3pCorrelatorVsPtSumDiffPro[sd] = s3pcvpsd;};
  TProfile* Get3pCorrelatorVsPtSumDiffPro(Int_t sd) const {return this->f3pCorrelatorVsPtSumDiffPro[sd];};
  void Set3pCorrelatorVsEtaSumDiffPro(TProfile* const s3pcvpsd, Int_t const sd) {this->f3pCorrelatorVsEtaSumDiffPro[sd] = s3pcvpsd;};
  TProfile* Get3pCorrelatorVsEtaSumDiffPro(Int_t sd) const {return this->f3pCorrelatorVsEtaSumDiffPro[sd];};

  //void Set2pCorrelatorHist(TH1D* const s2pHist) {this->f2pCorrelatorHist = s2pHist;};
  //TH1D* Get2pCorrelatorHist() const {return this->f2pCorrelatorHist;};    

 private:
  AliFlowAnalysisWithMixedHarmonics(const AliFlowAnalysisWithMixedHarmonics& afawQc);
  AliFlowAnalysisWithMixedHarmonics& operator=(const AliFlowAnalysisWithMixedHarmonics& afawQc); 
  
  // 0.) Base:
  TList *fHistList; // base list to hold all output objects
  TString *fHistListName; // name of base list
  Int_t fHarmonic; // harmonic n in cos[n*(phi1+phi2-2phi3)] and cos[n*(psi1+psi2-2phi3)]
  TString *fAnalysisLabel; // analysis label 
  TProfile *fAnalysisSettings; // profile to hold analysis settings
  Int_t fNoOfMultipicityBins; // number of multiplicity bins
  Double_t fMultipicityBinWidth; // width of multiplicity bin
  Double_t fMinMultiplicity; // minimal multiplicity
  Bool_t fOppositeChargesPOI; // two POIs, psi1 and psi2, in correlator <<cos[psi1+psi2-2phi3)]>> will be taken with opposite charges
  Bool_t fEvaluateDifferential3pCorrelator; // evaluate <<cos[psi1+psi2-2phi3)]>>, where psi1 and psi2 are two POIs 
  Bool_t fCorrectForDetectorEffects; // correct 3-p correlator for detector effects
  Bool_t fPrintOnTheScreen; // print or not the final results on the screen
  Bool_t fCalculateVsM; // calculate correlators vs multiplicity
  Bool_t fShowBinLabelsVsM; // in histograms holding results vs multiplicity show bin labels in the format M_lowEdge \leq M < M_upperEdge
  
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
  TProfile *fRePEBE[2]; // real part of p_n vs [(p1+p2)/2,|p1-p2|]
  TProfile *fImPEBE[2]; // imaginary part of p_n vs [(p1+p2)/2,|p1-p2|]
  TProfile *fOverlapEBE[2][2]; // cos[n(psi-phi)] vs [(p1+p2)/2,|p1-p2|], where phi stands for 1st/2nd POI which is also RP 
  TProfile *fReEtaEBE[2]; // real part of p_n vs [(eta1+eta2)/2,|eta1-eta2|]
  TProfile *fImEtaEBE[2]; // imaginary part of p_n vs [(eta1+eta2)/2,|eta1-eta2|]
  TProfile *fOverlapEBE2[2][2]; // cos[n(psi-phi)] vs [(eta1+eta2)/2,|eta1-eta2|], where phi stands for 1st/2nd POI which is also RP 
  
  // 4.) Profiles:
  TList *fProfileList; // list holding all all-event profiles 
  TProfile *f3pCorrelatorPro; // 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> (not corrected for detector effects)
  TProfile *f5pCorrelatorPro; // 5-p correlator <<cos[n*(2.*phi1+2.*phi2+2.*phi3-3.*phi4-3.*phi5)]>> (not corrected for detector effects)
  TProfile *fNonIsotropicTermsPro; // non-isotropic terms in the decomposition of 3-p correlator <<cos[n(phi1+phi2-2phi3)]>>
  TProfile *f3pCorrelatorVsMPro; // 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> vs multiplicity
  TProfile *f3pPOICorrelatorVsM; // 3-p correlator <<cos[n(psi1+psi2-2phi3)]>> vs multiplicity
  TProfile2D *fNonIsotropicTermsVsMPro; // non-isotropic terms in the decomposition of <cos[n(phi1+phi2-2phi3))]> vs multiplicity
  TProfile *f3pCorrelatorVsPtSumDiffPro[2]; // differential 3-p correlator <<cos[psi1+psi2-2phi3)]>> vs [(p1+p2)/2,|p1-p2|]
  TProfile *f3pCorrelatorVsEtaSumDiffPro[2]; // differential 3-p correlator <<cos[psi1+psi2-2phi3)]>> vs [(eta1+eta2)/2,|eta1-eta2|]

  //2p correlators vs |Pt1 - Pt2|
  TProfile *f2pCorrelatorCosPsiDiffPtDiff; // <<cos[n(psi1-psi2)] vs pt diff 
  TProfile *f2pCorrelatorCosPsiSumPtDiff; // <<cos[n(psi1+psi2)]  vs pt diff 
  TProfile *f2pCorrelatorSinPsiDiffPtDiff; // <<sin[n(psi1-psi2)]  vs pt diff 
  TProfile *f2pCorrelatorSinPsiSumPtDiff; // <<sin[n(psi1+psi2)]  vs pt diff 

  //2p correlators vs (Pt1 + Pt2)/2
  TProfile *f2pCorrelatorCosPsiDiffPtSum; // <<cos[n(psi1-psi2)] vs pt sum 
  TProfile *f2pCorrelatorCosPsiSumPtSum; // <<cos[n(psi1+psi2)]  vs pt sum 
  TProfile *f2pCorrelatorSinPsiDiffPtSum; // <<sin[n(psi1-psi2)]  vs pt sum 
  TProfile *f2pCorrelatorSinPsiSumPtSum; // <<sin[n(psi1+psi2)]  vs pt sum 

  //2p correlators vs |eta1 - eta2|
  TProfile *f2pCorrelatorCosPsiDiffEtaDiff; // <<cos[n(psi1-psi2)] vs eta diff 
  TProfile *f2pCorrelatorCosPsiSumEtaDiff; // <<cos[n(psi1+psi2)]  vs eta diff 
  TProfile *f2pCorrelatorSinPsiDiffEtaDiff; // <<sin[n(psi1-psi2)]  vs eta diff 
  TProfile *f2pCorrelatorSinPsiSumEtaDiff; // <<sin[n(psi1+psi2)]  vs eta diff 

  //2p correlators vs (eta1 + eta2)/2
  TProfile *f2pCorrelatorCosPsiDiffEtaSum; // <<cos[n(psi1-psi2)] vs eta sum 
  TProfile *f2pCorrelatorCosPsiSumEtaSum; // <<cos[n(psi1+psi2)]  vs eta sum 
  TProfile *f2pCorrelatorSinPsiDiffEtaSum; // <<sin[n(psi1-psi2)]  vs eta sum 
  TProfile *f2pCorrelatorSinPsiSumEtaSum; // <<sin[n(psi1+psi2)]  vs eta sum 

  // 5.) Final results:
  TList *fResultsList; // list holding objects with final results 
  TH1D *f3pCorrelatorHist; // 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> corrected for detector effects
  TH1D *fDetectorBiasHist; // bias coming from detector inefficiencies to 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> (in %)
  TH1D *f3pCorrelatorVsMHist; // 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> vs multiplicity corrected for detector effects
  TH1D *fDetectorBiasVsMHist; // bias coming from detector inefficiencies to 3-p correlator <<cos[n(phi1+phi2-2phi3)]>> (in %) versus multiplicity
  //TH1D *f2pCorrelatorHist;//<<cos[(psi1-psi2)]>>

  ClassDef(AliFlowAnalysisWithMixedHarmonics, 0);

};

//================================================================================================================

#endif





