/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/*************************************************************** 
 * Only in this class nested loops are used for flow analysis. *
 * Nested loops are used to evaluate:                          *
 *                                                             *  
 *  a) Distribution of relative angle difference (phi1-phi2);  *
 *  b) Cross-check the results for mixed harmonics.            *
 *                                                             *
 *       Author: Ante Bilandzic (abilandzic@gmail.com)         *
 ***************************************************************/ 

#ifndef ALIFLOWANALYSISNESTEDLOOPS_H
#define ALIFLOWANALYSISNESTEDLOOPS_H

#include "AliFlowCommonConstants.h" // needed as include

class TList;
class TDirectoryFile;
class TH1F;
class TH1D;
class TProfile;

class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

//================================================================================================================

class AliFlowAnalysisWithNestedLoops
{
 public:
  AliFlowAnalysisWithNestedLoops();
  virtual ~AliFlowAnalysisWithNestedLoops(); 
  // 0.) Methods called in the constructor:
    virtual void InitializeArraysForMH();
  // 1.) Method Init() and methods called within Init():
  virtual void Init();
    virtual void CrossCheckSettings();
    virtual void AccessConstants();
    virtual void BookAndNestAllLists();
    virtual void BookAndFillProfileHoldingSettings();
    virtual void BookCommonHistograms();
    virtual void BookEverythingForRAD(); // RAD = relative angle distribution phi1-phi2
    virtual void BookEverythingForMH(); // MH = Mixed Harmonics
    virtual void BookAndFillWeightsHistograms();
  // 2.) Method Make() and methods called within Make():
  virtual void Make(AliFlowEventSimple *anEvent);
    virtual void CheckPointersUsedInMake();
    virtual void EvaluateNestedLoopsForRAD(AliFlowEventSimple *anEvent);
    virtual void EvaluateNestedLoopsForMH(AliFlowEventSimple *anEvent);
  // 3.) Method Finish() and methods called within Finish():
  virtual void Finish();  
    virtual void CheckPointersUsedInFinish(); 
    virtual void AccessSettings();  
    virtual void PrintOnTheScreen();     
  // 4.) Method GetOutputHistograms and method called within it:
  virtual void GetOutputHistograms(TList *outputListHistos);
    virtual void GetPointersForBaseHistograms();
    virtual void GetPointersForCommonHistograms();
    virtual void GetPointersForRAD();
    virtual void GetPointersForMH();
  // 5.) Other methods:   
  virtual void WriteHistograms(TString outputFileName);
  virtual void WriteHistograms(TDirectoryFile *outputFileName);  
  virtual void CheckPointersForRAD(TString where);
  virtual void CheckPointersForMH(TString where);
  // 6.) Setters and getters:
  void SetHistList(TList* const hl) {this->fHistList = hl;}
  TList* GetHistList() const {return this->fHistList;}  
  void SetHistListName(const char *hln) {this->fHistListName->Append(*hln);}; 
  TString *GetHistListName() const {return this->fHistListName;};
  void SetAnalysisLabel(const char *al) {this->fAnalysisLabel->Append(*al);}; 
  TString *GetAnalysisLabel() const {return this->fAnalysisLabel;};
  void SetAnalysisSettings(TProfile* const as) {this->fAnalysisSettings = as;};
  TProfile* GetAnalysisSettings() const {return this->fAnalysisSettings;};
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
  void SetListRAD(TList* const lRAD) {this->fListRAD = lRAD;}
  TList* GetListRAD() const {return this->fListRAD;}  
  void SetEvaluateNestedLoopsForRAD(Bool_t const enlfRAD) {this->fEvaluateNestedLoopsForRAD = enlfRAD;};
  Bool_t GetEvaluateNestedLoopsForRAD() const {return this->fEvaluateNestedLoopsForRAD;};
  void SetRelativeAngleDistribution(TH1D* const rad) {this->fRelativeAngleDistribution = rad;};
  TH1D* GetRelativeAngleDistribution() const {return this->fRelativeAngleDistribution;}; 
  // QC:
  // ...
  // MH:
  void SetListMH(TList* const lMH) {this->fListMH = lMH;}
  TList* GetListMH() const {return this->fListMH;}  
  void SetEvaluateNestedLoopsForMH(Bool_t const enlfMH) {this->fEvaluateNestedLoopsForMH = enlfMH;};
  Bool_t GetEvaluateNestedLoopsForMH() const {return this->fEvaluateNestedLoopsForMH;};
  void SetCorrelatorIntegerMH(Int_t const ciMH) {this->fCorrelatorIntegerMH = ciMH;};
  Int_t GetCorrelatorIntegerMH() const {return this->fCorrelatorIntegerMH;}; 
  void Set3pCorrelatorVsPtSumDiffDirectPro(TProfile* const s3pcvpsdd, Int_t const sd) {this->f3pCorrelatorVsPtSumDiffDirectPro[sd] = s3pcvpsdd;};
  TProfile* Get3pCorrelatorVsPtSumDiffDirectPro(Int_t sd) const {return this->f3pCorrelatorVsPtSumDiffDirectPro[sd];};
  void SetCrossCheckInPtSumBinNo(Int_t const ccipsbn) {this->fCrossCheckInPtSumBinNo = ccipsbn;};
  Int_t GetCrossCheckInPtSumBinNo() const {return this->fCrossCheckInPtSumBinNo;}; 
  void SetCrossCheckInPtDiffBinNo(Int_t const ccipdbn) {this->fCrossCheckInPtDiffBinNo = ccipdbn;};
  Int_t GetCrossCheckInPtDiffBinNo() const {return this->fCrossCheckInPtDiffBinNo;};  
  
 private:
  AliFlowAnalysisWithNestedLoops(const AliFlowAnalysisWithNestedLoops& afawQc);
  AliFlowAnalysisWithNestedLoops& operator=(const AliFlowAnalysisWithNestedLoops& afawQc); 
  // 0.) Base:
  TList *fHistList; // base list to hold all output objects
  TString *fHistListName; // name of base list
  TString *fAnalysisLabel; // analysis label 
  TProfile *fAnalysisSettings; // profile to hold analysis settings
  Bool_t fPrintOnTheScreen; // print or not on the screen
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
  // 3.) Relative angle distribution (RAD):
  TList *fListRAD; // list holding objects for calculation of relative angle distribution phi1-phi2 
  Bool_t fEvaluateNestedLoopsForRAD; // evaluate nested loops for relative angle distribution
  TH1D *fRelativeAngleDistribution; // distribution of phi1-phi2 for all distinct pairs of particles
  // 4.) Debugging and cross-checking QC:
  // ...
  // 5.) Debugging and cross-checking MH:
  TList *fListMH; // list holding objects relevant for debugging and cross-checking of MH class
  Bool_t fEvaluateNestedLoopsForMH; // evaluate nested loops for mixed harmonics
  Int_t fCorrelatorIntegerMH; // integer n in cos[n(2phi1-psi2-psi3)]
  TProfile *f3pCorrelatorVsPtSumDiffDirectPro[2]; // differential 3-p correlator cos[n(2phi1-psi2-psi3)] vs [(p1+p2)/2,|p1-p2|]
  Int_t fCrossCheckInPtSumBinNo; // print on the screen result of cos[n(2phi1-psi2-psi3)] vs (p1+p2)/2 in this bin number 
  Int_t fCrossCheckInPtDiffBinNo; // print on the screen result of cos[n(2phi1-psi2-psi3)] vs |p1-p2| in this bin number 
  
  ClassDef(AliFlowAnalysisWithNestedLoops, 0);
};

//================================================================================================================

#endif





