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
class AliFlowTrackSimple;
class AliFlowVector;

class AliFlowCommonHist;
class AliFlowCommonHistResults;

//================================================================================================================

class AliFlowAnalysisWithQCumulants{
 public:
  AliFlowAnalysisWithQCumulants();
  virtual ~AliFlowAnalysisWithQCumulants(); 
  
  virtual void Init();
  virtual void Make(AliFlowEventSimple* anEvent);
  
  virtual void CalculateCorrelationsForIntegratedFlow(); // everything cross-checked (2-8)
  virtual void CalculateWeightedCorrelationsForIntegratedFlow();
  virtual void CalculateCorrelationsForDifferentialFlow(TString type="POI");
  virtual void CalculateWeightedCorrelationsForDifferentialFlow(TString type="POI");
  
  virtual void EvaluateNestedLoopsForIntegratedFlow(AliFlowEventSimple* anEvent); 
  virtual void EvaluateNestedLoopsForDifferentialFlow(AliFlowEventSimple* anEvent); 
  
  virtual void Finish();
  
  TProfile* MakePtProjection(TProfile2D *profilePtEta) const;
  TProfile* MakeEtaProjection(TProfile2D *profilePtEta) const;
  
  virtual void CalculateFinalResultsForNoNameIntegratedFlow(Bool_t useWeights=kFALSE);
  virtual void CalculateFinalResultsForRPandPOIIntegratedFlow(Bool_t useWeights, TString type);
  virtual void CalculateFinalResultsForDifferentialFlow(TH2D *flowPtEta, TH1D *flowPt, TH1D *flowEta, 
                                                        TProfile2D *profile2ndPtEta, TProfile2D *profile4thPtEta = NULL, 
                                                        TProfile2D *profile6thPtEta = NULL, TProfile2D *profile8thPtEta = NULL);
  
  virtual void PrintFinalResultsForIntegratedFlow(Bool_t useWeights=kTRUE, TString type="NONAME");
    
  virtual void CompareDirectAndQCorrelationsForIntegratedFlow(Bool_t useWeights);
  virtual void CompareDirectAndQCorrelationsForDifferentialFlow(Bool_t useWeights);

  virtual void WriteHistograms(TString outputFileName);
  
  virtual void TempDeleteMe();
 
//----------------------------------------------------------------------------------------------------------------
//                                            setters and getters                                                 
//----------------------------------------------------------------------------------------------------------------
  TList* GetHistList() const {return this->fHistList;} 
  
  void SetWeightsList(TList* wlist) {this->fWeightsList = wlist;}
  TList* GetWeightsList() const {return this->fWeightsList;}  
  
  void SetResultsList(TList* rlist) {this->fResultsList = rlist;}
  TList* GetResultsList() const {return this->fResultsList;}  
 
  void SetIntFlowResults(TH1D* const ifr) {this->fIntFlowResultsQC = ifr;};
  TH1D* GetIntFlowResults() const {return this->fIntFlowResultsQC;};
  
  void SetIntFlowResultsW(TH1D* const ifrw) {this->fIntFlowResultsQCW = ifrw;};
  TH1D* GetIntFlowResultsW() const {return this->fIntFlowResultsQCW;};
  
  void SetIntFlowResultsPOI(TH1D* const ifrp) {this->fIntFlowResultsPOIQC = ifrp;};
  TH1D* GetIntFlowResultsPOI() const {return this->fIntFlowResultsPOIQC;};
  
  void SetIntFlowResultsPOIW(TH1D* const ifrpw) {this->fIntFlowResultsPOIQCW = ifrpw;};
  TH1D* GetIntFlowResultsPOIW() const {return this->fIntFlowResultsPOIQCW;};

  void SetIntFlowResultsRP(TH1D* const ifrr) {this->fIntFlowResultsRPQC = ifrr;};
  TH1D* GetIntFlowResultsRP() const {return this->fIntFlowResultsRPQC;};
  
  void SetIntFlowResultsRPW(TH1D* const ifrrw) {this->fIntFlowResultsRPQCW = ifrrw;};
  TH1D* GetIntFlowResultsRPW() const {return this->fIntFlowResultsRPQCW;};





  void SetDiffFlowResults2nd(TH1D* const diff2nd) {this->fDiffFlowResults2ndOrderQC = diff2nd;};
  TH1D* GetDiffFlowResults2nd() const {return this->fDiffFlowResults2ndOrderQC;};
  
  void SetDiffFlowResults4th(TH1D* const diff4th) {this->fDiffFlowResults4thOrderQC = diff4th;};
  TH1D* GetDiffFlowResults4th() const {return this->fDiffFlowResults4thOrderQC;};
  
  void SetCovariances(TH1D* const cov) {this->fCovariances = cov;};
  TH1D* GetCovariances() const {return this->fCovariances;};
  
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
  
  void SetAverageMultiplicity(TProfile* const am) {this->fAvMultIntFlowQC = am;};
  TProfile* GetAverageMultiplicity() const {return this->fAvMultIntFlowQC;};
  
  void SetQvectorForEachEventX(TProfile* const qvfeex) {this->fQvectorForEachEventX = qvfeex;};
  TProfile* GetQvectorForEachEventX() const {return this->fQvectorForEachEventX;};

  void SetQvectorForEachEventY(TProfile* const qvfeey) {this->fQvectorForEachEventY = qvfeey;};
  TProfile* GetQvectorForEachEventY() const {return this->fQvectorForEachEventY;};
        
  void SetQCorrelations(TProfile* const QCorr) {this->fQCorrelations = QCorr;};
  TProfile* GetQCorrelations() const {return this->fQCorrelations;};
  
  void SetQCorrelationsW(TProfile* const QCorrW) {this->fQCorrelationsW = QCorrW;};
  TProfile* GetQCorrelationsW() const {return this->fQCorrelationsW;};
  
  void SetQProduct(TProfile* const qp) {this->fQProduct = qp;};
  TProfile* GetQProduct() const {return this->fQProduct;};
  
  void SetQVectorComponents(TProfile* const qvc) {this->fQvectorComponents = qvc;};
  TProfile* GetQVectorComponents() const {return this->fQvectorComponents;};
  
  void SetTwo1n1nPerPtBinRP(TProfile* const pb2PerPtBin1n1nRP) {this->f2PerPtBin1n1nRP = pb2PerPtBin1n1nRP;};
  TProfile* GetTwo1n1nPerPtBinRP() const {return this->f2PerPtBin1n1nRP;};
  
  void SetFour1n1n1n1nPerPtBinRP(TProfile* const pb4PerPtBin1n1n1n1nRP) {this->f4PerPtBin1n1n1n1nRP = pb4PerPtBin1n1n1n1nRP;};
  TProfile* GetFour1n1n1n1nPerPtBinRP() const {return this->f4PerPtBin1n1n1n1nRP;}; 

  void SetTwo1n1nPerEtaBinRP(TProfile* const pb2PerEtaBin1n1nRP) {this->f2PerEtaBin1n1nRP = pb2PerEtaBin1n1nRP;};
  TProfile* GetTwo1n1nPerEtaBinRP() const {return this->f2PerEtaBin1n1nRP;};
  
  void SetFour1n1n1n1nPerEtaBinRP(TProfile* const pb4PerEtaBin1n1n1n1nRP) {this->f4PerEtaBin1n1n1n1nRP = pb4PerEtaBin1n1n1n1nRP;};
  TProfile* GetFour1n1n1n1nPerEtaBinRP() const {return this->f4PerEtaBin1n1n1n1nRP;}; 
  
  void SetTwo1n1nPerPtBinPOI(TProfile* const pb2PerPtBin1n1nPOI) {this->f2PerPtBin1n1nPOI = pb2PerPtBin1n1nPOI;};
  TProfile* GetTwo1n1nPerPtBinPOI() const {return this->f2PerPtBin1n1nPOI;};
  
  void SetFour1n1n1n1nPerPtBinPOI(TProfile* const pb4PerPtBin1n1n1n1nPOI) {this->f4PerPtBin1n1n1n1nPOI = pb4PerPtBin1n1n1n1nPOI;};
  TProfile* GetFour1n1n1n1nPerPtBinPOI() const {return this->f4PerPtBin1n1n1n1nPOI;}; 
  
  void SetTwo1n1nPerEtaBinPOI(TProfile* const pb2PerEtaBin1n1nPOI) {this->f2PerEtaBin1n1nPOI = pb2PerEtaBin1n1nPOI;};
  TProfile* GetTwo1n1nPerEtaBinPOI() const {return this->f2PerEtaBin1n1nPOI;};
  
  void SetFour1n1n1n1nPerEtaBinPOI(TProfile* const pb4PerEtaBin1n1n1n1nPOI) {this->f4PerEtaBin1n1n1n1nPOI = pb4PerEtaBin1n1n1n1nPOI;};
  TProfile* GetFour1n1n1n1nPerEtaBinPOI() const {return this->f4PerEtaBin1n1n1n1nPOI;}; 
  
  void SetTwo1n1nWPerPtBinPOI(TProfile* const pb2WPerPtBin1n1nPOI) {this->f2WPerPtBin1n1nPOI = pb2WPerPtBin1n1nPOI;};
  TProfile* GetTwo1n1nWPerPtBinPOI() const {return this->f2WPerPtBin1n1nPOI;};
  
  void SetFour1n1n1n1nWPerPtBinPOI(TProfile* const pb4WPerPtBin1n1n1n1nPOI) {this->f4WPerPtBin1n1n1n1nPOI = pb4WPerPtBin1n1n1n1nPOI;};
  TProfile* GetFour1n1n1n1nWPerPtBinPOI() const {return this->f4WPerPtBin1n1n1n1nPOI;}
  
  void SetTwo1n1nWPerEtaBinPOI(TProfile* const pb2WPerEtaBin1n1nPOI) {this->f2WPerEtaBin1n1nPOI = pb2WPerEtaBin1n1nPOI;};
  TProfile* GetTwo1n1nWPerEtaBinPOI() const {return this->f2WPerEtaBin1n1nPOI;};
  
  void SetFour1n1n1n1nWPerEtaBinPOI(TProfile* const pb4WPerEtaBin1n1n1n1nPOI) {this->f4WPerEtaBin1n1n1n1nPOI = pb4WPerEtaBin1n1n1n1nPOI;};
  TProfile* GetFour1n1n1n1nWPerEtaBinPOI() const {return this->f4WPerEtaBin1n1n1n1nPOI;}
  
  void SetTwo1n1nWPerPtBinRP(TProfile* const pb2WPerPtBin1n1nRP) {this->f2WPerPtBin1n1nRP = pb2WPerPtBin1n1nRP;};
  TProfile* GetTwo1n1nWPerPtBinRP() const {return this->f2WPerPtBin1n1nRP;};
  
  void SetFour1n1n1n1nWPerPtBinRP(TProfile* const pb4WPerPtBin1n1n1n1nRP) {this->f4WPerPtBin1n1n1n1nRP = pb4WPerPtBin1n1n1n1nRP;};
  TProfile* GetFour1n1n1n1nWPerPtBinRP() const {return this->f4WPerPtBin1n1n1n1nRP;}
  
  void SetTwo1n1nWPerEtaBinRP(TProfile* const pb2WPerEtaBin1n1nRP) {this->f2WPerEtaBin1n1nRP = pb2WPerEtaBin1n1nRP;};
  TProfile* GetTwo1n1nWPerEtaBinRP() const {return this->f2WPerEtaBin1n1nRP;};
  
  void SetFour1n1n1n1nWPerEtaBinRP(TProfile* const pb4WPerEtaBin1n1n1n1nRP) {this->f4WPerEtaBin1n1n1n1nRP = pb4WPerEtaBin1n1n1n1nRP;};
  TProfile* GetFour1n1n1n1nWPerEtaBinRP() const {return this->f4WPerEtaBin1n1n1n1nRP;}
  
  void SetDirectCorrelations(TProfile* const dc) {this->fDirectCorrelations = dc;};
  TProfile* GetDirectCorrelations() const {return this->fDirectCorrelations;};
  
  void SetDirectCorrelationsW(TProfile* const dcw) {this->fDirectCorrelationsW = dcw;};
  TProfile* GetDirectCorrelationsW() const {return this->fDirectCorrelationsW;};
  
  void SetDirectCorrelationsDiffFlow(TProfile* const dcdf) {this->fDirectCorrelationsDiffFlow = dcdf;};
  TProfile* GetDirectCorrelationsDiffFlow() const {return this->fDirectCorrelationsDiffFlow;};
  
  void SetDirectCorrelationsDiffFlowW(TProfile* const dcdfw) {this->fDirectCorrelationsDiffFlowW = dcdfw;};
  TProfile* GetDirectCorrelationsDiffFlowW() const {return this->fDirectCorrelationsDiffFlowW;};
  
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
  
  void SetUseWeights(Bool_t const uw) {this->fUseWeights = uw;};
  Bool_t GetUseWeights() const {return this->fUseWeights;};
  
  void SetUseWeightsBits(TBits* const uwb) {this->fUseWeightsBits = uwb;};
  TBits* GetUseWeightsBits() const {return this->fUseWeightsBits;};
  
  
  
  
  
  
  
  // .................................................................................................
  // non-weighted correlations for differential flow of POIs:
  void Set2pPtEtaPOI(TProfile2D* const tppep) {this->f2pPtEtaPOI = tppep;};
  TProfile2D* Get2pPtEtaPOI() const {return this->f2pPtEtaPOI;};
  void Set4pPtEtaPOI(TProfile2D* const fppep) {this->f4pPtEtaPOI = fppep;};
  TProfile2D* Get4pPtEtaPOI() const {return this->f4pPtEtaPOI;};
  void Set6pPtEtaPOI(TProfile2D* const sppep) {this->f6pPtEtaPOI = sppep;};
  TProfile2D* Get6pPtEtaPOI() const {return this->f6pPtEtaPOI;};
  void Set8pPtEtaPOI(TProfile2D* const eppep) {this->f8pPtEtaPOI = eppep;};
  TProfile2D* Get8pPtEtaPOI() const {return this->f8pPtEtaPOI;};
  
  // non-weighted final results for differential flow of POIs:
  // 3D (pt,eta):
  void Setvn2ndPtEtaPOI(TH2D* const v2pep) {this->fvn2ndPtEtaPOI = v2pep;};
  TH2D* Getvn2ndPtEtaPOI() const {return this->fvn2ndPtEtaPOI;};
  void Setvn4thPtEtaPOI(TH2D* const v4pep) {this->fvn4thPtEtaPOI = v4pep;};
  TH2D* Getvn4thPtEtaPOI() const {return this->fvn4thPtEtaPOI;};
  void Setvn6thPtEtaPOI(TH2D* const v6pep) {this->fvn6thPtEtaPOI = v6pep;};
  TH2D* Getvn6thPtEtaPOI() const {return this->fvn6thPtEtaPOI;};
  void Setvn8thPtEtaPOI(TH2D* const v8pep) {this->fvn8thPtEtaPOI = v8pep;};
  TH2D* Getvn8thPtEtaPOI() const {return this->fvn8thPtEtaPOI;};
  // 2D (pt):
  void Setvn2ndPtPOI(TH1D* const v2pp) {this->fvn2ndPtPOI = v2pp;};
  TH1D* Getvn2ndPtPOI() const {return this->fvn2ndPtPOI;};
  void Setvn4thPtPOI(TH1D* const v4pp) {this->fvn4thPtPOI = v4pp;};
  TH1D* Getvn4thPtPOI() const {return this->fvn4thPtPOI;};
  void Setvn6thPtPOI(TH1D* const v6pp) {this->fvn6thPtPOI = v6pp;};
  TH1D* Getvn6thPtPOI() const {return this->fvn6thPtPOI;};
  void Setvn8thPtPOI(TH1D* const v8pp) {this->fvn8thPtPOI = v8pp;};
  TH1D* Getvn8thPtPOI() const {return this->fvn8thPtPOI;};
  // 2D (eta):
  void Setvn2ndEtaPOI(TH1D* const v2ep) {this->fvn2ndEtaPOI = v2ep;};
  TH1D* Getvn2ndEtaPOI() const {return this->fvn2ndEtaPOI;};
  void Setvn4thEtaPOI(TH1D* const v4ep) {this->fvn4thEtaPOI = v4ep;};
  TH1D* Getvn4thEtaPOI() const {return this->fvn4thEtaPOI;};
  void Setvn6thEtaPOI(TH1D* const v6ep) {this->fvn6thEtaPOI = v6ep;};
  TH1D* Getvn6thEtaPOI() const {return this->fvn6thEtaPOI;};
  void Setvn8thEtaPOI(TH1D* const v8ep) {this->fvn8thEtaPOI = v8ep;};
  TH1D* Getvn8thEtaPOI() const {return this->fvn8thEtaPOI;};
  
  // weighted correlations for differential flow of POIs:
  void Set2pPtEtaPOIW(TProfile2D* const tppepw) {this->f2pPtEtaPOIW = tppepw;};
  TProfile2D* Get2pPtEtaPOIW() const {return this->f2pPtEtaPOIW;};
  void Set4pPtEtaPOIW(TProfile2D* const fppepw) {this->f4pPtEtaPOIW = fppepw;};
  TProfile2D* Get4pPtEtaPOIW() const {return this->f4pPtEtaPOIW;};
  void Set6pPtEtaPOIW(TProfile2D* const sppepw) {this->f6pPtEtaPOIW = sppepw;};
  TProfile2D* Get6pPtEtaPOIW() const {return this->f6pPtEtaPOIW;};
  void Set8pPtEtaPOIW(TProfile2D* const eppepw) {this->f8pPtEtaPOIW = eppepw;};
  TProfile2D* Get8pPtEtaPOIW() const {return this->f8pPtEtaPOIW;};
  
  // weighted final results for differential flow of POIs:
  // 3D (pt,eta):
  void Setvn2ndPtEtaPOIW(TH2D* const v2pepw) {this->fvn2ndPtEtaPOIW = v2pepw;};
  TH2D* Getvn2ndPtEtaPOIW() const {return this->fvn2ndPtEtaPOIW;};
  void Setvn4thPtEtaPOIW(TH2D* const v4pepw) {this->fvn4thPtEtaPOIW = v4pepw;};
  TH2D* Getvn4thPtEtaPOIW() const {return this->fvn4thPtEtaPOIW;};
  void Setvn6thPtEtaPOIW(TH2D* const v6pepw) {this->fvn6thPtEtaPOIW = v6pepw;};
  TH2D* Getvn6thPtEtaPOIW() const {return this->fvn6thPtEtaPOIW;};
  void Setvn8thPtEtaPOIW(TH2D* const v8pepw) {this->fvn8thPtEtaPOIW = v8pepw;};
  TH2D* Getvn8thPtEtaPOIW() const {return this->fvn8thPtEtaPOIW;};
  // 2D (pt):
  void Setvn2ndPtPOIW(TH1D* const v2ppw) {this->fvn2ndPtPOIW = v2ppw;};
  TH1D* Getvn2ndPtPOIW() const {return this->fvn2ndPtPOIW;};
  void Setvn4thPtPOIW(TH1D* const v4ppw) {this->fvn4thPtPOIW = v4ppw;};
  TH1D* Getvn4thPtPOIW() const {return this->fvn4thPtPOIW;};
  void Setvn6thPtPOIW(TH1D* const v6ppw) {this->fvn6thPtPOIW = v6ppw;};
  TH1D* Getvn6thPtPOIW() const {return this->fvn6thPtPOIW;};
  void Setvn8thPtPOIW(TH1D* const v8ppw) {this->fvn8thPtPOIW = v8ppw;};
  TH1D* Getvn8thPtPOIW() const {return this->fvn8thPtPOIW;};
  // 2D (eta):
  void Setvn2ndEtaPOIW(TH1D* const v2epw) {this->fvn2ndEtaPOIW = v2epw;};
  TH1D* Getvn2ndEtaPOIW() const {return this->fvn2ndEtaPOIW;};
  void Setvn4thEtaPOIW(TH1D* const v4epw) {this->fvn4thEtaPOIW = v4epw;};
  TH1D* Getvn4thEtaPOIW() const {return this->fvn4thEtaPOIW;};
  void Setvn6thEtaPOIW(TH1D* const v6epw) {this->fvn6thEtaPOIW = v6epw;};
  TH1D* Getvn6thEtaPOIW() const {return this->fvn6thEtaPOIW;};
  void Setvn8thEtaPOIW(TH1D* const v8epw) {this->fvn8thEtaPOIW = v8epw;};
  TH1D* Getvn8thEtaPOIW() const {return this->fvn8thEtaPOIW;};
       
  // non-weighted correlations for differential flow of RPs:
  void Set2pPtEtaRP(TProfile2D* const tpper) {this->f2pPtEtaRP = tpper;};
  TProfile2D* Get2pPtEtaRP() const {return this->f2pPtEtaRP;};
  void Set4pPtEtaRP(TProfile2D* const fpper) {this->f4pPtEtaRP = fpper;};
  TProfile2D* Get4pPtEtaRP() const {return this->f4pPtEtaRP;};
  void Set6pPtEtaRP(TProfile2D* const spper) {this->f6pPtEtaRP = spper;};
  TProfile2D* Get6pPtEtaRP() const {return this->f6pPtEtaRP;};
  void Set8pPtEtaRP(TProfile2D* const epper) {this->f8pPtEtaRP = epper;};
  TProfile2D* Get8pPtEtaRP() const {return this->f8pPtEtaRP;};
  
  // non-weighted final results for differential flow of RPs:
  // 3D (pt,eta):
  void Setvn2ndPtEtaRP(TH2D* const v2per) {this->fvn2ndPtEtaRP = v2per;};
  TH2D* Getvn2ndPtEtaRP() const {return this->fvn2ndPtEtaRP;};
  void Setvn4thPtEtaRP(TH2D* const v4per) {this->fvn4thPtEtaRP = v4per;};
  TH2D* Getvn4thPtEtaRP() const {return this->fvn4thPtEtaRP;};
  void Setvn6thPtEtaRP(TH2D* const v6per) {this->fvn6thPtEtaRP = v6per;};
  TH2D* Getvn6thPtEtaRP() const {return this->fvn6thPtEtaRP;};
  void Setvn8thPtEtaRP(TH2D* const v8per) {this->fvn4thPtEtaRP = v8per;};
  TH2D* Getvn8thPtEtaRP() const {return this->fvn8thPtEtaRP;};
  // 2D (pt):
  void Setvn2ndPtRP(TH1D* const v2pp) {this->fvn2ndPtRP = v2pp;};
  TH1D* Getvn2ndPtRP() const {return this->fvn2ndPtRP;};
  void Setvn4thPtRP(TH1D* const v4pp) {this->fvn4thPtRP = v4pp;};
  TH1D* Getvn4thPtRP() const {return this->fvn4thPtRP;};
  void Setvn6thPtRP(TH1D* const v6pp) {this->fvn6thPtRP = v6pp;};
  TH1D* Getvn6thPtRP() const {return this->fvn6thPtRP;};
  void Setvn8thPtRP(TH1D* const v8pp) {this->fvn8thPtRP = v8pp;};
  TH1D* Getvn8thPtRP() const {return this->fvn8thPtRP;};
  // 2D (eta):
  void Setvn2ndEtaRP(TH1D* const v2ep) {this->fvn2ndEtaRP = v2ep;};
  TH1D* Getvn2ndEtaRP() const {return this->fvn2ndEtaRP;};
  void Setvn4thEtaRP(TH1D* const v4ep) {this->fvn4thEtaRP = v4ep;};
  TH1D* Getvn4thEtaRP() const {return this->fvn4thEtaRP;};
  void Setvn6thEtaRP(TH1D* const v6ep) {this->fvn6thEtaRP = v6ep;};
  TH1D* Getvn6thEtaRP() const {return this->fvn6thEtaRP;};
  void Setvn8thEtaRP(TH1D* const v8ep) {this->fvn8thEtaRP = v8ep;};
  TH1D* Getvn8thEtaRP() const {return this->fvn8thEtaRP;};
  
  // weighted correlations for differential flow of RPs:
  void Set2pPtEtaRPW(TProfile2D* const tpperw) {this->f2pPtEtaRPW = tpperw;};
  TProfile2D* Get2pPtEtaRPW() const {return this->f2pPtEtaRPW;};
  void Set4pPtEtaRPW(TProfile2D* const fpperw) {this->f4pPtEtaRPW = fpperw;};
  TProfile2D* Get4pPtEtaRPW() const {return this->f4pPtEtaRPW;};
  void Set6pPtEtaRPW(TProfile2D* const spperw) {this->f6pPtEtaRPW = spperw;};
  TProfile2D* Get6pPtEtaRPW() const {return this->f6pPtEtaRPW;};
  void Set8pPtEtaRPW(TProfile2D* const epperw) {this->f8pPtEtaRPW = epperw;};
  TProfile2D* Get8pPtEtaRPW() const {return this->f8pPtEtaRPW;};
  
  // weighted final results for differential flow of RPs:
  // 3D (pt,eta):
  void Setvn2ndPtEtaRPW(TH2D* const v2perw) {this->fvn2ndPtEtaRPW = v2perw;};
  TH2D* Getvn2ndPtEtaRPW() const {return this->fvn2ndPtEtaRPW;}; 
  void Setvn4thPtEtaRPW(TH2D* const v4perw) {this->fvn4thPtEtaRPW = v4perw;};
  TH2D* Getvn4thPtEtaRPW() const {return this->fvn4thPtEtaRPW;};
  void Setvn6thPtEtaRPW(TH2D* const v6perw) {this->fvn6thPtEtaRPW = v6perw;};
  TH2D* Getvn6thPtEtaRPW() const {return this->fvn6thPtEtaRPW;};
  void Setvn8thPtEtaRPW(TH2D* const v8perw) {this->fvn4thPtEtaRPW = v8perw;};
  TH2D* Getvn8thPtEtaRPW() const {return this->fvn8thPtEtaRPW;};
  // 2D (pt):
  void Setvn2ndPtRPW(TH1D* const v2ppw) {this->fvn2ndPtRPW = v2ppw;};
  TH1D* Getvn2ndPtRPW() const {return this->fvn2ndPtRPW;};
  void Setvn4thPtRPW(TH1D* const v4ppw) {this->fvn4thPtRPW = v4ppw;};
  TH1D* Getvn4thPtRPW() const {return this->fvn4thPtRPW;};
  void Setvn6thPtRPW(TH1D* const v6ppw) {this->fvn6thPtRPW = v6ppw;};
  TH1D* Getvn6thPtRPW() const {return this->fvn6thPtRPW;};
  void Setvn8thPtRPW(TH1D* const v8ppw) {this->fvn8thPtRPW = v8ppw;};
  TH1D* Getvn8thPtRPW() const {return this->fvn8thPtRPW;};
  // 2D (eta):
  void Setvn2ndEtaRPW(TH1D* const v2epw) {this->fvn2ndEtaRPW = v2epw;};
  TH1D* Getvn2ndEtaRPW() const {return this->fvn2ndEtaRPW;};
  void Setvn4thEtaRPW(TH1D* const v4epw) {this->fvn4thEtaRPW = v4epw;};
  TH1D* Getvn4thEtaRPW() const {return this->fvn4thEtaRPW;};
  void Setvn6thEtaRPW(TH1D* const v6epw) {this->fvn6thEtaRPW = v6epw;};
  TH1D* Getvn6thEtaRPW() const {return this->fvn6thEtaRPW;};
  void Setvn8thEtaRPW(TH1D* const v8epw) {this->fvn8thEtaRPW = v8epw;};
  TH1D* Getvn8thEtaRPW() const {return this->fvn8thEtaRPW;};
  // .................................................................................................
  
  
  
  
//----------------------------------------------------------------------------------------------------------------
 
 private:
  AliFlowAnalysisWithQCumulants(const AliFlowAnalysisWithQCumulants& afawQc);
  AliFlowAnalysisWithQCumulants& operator=(const AliFlowAnalysisWithQCumulants& afawQc);
  
  AliFlowTrackSimple* fTrack;                           //track
  TList*              fHistList;                        //list to hold all output histograms
  TList*              fDiffFlowList;                    //list to hold all histograms and profiles needed for differential flow
  TList*              fWeightsList;                     //list to hold all histograms with weights
  TList*              fResultsList;                         // list to hold all histograms with results
    
  TProfile*           fAvMultIntFlowQC;                 //average selected multiplicity (for int. flow)
 
  TProfile*           fQvectorComponents;               //averages of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, ...)
            
  TH1D*               fDiffFlowResults2ndOrderQC;       //differential flow results from 2nd order Q-cumulant
  TH1D*               fDiffFlowResults4thOrderQC;       //differential flow results from 4th order Q-cumulant
  TH1D*               fCovariances;                     //final results for covariances: 1st bin: <2*4>-<2>*<4>, 2nd bin: <2*6>-<2>*<6>, ...
  
  TProfile*                  fQvectorForEachEventX;     //profile containing the x-components of Q-vectors from all events (to be removed)  
  TProfile*                  fQvectorForEachEventY;     //profile containing the y-components of Q-vectors from all events (to be removed)   
  TProfile*                  fQCorrelations;            //multi-particle correlations calculated from Q-vectors 
  TProfile*                  fQCorrelationsW;           //weighted multi-particle correlations calculated from Q-vectors 
  TProfile*                  fQProduct;                 //average of products: 1st bin: <2*4>, 2nd bin: <2*6>, ...
  
  TProfile*          fDirectCorrelations;               // multi-particle correlations calculated with nested loop needed for int. flow 
  TProfile*          fDirectCorrelationsW;              // multi-particle correlations calculated with nested loop needed for weighted int. flow
  TProfile*          fDirectCorrelationsDiffFlow;       // multi-particle correlations calculated with nested loop needed for diff. flow
  TProfile*          fDirectCorrelationsDiffFlowW;      // multi-particle correlations calculated with nested loop needed for weighted int. flow
  
  // POI (Particles Of Interest):
  // non-weighted correlations
  TProfile*                  f2PerPtBin1n1nPOI;         //<<2'>>_{n|n} per pt-bin
  TProfile*                  f4PerPtBin1n1n1n1nPOI;     //<<4'>>_{n,n|n,n} per pt-bin

  TProfile*                  f2PerEtaBin1n1nPOI;        //<<2'>>_{n|n} per eta-bin
  TProfile*                  f4PerEtaBin1n1n1n1nPOI;    //<<4'>>_{n,n|n,n} per eta-bin  
  // weighted correlations
  TProfile*                  f2WPerPtBin1n1nPOI;        //<<2'>>_{n|n} per eta-bin
  TProfile*                  f4WPerPtBin1n1n1n1nPOI;    //<<4'>>_{n,n|n,n} per eta-bin    
  
  TProfile*                  f2WPerEtaBin1n1nPOI;       //<<2'>>_{n|n} per eta-bin 
  TProfile*                  f4WPerEtaBin1n1n1n1nPOI;   //<<4'>>_{n,n|n,n} per eta-bin
  
  // RP (Reaction Plane particles)
  // non-weighted correlations
  TProfile*                  f2PerPtBin1n1nRP;          //<<2'>>_{n|n} per pt-bin
  TProfile*                  f4PerPtBin1n1n1n1nRP;      //<<4'>>_{n,n|n,n} per pt-bin

  TProfile*                  f2PerEtaBin1n1nRP;         //<<2'>>_{n|n} per eta-bin
  TProfile*                  f4PerEtaBin1n1n1n1nRP;     //<<4'>>_{n,n|n,n} per eta-bin  
  // weighted correlations
  TProfile*                  f2WPerPtBin1n1nRP;         //<<2'>>_{n|n} per eta-bin 
  TProfile*                  f4WPerPtBin1n1n1n1nRP;     //<<4'>>_{n,n|n,n} per eta-bin
  
  TProfile*                  f2WPerEtaBin1n1nRP;        //<<2'>>_{n|n} per eta-bin 
  TProfile*                  f4WPerEtaBin1n1n1n1nRP;    //<<4'>>_{n,n|n,n} per eta-bin
  
  AliFlowCommonHist*         fCommonHists2nd;           //common control histograms (taking into account only the events with 2 and more particles) 
  AliFlowCommonHist*         fCommonHists4th;           //common control histograms (taking into account only the events with 4 and more particles) 
  AliFlowCommonHist*         fCommonHists6th;           //common control histograms (taking into account only the events with 6 and more particles) 
  AliFlowCommonHist*         fCommonHists8th;           //common control histograms (taking into account only the events with 8 and more particles) 
  
  AliFlowCommonHistResults*  fCommonHistsResults2nd;    //final results for 2nd order int. and diff. flow stored in the common histograms 
  AliFlowCommonHistResults*  fCommonHistsResults4th;    //final results for 4th order int. and diff. flow stored in the common histograms 
  AliFlowCommonHistResults*  fCommonHistsResults6th;    //final results for 6th order int. and diff. flow stored in the common histograms
  AliFlowCommonHistResults*  fCommonHistsResults8th;    //final results for 8th order int. and diff. flow stored in the common histograms
      
  TH1D*                      f2pDistribution;           //distribution of <2>_{n|n}
  TH1D*                      f4pDistribution;           //distribution of <4>_{n,n|n,n}
  TH1D*                      f6pDistribution;           //distribution of <6>_{n,n,n|n,n,n} 
  TH1D*                      f8pDistribution;           //distribution of <8>_{n,n,n,n|n,n,n,n}
 
  Int_t                      fnBinsPt;                  //number of pt bins
  Double_t                   fPtMin;                    //minimum pt   
  Double_t                   fPtMax;                    //maximum pt    
  
  Int_t                      fnBinsEta;                 //number of eta bins
  Double_t                   fEtaMin;                   //minimum eta   
  Double_t                   fEtaMax;                   //maximum eta
  Int_t                      fEventCounter;             //counting the number of events    
   
  Bool_t                     fUsePhiWeights;            // phi weights
  Bool_t                     fUsePtWeights;             // pt weights
  Bool_t                     fUseEtaWeights;            // eta weights
  Bool_t                     fUseWeights;               // use phi || pt || eta weights
  TBits*                     fUseWeightsBits;           // use phi || pt || eta weights 
  
  // ...................................................................................................................  
  // Q_{n,k} and S^M_{n,k}:        
  TMatrixD *fReQ;  // real part of the Q-vectors stored in matrix fReQ[n][k] = sum_{i=1}^{M} w_{i}^{k} cos(n phi_{i})
  TMatrixD *fImQ;  // imaginary part of the Q-vectors stored in matrix fImQ[n][k] = sum_{i=1}^{M} w_{i}^{k} sin(n phi_{i})
  TMatrixD *fSMpk; // fSM[p][k] = (sum_{i=1}^{M} w_{i}^{k})^{p}
  
  // q_{n} (POIs):
  TH2D *fReqnPtEta; // real part of q_n (q_n is a Q-vector evaluated only for POIs in harmonic n for each (pt,eta) bin)  
  TH2D *fImqnPtEta; // imaginary part of q_n (q_n is a Q-vector evaluated only for POIs in harmonic n for each (pt,eta) bin)
  TH2D *fmPtEta;    // # of POIs (m) for each (pt,eta) bin  
  
  // non-weighted q''_{n} and q''_{2n} (both POIs and RPs)
  TH2D *fReqPrimePrime1nPtEta; // real part of q''_{n} for each (pt,eta) bin  
  TH2D *fImqPrimePrime1nPtEta; // imaginary part of q''_{n} for each (pt,eta) bin 
  TH2D *fReqPrimePrime2nPtEta; // real part of q''_{2n} for each (pt,eta) bin
  TH2D *fImqPrimePrime2nPtEta; // imaginary part of q''_{2n} for each (pt,eta) bin
  
  // weighted q''_{n,2k} and q''_{2n,k} (both POIs and RPs)
  TH2D *fReqPrimePrime1n2kPtEta; // real part of q''_{n,2k} for each (pt,eta) bin  
  TH2D *fImqPrimePrime1n2kPtEta; // imaginary part of q''_{n,2k} for each (pt,eta) bin  
  TH2D *fReqPrimePrime2n1kPtEta; // real part of q''_{2n,k} for each (pt,eta) bin 
  TH2D *fImqPrimePrime2n1kPtEta; // imaginary part of q''_{2n,k} for each (pt,eta) bin 
  
  // m'' (both POIs and RPs) :
  TH2D *fmPrimePrimePtEta; // # of particles which are both POIs and RPs for each (pt,eta) bin
  
  // S^{m''}_{p,k} (both POIs and RPs):
  TH2D *fSmPrimePrime1p1kPtEta; // pow(sum_{i=1}^{m''} w_{i} cos(n phi_{i}), 1)
  TH2D *fSmPrimePrime1p2kPtEta; // pow(sum_{i=1}^{m''} w_{i}^{2} cos(n phi_{i}), 1)
  TH2D *fSmPrimePrime1p3kPtEta; // pow(sum_{i=1}^{m''} w_{i}^{3} cos(n phi_{i}), 1)
  
  // non-weighted q_RP{n} and q_RP{2n} (for each (pt,eta) bin for RPs)
  TH2D *fReqRP1nPtEta; // real part of q_RP{n} (q_RP{n} is a Q-vector evaluated only for RPs in harmonic n for each (pt,eta) bin)
  TH2D *fImqRP1nPtEta; // imaginary part of q_RP{n} (q_RP{n} is a Q-vector evaluated only for RPs in harmonic n for each (pt,eta) bin)
  TH2D *fReqRP2nPtEta; // real part of q_RP{2n} (q_RP{2n} is a Q-vector evaluated only for RPs in harmonic 2n for each (pt,eta) bin)
  TH2D *fImqRP2nPtEta; // imaginary part of q_RP{2n} (q_RP{2n} is a Q-vector evaluated only for RPs in harmonic 2n for each (pt,eta) bin)
  
  // weighted q_RP{n,2k} and q_RP{2n,k} (for each (pt,eta) bin for RPs)
  TH2D *fReqRP1n2kPtEta; // real part of q_RP{n,2k} for each (pt,eta) bin  
  TH2D *fImqRP1n2kPtEta; // imaginary part of q_RP{n,2k} for each (pt,eta) bin  
  TH2D *fReqRP2n1kPtEta; // real part of q_RP{2n,k} for each (pt,eta) bin 
  TH2D *fImqRP2n1kPtEta; // imaginary part of q_RP{2n,k} for each (pt,eta) bin 
  
  // m_RP:
  TH2D *fmRPPtEta; // # of particles which are RPs for each (pt,eta) bin
  
  // S^{m_RP}_{p,k} (for each (pt,eta) bin for RPs):
  TH2D *fSmRP1p1kPtEta; // pow(sum_{i=1}^{m_RP} w_{i} cos(n phi_{i}), 1)
  TH2D *fSmRP1p2kPtEta; // pow(sum_{i=1}^{m_RP} w_{i}^{2} cos(n phi_{i}), 1)
  TH2D *fSmRP1p3kPtEta; // pow(sum_{i=1}^{m_RP} w_{i}^{3} cos(n phi_{i}), 1)
  
  // ----- RESULTS ----
  
  // non-weighted integrated flow:
  TH1D *fIntFlowResultsQC;     // final results for non-weighted no-name integrated flow
  TH1D *fIntFlowResultsPOIQC;  // final results for non-weighted POIs integrated flow
  TH1D *fIntFlowResultsRPQC;   // final results for non-weighted RPs integrated flow
  
  // weighted integrated flow:
  TH1D *fIntFlowResultsQCW;    // final results for weighted no-name integrated flow
  TH1D *fIntFlowResultsPOIQCW; // final results for weighted POIs integrated flow
  TH1D *fIntFlowResultsRPQCW;  // final results for weighted RPs integrated flow
  
  // non-weighted correlations for each (pt,eta) bin for POIs:
  TProfile2D *f2pPtEtaPOI; // <cos n(psi1-phi2)> for POIs
  TProfile2D *f4pPtEtaPOI; // <cos n(psi1+phi2-phi3-phi4)> for POIs 
  TProfile2D *f6pPtEtaPOI; // <cos n(psi1+phi2+phi3-phi4-phi5-phi6)> for POIs 
  TProfile2D *f8pPtEtaPOI; // <cos n(psi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)> for POIs 
  
  // non-weighted final results for differential flow for POIs:
  // 3D (pt,eta):
  TH2D *fvn2ndPtEtaPOI; // v'_{n}{2,QC} (pt,eta) for POIs
  TH2D *fvn4thPtEtaPOI; // v'_{n}{4,QC} (pt,eta) for POIs
  TH2D *fvn6thPtEtaPOI; // v'_{n}{6,QC} (pt,eta) for POIs
  TH2D *fvn8thPtEtaPOI; // v'_{n}{8,QC} (pt,eta) for POIs
  // 2D (pt):
  TH1D *fvn2ndPtPOI; // v'_{n}{2,QC} (pt) for POIs
  TH1D *fvn4thPtPOI; // v'_{n}{4,QC} (pt) for POIs
  TH1D *fvn6thPtPOI; // v'_{n}{6,QC} (pt) for POIs
  TH1D *fvn8thPtPOI; // v'_{n}{8,QC} (pt) for POIs
  // 2D (eta):
  TH1D *fvn2ndEtaPOI; // v'_{n}{2,QC} (eta) for POIs
  TH1D *fvn4thEtaPOI; // v'_{n}{4,QC} (eta) for POIs
  TH1D *fvn6thEtaPOI; // v'_{n}{6,QC} (eta) for POIs
  TH1D *fvn8thEtaPOI; // v'_{n}{8,QC} (eta) for POIs

  // weighted correlations for each (pt,eta) bin for POIs:
  TProfile2D *f2pPtEtaPOIW; // <w2 cos n(psi1-phi2)> for POIs
  TProfile2D *f4pPtEtaPOIW; // <w2 w3 w4 cos n(psi1+phi2-phi3-phi4)> for POIs 
  TProfile2D *f6pPtEtaPOIW; // <w2 w3 w4 w5 w6 cos n(psi1+phi2+phi3-phi4-phi5-phi6)> for POIs 
  TProfile2D *f8pPtEtaPOIW; // <w2 w3 w4 w5 w6 w7 w8 cos n(psi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)> for POIs 
  
  // weighted final results for differential flow  for POIs:
  // 3D (pt,eta):
  TH2D *fvn2ndPtEtaPOIW; // v'_{n}{2,QC} (pt,eta) for POIs
  TH2D *fvn4thPtEtaPOIW; // v'_{n}{4,QC} (pt,eta) for POIs
  TH2D *fvn6thPtEtaPOIW; // v'_{n}{6,QC} (pt,eta) for POIs
  TH2D *fvn8thPtEtaPOIW; // v'_{n}{8,QC} (pt,eta) for POIs
  // 2D (pt):
  TH1D *fvn2ndPtPOIW; // v'_{n}{2,QC} (pt) for POIs
  TH1D *fvn4thPtPOIW; // v'_{n}{4,QC} (pt) for POIs
  TH1D *fvn6thPtPOIW; // v'_{n}{6,QC} (pt) for POIs
  TH1D *fvn8thPtPOIW; // v'_{n}{8,QC} (pt) for POIs
  // 2D (eta):
  TH1D *fvn2ndEtaPOIW; // v'_{n}{2,QC} (eta) for POIs
  TH1D *fvn4thEtaPOIW; // v'_{n}{4,QC} (eta) for POIs
  TH1D *fvn6thEtaPOIW; // v'_{n}{6,QC} (eta) for POIs
  TH1D *fvn8thEtaPOIW; // v'_{n}{8,QC} (eta) for POIs
  
  // non-weighted correlations for each (pt,eta) bin for RPs:
  TProfile2D *f2pPtEtaRP; // <cos n(psi1-phi2)> for RPs
  TProfile2D *f4pPtEtaRP; // <cos n(psi1+phi2-phi3-phi4)> for RPs 
  TProfile2D *f6pPtEtaRP; // <cos n(psi1+phi2+phi3-phi4-phi5-phi6)> for RPs 
  TProfile2D *f8pPtEtaRP; // <cos n(psi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)> for RPs
  
  // non-weighted final results for differential flow for RPs:
  // 3D (pt,eta):
  TH2D *fvn2ndPtEtaRP; // v'_{n}{2,QC} (pt,eta) for RPs
  TH2D *fvn4thPtEtaRP; // v'_{n}{4,QC} (pt,eta) for RPs
  TH2D *fvn6thPtEtaRP; // v'_{n}{6,QC} (pt,eta) for RPs
  TH2D *fvn8thPtEtaRP; // v'_{n}{8,QC} (pt,eta) for RPs
  // 2D (pt):
  TH1D *fvn2ndPtRP; // v'_{n}{2,QC} (pt) for RPs
  TH1D *fvn4thPtRP; // v'_{n}{4,QC} (pt) for RPs
  TH1D *fvn6thPtRP; // v'_{n}{6,QC} (pt) for RPs
  TH1D *fvn8thPtRP; // v'_{n}{8,QC} (pt) for RPs
  // 2D (eta):
  TH1D *fvn2ndEtaRP; // v'_{n}{2,QC} (eta) for RPs
  TH1D *fvn4thEtaRP; // v'_{n}{4,QC} (eta) for RPs
  TH1D *fvn6thEtaRP; // v'_{n}{6,QC} (eta) for RPs
  TH1D *fvn8thEtaRP; // v'_{n}{8,QC} (eta) for RPs
 
  // weighted correlations for each (pt,eta) bin for RPs:
  TProfile2D *f2pPtEtaRPW; // <w2 cos n(psi1-phi2)> for RPs
  TProfile2D *f4pPtEtaRPW; // <w2 w3 w4 cos n(psi1+phi2-phi3-phi4)> for RPs 
  TProfile2D *f6pPtEtaRPW; // <w2 w3 w4 w5 w6 cos n(psi1+phi2+phi3-phi4-phi5-phi6)> for RPs 
  TProfile2D *f8pPtEtaRPW; // <w2 w3 w4 w5 w6 w7 w8 cos n(psi1+phi2+phi3+phi4-phi5-phi6-phi7-phi8)> for RPs 
  
  // weighted final results for differential flow for RPs:
  // 3D (pt,eta):
  TH2D *fvn2ndPtEtaRPW; // v'_{n}{2,QC} (pt,eta) for RPs
  TH2D *fvn4thPtEtaRPW; // v'_{n}{4,QC} (pt,eta) for RPs
  TH2D *fvn6thPtEtaRPW; // v'_{n}{6,QC} (pt,eta) for RPs
  TH2D *fvn8thPtEtaRPW; // v'_{n}{8,QC} (pt,eta) for RPs
  // 2D (pt):
  TH1D *fvn2ndPtRPW; // v'_{n}{2,QC} (pt) for RPs
  TH1D *fvn4thPtRPW; // v'_{n}{4,QC} (pt) for RPs
  TH1D *fvn6thPtRPW; // v'_{n}{6,QC} (pt) for RPs
  TH1D *fvn8thPtRPW; // v'_{n}{8,QC} (pt) for RPs
  // 2D (eta):
  TH1D *fvn2ndEtaRPW; // v'_{n}{2,QC} (eta) for RPs
  TH1D *fvn4thEtaRPW; // v'_{n}{4,QC} (eta) for RPs
  TH1D *fvn6thEtaRPW; // v'_{n}{6,QC} (eta) for RPs
  TH1D *fvn8thEtaRPW; // v'_{n}{8,QC} (eta) for RPs
  // ...................................................................................................................
  
  
  
                       
  ClassDef(AliFlowAnalysisWithQCumulants, 0);
};

//================================================================================================================

#endif





