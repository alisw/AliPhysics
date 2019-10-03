#ifndef ALIANALYSISTASKEFFICIENCYBFPSI_cxx
#define ALIANALYSISTASKEFFICIENCYBFPSI_cxx

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency of the Balance Function 
// for single particles and pairs: multi-D analysis
//
// Authors: Panos Christakoglou, Michael Weber
// 
// ---------------------------------------------------------------------

class TH1F;
class TH3D;
class TH2F;
class TString;
class AliESDEvent;
class AliESDtrackCuts;
#include "AliTHn.h"

#include "AliAnalysisTaskSE.h"

const Int_t kVariablesSingle = 2;       // track variables in histogram (phi-Psi2, pTtrig)
const Int_t kVariablesPair   = 5;       // track variables in histogram (phi-Psi2, dEta, dPhi, pTtrig, ptAssociated)

class AliAnalysisTaskEfficiencyBFPsi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEfficiencyBFPsi() : AliAnalysisTaskSE(), 
    fESD(0), fQAList(0), fOutputList(0), 
    fHistEventStats(0), fHistCentrality(0), fHistNMult(0), 
    fHistGeneratedPlus(0), fHistSurvivedPlus(0),
    fHistGeneratedMinus(0),fHistSurvivedMinus(0),
    fHistGeneratedPlusPlus(0), fHistSurvivedPlusPlus(0),
    fHistGeneratedMinusMinus(0), fHistSurvivedMinusMinus(0),
    fHistGeneratedPlusMinus(0), fHistSurvivedPlusMinus(0),
    fHistGeneratedMinusPlus(0), fHistSurvivedMinusPlus(0),
    fESDtrackCuts(0), fAnalysisMode(0), 
    fCentralityEstimator("V0M"), 
    fCentralityPercentileMin(0.0), fCentralityPercentileMax(10.0), 
    fVxMax(3.0), fVyMax(3.0), fVzMax(10.), 
    fMinNumberOfTPCClusters(80), fMaxChi2PerTPCCluster(4.0), 
    fMaxDCAxy(3.0), fMaxDCAz(3.0),
    fMinPt(0.3), fMaxPt(1.5), fMaxEta(0.8), fEtaRangeMax(0.8), 
    fPtRangeMin(0.1), fPtRangeMax(5.0), fPhiRangeMin(0.0),fPhiRangeMax(6.28){}
  AliAnalysisTaskEfficiencyBFPsi(const char *name);
  virtual ~AliAnalysisTaskEfficiencyBFPsi() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  Bool_t   IsLabelUsed(TArrayI array, Int_t label);

  void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
    fESDtrackCuts = trackCuts;}
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;
    fVyMax = vy;
    fVzMax = vz;
  }
 
  //Centrality
  void SetCentralityEstimator(const char* centralityEstimator) {
    fCentralityEstimator = centralityEstimator;}
  void SetCentralityPercentileRange(Float_t min, Float_t max) { 
    fCentralityPercentileMin=min;
    fCentralityPercentileMax=max;
  }

  void SetAnalysisMode(const char* analysisMode) {
    fAnalysisMode = analysisMode;}

  //Track cuts
  void SetMinNumberOfTPCClusters(Double_t min) {
    fMinNumberOfTPCClusters = min;}
  void SetMaxChi2PerTPCCluster(Double_t max) {
    fMaxChi2PerTPCCluster = max;}
  void SetMaxDCAxy(Double_t max) {
    fMaxDCAxy = max;}
  void SetMaxDCAz(Double_t max) {
    fMaxDCAz = max;}
  void SetMinPt(Double_t minPt) {
    fMinPt = minPt;}
  void SetMaxPt(Double_t maxPt) {
    fMaxPt = maxPt;}
  void SetMaxEta(Double_t maxEta) {
    fMaxEta = maxEta;}

  void SetEtaRangeMax(Double_t maxRangeEta){
    fEtaRangeMax = maxRangeEta;}//
  void SetPtRangeMin(Double_t minRangePt){
    fPtRangeMin = minRangePt;} // 
  void SetPtRangeMax(Double_t maxRangePt){
    fPtRangeMax = maxRangePt;} //
  void SetPhiRangeMin(Double_t minRangePhi){
    fPhiRangeMin = minRangePhi;} //
  void SetPhiRangeMax(Double_t maxRangePhi){
    fPhiRangeMax = maxRangePhi;} //
  
   
 private:
  AliESDEvent *fESD;    //! ESD object
  TList       *fQAList; //! QA list
  TList       *fOutputList; //! Output list
  
  // QA histograms
  TH1F *fHistEventStats; //!event stats
  TH1F *fHistCentrality; //!centrality
  TH1F *fHistNMult; //! nmult   

  // output histograms (single particles)
  AliTHn *fHistGeneratedPlus;//!correction map for positives (generated)
  AliTHn *fHistSurvivedPlus;//!correction map positives (survived)
  AliTHn *fHistGeneratedMinus;//!correction map for negatives (generated)
  AliTHn *fHistSurvivedMinus;//!correction map negatives (survived)

  // output histograms (pairs)
  AliTHn *fHistGeneratedPlusPlus;//!correction map for ++ (generated)
  AliTHn *fHistSurvivedPlusPlus;//!correction map ++ (survived)
  AliTHn *fHistGeneratedMinusMinus;//!correction map for -- (generated)
  AliTHn *fHistSurvivedMinusMinus;//!correction map -- (survived)
  AliTHn *fHistGeneratedPlusMinus;//!correction map for +- (generated)
  AliTHn *fHistSurvivedPlusMinus;//!correction map +- (survived)
  AliTHn *fHistGeneratedMinusPlus;//!correction map for -+ (generated)
  AliTHn *fHistSurvivedMinusPlus;//!correction map -+ (survived)
  //============//

  AliESDtrackCuts *fESDtrackCuts; //ESD track cuts

  TString fAnalysisMode;//"TPC", "Global"

  TString fCentralityEstimator;//"V0M","TRK","TKL","ZDC","FMD"
  Float_t fCentralityPercentileMin, fCentralityPercentileMax; //min-max centrality percentile

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  Double_t fMinNumberOfTPCClusters;
  Double_t fMaxChi2PerTPCCluster;
  Double_t fMaxDCAxy, fMaxDCAz;
  Double_t fMinPt, fMaxPt;
  Double_t fMaxEta;
  Double_t fEtaRangeMax; // acceptance cuts
  Double_t fPtRangeMin;  //acceptance cuts
  Double_t fPtRangeMax;  //acceptance cuts
  Double_t fPhiRangeMin; //acceptance cuts
  Double_t fPhiRangeMax; // acceptance cuts

  AliAnalysisTaskEfficiencyBFPsi(const AliAnalysisTaskEfficiencyBFPsi&); // not implemented
  AliAnalysisTaskEfficiencyBFPsi& operator=(const AliAnalysisTaskEfficiencyBFPsi&); // not implemented
  
  ClassDef(AliAnalysisTaskEfficiencyBFPsi, 1); // example of analysis
};

#endif
