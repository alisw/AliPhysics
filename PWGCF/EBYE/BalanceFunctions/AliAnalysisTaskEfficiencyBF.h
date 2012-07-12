#ifndef ALIANALYSISTASKEFFICIENCYBF_cxx
#define ALIANALYSISTASKEFFICIENCYBF_cxx

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency of the Balance Function 
// for single particles and pairs
//
// Authors: Panos Christakoglou, Michael Weber
// 
// ---------------------------------------------------------------------

class TH1F;
class TH2F;
class TString;
class AliESDEvent;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEfficiencyBF : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEfficiencyBF() : AliAnalysisTaskSE(), 
    fESD(0), fQAList(0), fOutputList(0), 
    fHistEventStats(0), fHistCentrality(0), fHistNMult(0), 
    fHistGeneratedEtaPtPlus(0), fHistFindableEtaPtPlus(0), 
    fHistReconstructedEtaPtPlus(0), fHistSurvivedEtaPtPlus(0),
    fHistGeneratedEtaPtMinus(0), fHistFindableEtaPtMinus(0), 
    fHistReconstructedEtaPtMinus(0), fHistSurvivedEtaPtMinus(0),
    fHistGeneratedEtaPtPlusControl(0), fHistFindableEtaPtPlusControl(0), 
    fHistReconstructedEtaPtPlusControl(0), fHistSurvivedEtaPtPlusControl(0),
    fHistGeneratedEtaPtMinusControl(0), fHistFindableEtaPtMinusControl(0), 
    fHistReconstructedEtaPtMinusControl(0), fHistSurvivedEtaPtMinusControl(0),
    fHistGeneratedEtaPtPlusPlus(0), fHistFindableEtaPtPlusPlus(0), 
    fHistReconstructedEtaPtPlusPlus(0), fHistSurvivedEtaPtPlusPlus(0),
    fHistGeneratedEtaPtMinusMinus(0), fHistFindableEtaPtMinusMinus(0), 
    fHistReconstructedEtaPtMinusMinus(0), fHistSurvivedEtaPtMinusMinus(0),
    fHistGeneratedEtaPtPlusMinus(0), fHistFindableEtaPtPlusMinus(0), 
    fHistReconstructedEtaPtPlusMinus(0), fHistSurvivedEtaPtPlusMinus(0),
    fESDtrackCuts(0), fAnalysisMode(0), 
    fCentralityEstimator("V0M"), fCentralityPercentileMin(0.0), fCentralityPercentileMax(5.0), 
    fVxMax(3.0), fVyMax(3.0), fVzMax(10.), 
    fMinNumberOfTPCClusters(80), fMaxChi2PerTPCCluster(4.0), fMaxDCAxy(3.0), fMaxDCAz(3.0),
    fMinPt(0.3), fMaxPt(1.5), fMaxEta(0.8) {}
  AliAnalysisTaskEfficiencyBF(const char *name);
  virtual ~AliAnalysisTaskEfficiencyBF() {}
  
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

  
 private:
  AliESDEvent *fESD;    //! ESD object
  TList       *fQAList; //! QA list
  TList       *fOutputList; //! Output list
  
  // QA histograms
  TH1F        *fHistEventStats; //!event stats
  TH1F        *fHistCentrality; //!centrality
  TH1F        *fHistNMult; //! nmult   

  // output histograms (single particles)
  TH2F        *fHistGeneratedEtaPtPlus;//!correction map for positives (generated)
  TH2F        *fHistFindableEtaPtPlus;//!correction map for positives (findable)
  TH2F        *fHistReconstructedEtaPtPlus;//!correction map for positives (reconstructed)
  TH2F        *fHistSurvivedEtaPtPlus;//!correction map positives (survived)

  TH2F        *fHistGeneratedEtaPtMinus;//!correction map for negatives (generated)
  TH2F        *fHistFindableEtaPtMinus;//!correction map for negatives (findable)
  TH2F        *fHistReconstructedEtaPtMinus;//!correction map for negatives (reconstructed)
  TH2F        *fHistSurvivedEtaPtMinus;//!correction map negatives (survived)

  TH2F        *fHistGeneratedEtaPtPlusControl;//!correction map for positives (generated)
  TH2F        *fHistFindableEtaPtPlusControl;//!correction map for positives (findable)
  TH2F        *fHistReconstructedEtaPtPlusControl;//!correction map for positives (reconstructed)
  TH2F        *fHistSurvivedEtaPtPlusControl;//!correction map positives (survived)

  TH2F        *fHistGeneratedEtaPtMinusControl;//!correction map for negatives (generated)
  TH2F        *fHistFindableEtaPtMinusControl;//!correction map for negatives (findable)
  TH2F        *fHistReconstructedEtaPtMinusControl;//!correction map for negatives (reconstructed)
  TH2F        *fHistSurvivedEtaPtMinusControl;//!correction map negatives (survived)

  // output histograms (pairs)
  TH2F        *fHistGeneratedEtaPtPlusPlus;//!correction map for ++ (generated)
  TH2F        *fHistFindableEtaPtPlusPlus;//!correction map for ++ (findable)
  TH2F        *fHistReconstructedEtaPtPlusPlus;//!correction map for ++ (reconstructed)
  TH2F        *fHistSurvivedEtaPtPlusPlus;//!correction map ++ (survived)

  TH2F        *fHistGeneratedEtaPtMinusMinus;//!correction map for -- (generated)
  TH2F        *fHistFindableEtaPtMinusMinus;//!correction map for -- (findable)
  TH2F        *fHistReconstructedEtaPtMinusMinus;//!correction map for -- (reconstructed)
  TH2F        *fHistSurvivedEtaPtMinusMinus;//!correction map -- (survived)

  TH2F        *fHistGeneratedEtaPtPlusMinus;//!correction map for +- (generated)
  TH2F        *fHistFindableEtaPtPlusMinus;//!correction map for +- (findable)
  TH2F        *fHistReconstructedEtaPtPlusMinus;//!correction map for +- (reconstructed)
  TH2F        *fHistSurvivedEtaPtPlusMinus;//!correction map +- (survived)



  AliESDtrackCuts *fESDtrackCuts; //ESD track cuts

  TString fAnalysisMode;//"TPC", "Global"

  TString fCentralityEstimator;//"V0M","TRK","TKL","ZDC","FMD"
  Float_t fCentralityPercentileMin, fCentralityPercentileMax; //min-max centrality percentile

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  Double_t fMinNumberOfTPCClusters; //
  Double_t fMaxChi2PerTPCCluster; //
  Double_t fMaxDCAxy, fMaxDCAz;//
  Double_t fMinPt, fMaxPt;
  Double_t fMaxEta;

  AliAnalysisTaskEfficiencyBF(const AliAnalysisTaskEfficiencyBF&); // not implemented
  AliAnalysisTaskEfficiencyBF& operator=(const AliAnalysisTaskEfficiencyBF&); // not implemented
  
  ClassDef(AliAnalysisTaskEfficiencyBF, 1); // example of analysis
};

#endif
