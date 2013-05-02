#ifndef ALIANALYSISTASKEFFCONTBF_cxx
#define ALIANALYSISTASKEFFCONTBF_cxx

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency and contamination of the Balance 
// Function for single particles and pairs
// 
// ---------------------------------------------------------------------

class TList;
class TH1F;
class TH3D;
class TH2F;
class TString;
class AliAODEvent;
class AliAODInputHandler;
class TH2D;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEffContBF : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEffContBF() : AliAnalysisTaskSE(), 
    fAOD(0), 
    fAODtrackCutBit(128), 
    fArrayMC(0),  
    fQAList(0), 
    fOutputList(0),
    fHistEventStats(0),
    fHistCentrality(0), 
    fHistNMult(0), 
    fHistVz(0),
    fHistContaminationSecondaries(0),
    fHistContaminationPrimaries(0),
    fHistGeneratedEtaPtPhiPlus(0),
    fHistSurvivedEtaPtPhiPlus(0),
    fHistGeneratedEtaPtPhiMinus(0),
    fHistSurvivedEtaPtPhiMinus(0),
    fHistGeneratedEtaPtPlusControl(0), 
    fHistSurvivedEtaPtPlusControl(0),
    fHistGeneratedEtaPtMinusControl(0), 
    fHistSurvivedEtaPtMinusControl(0), 
    fHistGeneratedEtaPtPlusPlus(0), 
    fHistSurvivedEtaPtPlusPlus(0),
    fHistGeneratedEtaPtMinusMinus(0), 
    fHistSurvivedEtaPtMinusMinus(0),
    fHistGeneratedEtaPtPlusMinus(0), 
    fHistSurvivedEtaPtPlusMinus(0),
    fHistGeneratedPhiEtaPlusPlus(0),
    fHistSurvivedPhiEtaPlusPlus(0),
    fHistGeneratedPhiEtaMinusMinus(0), 
    fHistSurvivedPhiEtaMinusMinus(0),
    fHistGeneratedPhiEtaPlusMinus(0),
    fHistSurvivedPhiEtaPlusMinus(0), 
    fAnalysisMode(0), 
    fCentralityEstimator("V0M"), 
    fCentralityPercentileMin(0.0), 
    fCentralityPercentileMax(5.0), 
    fVxMax(3.0), 
    fVyMax(3.0), 
    fVzMax(10.), 
    fMinNumberOfTPCClusters(80), 
    fMaxChi2PerTPCCluster(4.0),
    fMaxDCAxy(3.0), 
    fMaxDCAz(3.0),
    fMinPt(0.0), 
    fMaxPt(20.0), 
    fMinEta(-0.8),
    fMaxEta(0.8), 
    fEtaRangeMin(0.0), 
    fEtaRangeMax(1.6), 
    fPtRangeMin(0.0), 
    fPtRangeMax(20.0), 
    fPhiRangeMin(0.0),
    fPhiRangeMax(360.), 
    fdPhiRangeMax(180.), 
    fEtaBin(100),
    fdEtaBin(64),
    fPtBin(100),
    fPhiBin(100),
    fdPhiBin(90), 
    fHistSurvived4EtaPtPhiPlus(0), 
    fHistSurvived8EtaPtPhiPlus(0){}
    AliAnalysisTaskEffContBF(const char *name);
    virtual ~AliAnalysisTaskEffContBF() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  Bool_t   IsLabelUsed(TArrayI array, Int_t label);

  void SetAODtrackCutBit(Int_t bit){
    fAODTrackCutBit = bit;
  }
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
 
  void SetEtaRange(Double_t minEta, Double_t maxEta, Int_t binEta, Double_t minRangeEta, Double_t maxRangeEta, Int_t bindEta){
    fMinEta = minEta;
    fMaxEta = maxEta;
    fEtaBin = binEta;
    fEtaRangeMax = maxRangeEta;
    fEtaRangeMin = minRangeEta;
    fdEtaBin = bindEta;
  }
  void SetPtRange(Double_t minRangePt, Double_t maxRangePt,Int_t binPt){
    fPtRangeMin = minRangePt;
    fPtRangeMax = maxRangePt;
    fPtBin = binPt;}  
  void SetPhiRange(Double_t minRangePhi, Double_t maxRangePhi,Int_t binPhi,Double_t maxRangedPhi,Int_t bindPhi ){
    fPhiRangeMin = minRangePhi;
    fPhiRangeMax = maxRangePhi;
    fPhiBin = binPhi;
    fdPhiRangeMax = maxRangedPhi;
    fdPhiBin = bindPhi;
  } 

 private:
  AliAODEvent* fAOD; //! AOD object  
  Int_t fAODtrackCutBit;  
  TClonesArray *fArrayMC; //! array of MC particles  
  TList       *fQAList; //! QA list
  TList       *fOutputList; //! Output list
  
  // QA histograms
  TH1F        *fHistEventStats; //!event stats
  TH1F        *fHistCentrality; //!centrality
  TH1F        *fHistNMult; //! nmult  
  TH1F        *fHistVz;//!

  // output histograms
  TH3D        *fHistContaminationSecondaries;//!
  TH3D        *fHistContaminationPrimaries;//!
  // output histograms (single particles)
  TH3D        *fHistGeneratedEtaPtPhiPlus;//!correction map for positives (generated)
  TH3D        *fHistSurvivedEtaPtPhiPlus;//!correction map positives (survived)
 
  TH3D        *fHistGeneratedEtaPtPhiMinus;//!correction map for negatives (generated)
  TH3D        *fHistSurvivedEtaPtPhiMinus;//!correction map negatives (survived)
 
  TH2F        *fHistGeneratedEtaPtPlusControl;//!correction map for positives (generated)
  TH2F        *fHistSurvivedEtaPtPlusControl;//!correction map positives (survived)
 
  TH2F        *fHistGeneratedEtaPtMinusControl;//!correction map for negatives (generated)
  TH2F        *fHistSurvivedEtaPtMinusControl;//!correction map negatives (survived)
 
  // output histograms (pairs)
  TH2F        *fHistGeneratedEtaPtPlusPlus;//!correction map for ++ (generated)
  TH2F        *fHistSurvivedEtaPtPlusPlus;//!correction map ++ (survived)
 
  TH2F        *fHistGeneratedEtaPtMinusMinus;//!correction map for -- (generated)
  TH2F        *fHistSurvivedEtaPtMinusMinus;//!correction map -- (survived)
 
  TH2F        *fHistGeneratedEtaPtPlusMinus;//!correction map for +- (generated)
  TH2F        *fHistSurvivedEtaPtPlusMinus;//!correction map +- (survived)

  TH2F        *fHistGeneratedPhiEtaPlusPlus;//!correction map for ++ (generated)
  TH2F        *fHistSurvivedPhiEtaPlusPlus;//!correction map ++ (survived)
 
  TH2F        *fHistGeneratedPhiEtaMinusMinus;//!correction map for -- (generated)
  TH2F        *fHistSurvivedPhiEtaMinusMinus;//!correction map -- (survived)
 
  TH2F        *fHistGeneratedPhiEtaPlusMinus;//!correction map for +- (generated)
  TH2F        *fHistSurvivedPhiEtaPlusMinus;//!correction map +- (survived)

  TString fAnalysisMode;//"TPC", "Global"

  TString fCentralityEstimator;//"V0M","TRK","TKL","ZDC","FMD"
  Float_t fCentralityPercentileMin, fCentralityPercentileMax; //min-max centrality percentile

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax
  
  Int_t fAODTrackCutBit;//track cut bit from track selection (only used for AODs)

  Double_t fMinNumberOfTPCClusters;//!
  Double_t fMaxChi2PerTPCCluster;//!
  Double_t fMaxDCAxy, fMaxDCAz;//!
  Double_t fMinPt, fMaxPt;//!
  Double_t fMinEta, fMaxEta;//!
  Double_t fEtaRangeMin;// acceptance cuts 
  Double_t fEtaRangeMax; // acceptance cuts
  Double_t fPtRangeMin;  //acceptance cuts
  Double_t fPtRangeMax;  //acceptance cuts
  Double_t fPhiRangeMin; //acceptance cuts
  Double_t fPhiRangeMax; // acceptance cuts
  Double_t fdPhiRangeMax; // acceptance cuts
  
  Int_t fEtaBin;  //acceptance cuts
  Int_t fdEtaBin;  //acceptance cuts
  Int_t fPtBin; //acceptance cuts
  Int_t fPhiBin; // acceptance cuts
  Int_t fdPhiBin; // acceptance cuts

  TH3F        *fHistSurvived4EtaPtPhiPlus;//!
  TH3F        *fHistSurvived8EtaPtPhiPlus;//!

  AliAnalysisTaskEffContBF(const AliAnalysisTaskEffContBF&); // not implemented
  AliAnalysisTaskEffContBF& operator=(const AliAnalysisTaskEffContBF&); // not implemented
  
  ClassDef(AliAnalysisTaskEffContBF, 1); // example of analysis
};

#endif
