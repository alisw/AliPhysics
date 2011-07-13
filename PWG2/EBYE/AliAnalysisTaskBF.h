#ifndef ALIANALYSISTASKBF_CXX
#define ALIANALYSISTASKBF_CXX

// Analysis task for the BF code
// Authors: Panos Cristakoglou@cern.ch

class TList;
class TH1F;
class TH2F;

class AliBalance;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"
#include "AliBalance.h"


class AliAnalysisTaskBF : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskBF(const char *name = "AliAnalysisTaskBF");
  virtual ~AliAnalysisTaskBF() {}
  
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);

  void SetAnalysisObject(AliBalance *const analysis) {
    fBalance = analysis;}
  void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
    fESDtrackCuts = trackCuts;}
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;
    fVyMax = vy;
    fVzMax = vz;
  }
  void SetKinematicsCutsAOD(Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax){
    fPtMin  = ptmin;
    fPtMax  = ptmax;
    fEtaMin = etamin;
    fEtaMax = etamax;

  }
 
  //Centrality
  void SetCentralityEstimator(const char* centralityEstimator) {fCentralityEstimator = centralityEstimator;}
  const char* GetCentralityEstimator(void)                     {return fCentralityEstimator;}
  void SetCentralityPercentileRange(Double_t min, Double_t max) { 
    fCentralityPercentileMin=min;
    fCentralityPercentileMax=max;
  }

  void UseOfflineTrigger() {fUseOfflineTrigger = kTRUE;}

 private:
  AliBalance *fBalance; //BF object
  TList *fList; //fList object
  TList *fListBF; //fList object

  TH1F *fHistEventStats; //event stats
  TH1F *fHistTrackStats; //Track filter bit stats
  TH1F *fHistVx; //x coordinate of the primary vertex
  TH1F *fHistVy; //y coordinate of the primary vertex
  TH1F *fHistVz; //z coordinate of the primary vertex

  TH2F *fHistClus;
  TH2F *fHistDCA;
  TH1F *fHistChi2;
  TH1F *fHistPt;
  TH1F *fHistEta;
  TH1F *fHistPhi;
  TH2F *fHistV0M;

  TH1F *fHistBF[NUMBER_OF_ANALYSES][3];
  TH1F *fHistN;

  AliESDtrackCuts *fESDtrackCuts; //ESD track cuts

  TString fCentralityEstimator;      //"V0M","TRK","TKL","ZDC","FMD"
  Double_t fCentralityPercentileMin;
  Double_t fCentralityPercentileMax;

  Bool_t fUseOfflineTrigger;//Usage of the offline trigger selection

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  Double_t fPtMin;//only used for AODs
  Double_t fPtMax;//only used for AODs
  Double_t fEtaMin;//only used for AODs
  Double_t fEtaMax;//only used for AODs

  AliAnalysisTaskBF(const AliAnalysisTaskBF&); // not implemented
  AliAnalysisTaskBF& operator=(const AliAnalysisTaskBF&); // not implemented
  
  ClassDef(AliAnalysisTaskBF, 5); // example of analysis
};

#endif
