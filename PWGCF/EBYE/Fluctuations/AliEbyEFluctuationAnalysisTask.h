#ifndef AliEbyEFluctuationAnalysisTask_cxx
#define AliEbyEFluctuationAnalysisTask_cxx

// Event by event charge fluctuation analysis
// Authors: Satyajit Jena and Panos Cristakoglou

class TH1F;
class TH2F;
class TString;
class AliESDEvent;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

const Int_t nCentralityBins = 20;

class AliEbyEFluctuationAnalysisTask : public AliAnalysisTaskSE {
 public:
  AliEbyEFluctuationAnalysisTask() : AliAnalysisTaskSE(), fESD(0), fOutputList(0), fHistEventStats(0), fHistCentrality(0), fHistNMultMC(0), fHistNPlusNMinusMC(0), fESDtrackCuts(0), fAnalysisType(0), fAnalysisMode(0), fCentralityEstimator("V0M"), fCentralityBins20(kFALSE), fVxMax(3.0),fVyMax(3.0), fVzMax(10.) {
    for(Int_t iBin = 0; iBin < nCentralityBins; iBin++) {
      fHistNMult[iBin] = NULL;
      fHistNPlusNMinus[iBin] = NULL;
    }
  }
  AliEbyEFluctuationAnalysisTask(const char *name);
  virtual ~AliEbyEFluctuationAnalysisTask() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
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
  void SetCentralityBins20() {fCentralityBins20 = kTRUE;}

  void SetAnalysisType(const char* analysisType) {
    fAnalysisType = analysisType;}
  void SetAnalysisMode(const char* analysisMode) {
    fAnalysisMode = analysisMode;}

 private:
  AliESDEvent *fESD;    //! ESD object
  TList       *fOutputList; //! Output list
  TH1F        *fHistEventStats; //!event stats
  TH1F        *fHistCentrality; //!centrality
  TH1F        *fHistNMult[nCentralityBins]; //! nmult
  TH2F        *fHistNPlusNMinus[nCentralityBins];//!nplus vs nminus correlation
  TH1F        *fHistNMultMC; //!nmult MC
  TH2F        *fHistNPlusNMinusMC;//!nplus vs nminus correlation

  AliESDtrackCuts *fESDtrackCuts; //ESD track cuts

  TString fAnalysisType;//"MC", "ESD", "AOD"
  TString fAnalysisMode;//"TPC", "Global"

  TString fCentralityEstimator;//"V0M","TRK","TKL","ZDC","FMD"
  Bool_t fCentralityBins20;//centrality bins of 5% width

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  AliEbyEFluctuationAnalysisTask(const AliEbyEFluctuationAnalysisTask&); // not implemented
  AliEbyEFluctuationAnalysisTask& operator=(const AliEbyEFluctuationAnalysisTask&); // not implemented
  
  ClassDef(AliEbyEFluctuationAnalysisTask, 1); // example of analysis
};

#endif
