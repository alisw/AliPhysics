#ifndef ALIANALYSISTASKBF_CXX
#define ALIANALYSISTASKBF_CXX

// Analysis task for the BF code
// Authors: Panos Cristakoglou@cern.ch

class TList;
class TH1F;

class AliBalance;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskBF : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskBF(const char *name = "AliAnalysisTaskBF");
  virtual ~AliAnalysisTaskBF() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
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

  void UseOfflineTrigger() {fUseOfflineTrigger = kTRUE;}

 private:
  AliBalance *fBalance; //BF object
  TList *fList; //fList object
  TH1F *fHistEventStats; //event stats
  TH1F *fHistVx; //x coordinate of the primary vertex
  TH1F *fHistVy; //y coordinate of the primary vertex
  TH1F *fHistVz; //z coordinate of the primary vertex

  AliESDtrackCuts *fESDtrackCuts; //ESD track cuts

  Bool_t fUseOfflineTrigger;//Usage of the offline trigger selection

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  AliAnalysisTaskBF(const AliAnalysisTaskBF&); // not implemented
  AliAnalysisTaskBF& operator=(const AliAnalysisTaskBF&); // not implemented
  
  ClassDef(AliAnalysisTaskBF, 3); // example of analysis
};

#endif
