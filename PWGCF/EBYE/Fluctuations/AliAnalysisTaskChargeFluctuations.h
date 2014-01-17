#ifndef ALIANALYSISTASKCHARGEFLUCTUATIONS_CXX
#define ALIANALYSISTASKCHARGEFLUCTUATIONS_CXX

// Analysis task for the charge fluctuations studies
// Authors: Panos Cristakoglou@cern.ch

class TList;
class TH1F;

class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskChargeFluctuations : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskChargeFluctuations(const char *name = "AliAnalysisTaskChargeFluctuations");
  virtual ~AliAnalysisTaskChargeFluctuations() {}
  
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

  void UseOfflineTrigger() {fUseOfflineTrigger = kTRUE;}

 private:
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

  AliAnalysisTaskChargeFluctuations(const AliAnalysisTaskChargeFluctuations&); // not implemented
  AliAnalysisTaskChargeFluctuations& operator=(const AliAnalysisTaskChargeFluctuations&); // not implemented
  
  ClassDef(AliAnalysisTaskChargeFluctuations, 1); // example of analysis
};

#endif
