#ifndef ALIANALYSISTASKAODFILTERBITQA_H
#define ALIANALYSISTASKAODFILTERBITQA_H

// Analysis task for the QA of AOD track filter bits
// Authors: m.weber@cern.ch

#define gBitMax 16// number of maximum filter bits
#define gNCharge 2// number of charges

class TList;
class TH2D;
class TH3D;


class AliAnalysisTaskAODFilterBitQA : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskAODFilterBitQA(const char *name = "AliAnalysisTaskAODFilterBitQA");
  virtual ~AliAnalysisTaskAODFilterBitQA(); 
   
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);

  void SetFillOnlySecondaries(){
    fillOnlySecondaries = kTRUE;
  }

  void SetFillHFVertexingTracks(TString HFBranchName="D0toKpi", Int_t bitIgnore1 = -1, Int_t bitIgnore2 = -1){
    fillHFVertexingTracks = kTRUE;
    fHFBranchName         = HFBranchName;
    fBitIgnore1           = bitIgnore1;
    fBitIgnore2           = bitIgnore2;
  }
  
  void SetCentralityPercentileRange(Double_t min, Double_t max){
    useCentrality            = kTRUE;
    fCentralityPercentileMin = min;
    fCentralityPercentileMax = max;
  }
  void SetUseMultSelectionFramework(Bool_t multFramework) { fUseMultSelectionFramework = multFramework;} 
  void SetUseUncheckedCentrality(Bool_t unchecked) { fUseUncheckedCentrality = unchecked; }

  void SetPtRange(Double_t min, Double_t max){
    fPtMin = min;
    fPtMax = max;
  }

  void SetEtaRange(Double_t min, Double_t max){
    fEtaMin = min;
    fEtaMax = max;
  }

 private:
  Double_t IsEventAccepted(AliVEvent *event);
  void GetAcceptedTracks(AliVEvent *event, Double_t gCentrality);
  void GetAcceptedHFVertexingTracks(AliVEvent *event, Double_t gCentrality);
  
  AliAnalysisTaskAODFilterBitQA(const AliAnalysisTaskAODFilterBitQA&); // not implemented
  AliAnalysisTaskAODFilterBitQA& operator=(const AliAnalysisTaskAODFilterBitQA&); // not implemented

  TClonesArray *fArrayMC;//MC track array for AODs

  TList *fListQA;//output list for QA histograms

  Bool_t useCentrality;// use centrality as event class estimator (default = OFF)
  Bool_t fUseMultSelectionFramework;// use the AliMultSelection framework; default: kFALSE
  Bool_t fUseUncheckedCentrality; // use unchecked centrality; default: kFALSE
  Bool_t fillOnlySecondaries;//fill only secondary particles (only for MC running)
  Bool_t fillHFVertexingTracks;//fill HF vertexing tracks

  TString fHFBranchName;// name of HF branch: D0toKpi, Dstar, ...

  Int_t fBitIgnore1;//filter bit to be skipped when looping over tracks (-1 = OFF)
  Int_t fBitIgnore2;//filter bit to be skipped when looping over tracks (-1 = OFF)

  Double_t fCentralityPercentileMin;// minimum centrality threshold (default = 0)
  Double_t fCentralityPercentileMax;// maximum centrality threshold (default = 80)
  Double_t fPtMin;// minimum pT threshold (default = 0)
  Double_t fPtMax;// maximum pT threshold (default = 1000)
  Double_t fEtaMin;// minimum eta threshold (default = -10)
  Double_t fEtaMax;// maximum eta threshold (default = 10)

  TH2D* fHistEventStats;//QA histogram for event trigger bit statistics vs. centrality
  TH2D* fHistTrackStats;//QA histogram for track filter bit statistics vs. centrality
  TH3D* fHistKinematics[gNCharge][gBitMax];//QA histograms for kinematics (eta, phi, pT) for different filter bits
  TH2D* fHistDCAconstrained[gNCharge][gBitMax];//QA histograms for DCA (xy,z) for different filter bits for constrained tracks (stored in DCA methods)
  TH3D* fHistDCAglobal[gNCharge][gBitMax];//QA histograms for DCA (xy,z) for different filter bits for global tracks (stored in Position methods)
  TH2D* fHistChiClus[gNCharge][gBitMax];//QA histograms for Chi2 and number of TPC clusters for different filter bits


  
  ClassDef(AliAnalysisTaskAODFilterBitQA, 0); //
};



#endif
