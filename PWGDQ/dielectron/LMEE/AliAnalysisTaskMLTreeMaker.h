

#ifndef ALIANALYSISTASKMLTREEMAKER_H
#define ALIANALYSISTASKMLTREEMAKER_H
class TH1F;
class TList;
class TH2D;
class TH3D;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"
#ifndef ALIANALYSISTASKSE_H
#endif

// Authors: Sebastian Lehner (SMI Vienna) - selehner@cern.ch


class AliAnalysisTaskMLTreeMaker : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMLTreeMaker(const char *name);
  AliAnalysisTaskMLTreeMaker();
  virtual ~AliAnalysisTaskMLTreeMaker(){} 
   
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);
//~ 

  void SetCentralityPercentileRange(Double_t min, Double_t max){
    fCentralityPercentileMin = min;
    fCentralityPercentileMax = max;

    Printf("Thresholds Set");
    Printf("pT = %f - %f",fPtMin,fPtMax);
    Printf("eta = %f - %f",fEtaMin,fEtaMax);
    Printf("cent = %f - %f",fCentralityPercentileMin,fCentralityPercentileMax);
  }

  void SetPtRange(Double_t min, Double_t max){
    fPtMin = min;
    fPtMax = max;

    Printf("Thresholds Set");
    Printf("pT = %f - %f",fPtMin,fPtMax);
    Printf("eta = %f - %f",fEtaMin,fEtaMax);
    Printf("cent = %f - %f",fCentralityPercentileMin,fCentralityPercentileMax);
  }

  void SetEtaRange(Double_t min, Double_t max){
    fEtaMin = min;
    fEtaMax = max;

    Printf("Thresholds Set");
    Printf("pT = %f - %f",fPtMin,fPtMax);
    Printf("eta = %f - %f",fEtaMin,fEtaMax);
    Printf("cent = %f - %f",fCentralityPercentileMin,fCentralityPercentileMax);
  }

  void SetFilterBit(Int_t filterBit){
    fFilterBit = filterBit;
  }
  
    void SetLoCuts(Bool_t x){       //use loose cuts?
    loCuts = x;
  }

  void SetESigRangeITS(Double_t min, Double_t max){
    fESigITSMin = min;
    fESigITSMax = max;
    Printf("ITS electron nSigma values (inclusive)");
    Printf(" Min = %f, Max = %f", fESigITSMin, fESigITSMax);
  }

  void SetESigRangeTPC(Double_t min, Double_t max){
    fESigTPCMin = min;
    fESigTPCMax = max;
    Printf("TPC electron nSigma values (inclusive)");
    Printf(" Min = %f, Max = %f", fESigTPCMin, fESigTPCMax);
  }

  void SetESigRangeTOF(Double_t min, Double_t max){
    fESigTOFMin = min;
    fESigTOFMax = max;
    Printf("TOF electron nSigma values (inclusive)");
    Printf(" Min = %f, Max = %f", fESigTOFMin, fESigTOFMax);
  }

  
  void SetPSigRangeTPC(Double_t min, Double_t max){
    fPSigTPCMin = min;
    fPSigTPCMax = max;
    Printf("TPC pion nSigma values (exclusive)");
    Printf(" Min = %f, Max = %f", fPSigTPCMin, fPSigTPCMax);
  }
  
  void TPCPionPID( Bool_t answer = kTRUE){ 
    fUsePionPIDTPC = answer;
    Printf("No pion PID in TPC applied");
  }

  void StorePionSigmaValues( Bool_t answer = kFALSE){
    fPionSigmas = answer;
    Printf("Pion nSigma values for ITS, TPC and TOF will be written to tree");
  }

  void StoreKaonSigmaValues( Bool_t answer = kFALSE){
    fKaonSigmas = answer;
    Printf("Kaon nSigma values for ITS, TPC and TOF will be written to tree");
  }
 private:
 
  AliPIDResponse *fPIDResponse;     //! PID response object

  //AliPIDCombined *fPIDCombined;    
     
  std::vector<Double_t> eta;
  std::vector<Double_t> phi;
  std::vector<Double_t> pt;
  std::vector<Int_t> charge;
  std::vector<Int_t> enh;  

//  std::vector<Int_t> NClustersITS;
  std::vector<Float_t> NCrossedRowsTPC;
  std::vector<Int_t> NClustersTPC;
  std::vector<Bool_t> HasSPDfirstHit; 
  std::vector<Double_t> RatioCrossedRowsFindableClusters;  
  std::vector<Int_t> NTPCSignal; 
  
  Bool_t loCuts;        //loose cuts?
  
//  std::vector<Int_t> IsBG;
  
  Int_t runn;
  Int_t n;
  Double_t cent;
  
  Double_t IsEventAccepted(AliVEvent *event);
  Int_t GetAcceptedTracks(AliVEvent *event, Double_t gCentrality);
  Bool_t GetDCA(const AliVEvent* event, const AliAODTrack *track, Double_t* d0z0, Double_t* covd0z0);
  
  AliAnalysisTaskMLTreeMaker(const AliAnalysisTaskMLTreeMaker&); // not implemented

  AliAnalysisTaskMLTreeMaker& operator=(const AliAnalysisTaskMLTreeMaker&); // not implemented


  TList *fList;//output list for QA histograms

  Double_t fCentralityPercentileMin;// minimum centrality threshold (default = 0)
  Double_t fCentralityPercentileMax;// maximum centrality threshold (default = 80)

  Double_t fPtMin;// minimum pT threshold (default = 0)
  Double_t fPtMax;// maximum pT threshold (default = 1000)
  Double_t fEtaMin;// minimum eta threshold (default = -10)
  Double_t fEtaMax;// maximum eta threshold (default = 10)

  //Electron nSigma values, initilaised to common cut values but can be set manually 
  Double_t fESigITSMin;
  Double_t fESigITSMax;
  Double_t fESigTPCMin;
  Double_t fESigTPCMax;
  Double_t fESigTOFMin;
  Double_t fESigTOFMax;
  //Pion nSigma values, initilaised to common cut values but can be set manually 
  Double_t fPSigTPCMin; 
  Double_t fPSigTPCMax; 

  Bool_t fUsePionPIDTPC; //Use Pion nSigma information in TPC

  //Flags to write extra particle PID information
  //Initiliased to kFALSE
  Bool_t fPionSigmas; 
  Bool_t fKaonSigmas;

  Int_t fFilterBit;// track cut bit from track selection (default = 96)

  AliESDtrackCuts* fESDTrackCuts;
  
  Int_t gMultiplicity;
  Int_t mcTrackIndex;
  AliMCEvent* fMcArray; 

  std::vector<Double_t> EsigTPC;
  std::vector<Double_t> EsigTOF;
  std::vector<Double_t> EsigITS;
  
  std::vector<Double_t> PsigTPC;
  std::vector<Double_t> PsigTOF;
  std::vector<Double_t> PsigITS;

  std::vector<Double_t> KsigTPC;
  std::vector<Double_t> KsigTOF;
  std::vector<Double_t> KsigITS;

  Bool_t hasMC;
  Bool_t Rej;
 
  std::vector<Double_t> MCpt;
  std::vector<Double_t> MCeta;
  std::vector<Double_t> MCphi;
  
  std::vector<Float_t> dcar;    //DCA
  std::vector<Float_t> dcaz;
  
  Double_t vertx;
  Double_t verty;
  Double_t vertz; 
  
  std::vector<Int_t> nITS;
  std::vector<Double_t> nITSshared;
  std::vector<Double_t> chi2ITS;
  std::vector<Double_t> chi2TPC;
  std::vector<Double_t> chi2Global;
  std::vector<Double_t> chi2GlobalvsTPC;
  Int_t	fCutMaxChi2TPCConstrainedVsGlobalVertexType;
  
  std::vector<Int_t> pdg;
  std::vector<Int_t> pdgmother;
  std::vector<Int_t> hasmother;
  std::vector<Int_t> motherlabel;
  
  //TBits*            fUsedVars;                // used variables by AliDielectronVarManager
  
//  Double_t probs[AliPID::kSPECIESC];

  TTree* fTree;
  TH1F* fQAHist;
  
//  TH2D* fHistTrackStats;//QA histogram for track filter bit statistics vs. centrality

//  TH3D* fHistEtaPhiPt;//QA histogram for eta/phi/pt distribution


  ClassDef(AliAnalysisTaskMLTreeMaker, 1); //

};



#endif

