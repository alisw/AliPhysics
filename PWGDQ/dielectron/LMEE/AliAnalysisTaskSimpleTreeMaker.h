#ifndef ALIANALYSISTASKSIMPLETREEMAKER_H
#define ALIANALYSISTASKSIMPLETREEMAKER_H
class TH1F;
class TList;
class TH2D;
class TH3D;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "TFile.h"
#include "TTreeStream.h"
#include "AliAnalysisTaskSE.h"
#ifndef ALIANALYSISTASKSE_H
#endif


class AliAnalysisTaskSimpleTreeMaker : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSimpleTreeMaker(const char *name);
  AliAnalysisTaskSimpleTreeMaker();
  virtual ~AliAnalysisTaskSimpleTreeMaker(){} 
   
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);
//~ 

  void SetCentralityPercentileRange(Double_t min, Double_t max){
    fCentralityPercentileMin = min;
    fCentralityPercentileMax = max;

    Printf("Set: cent = %f - %f",fCentralityPercentileMin,fCentralityPercentileMax);
  }

  void SetPtRange(Double_t min, Double_t max){
    fPtMin = min;
    fPtMax = max;

    Printf("Set: pT = %f - %f",fPtMin,fPtMax);
  }

  void SetEtaRange(Double_t min, Double_t max){
    fEtaMin = min;
    fEtaMax = max;

    Printf("Set: eta = %f - %f",fEtaMin,fEtaMax);
  }

  void SetESigRangeTPC(Double_t min, Double_t max){
    fESigTPCMin = min;
    fESigTPCMax = max;
    Printf("TPC electron nSigma values (inclusive)");
    Printf(" Min = %f, Max = %f", fESigTPCMin, fESigTPCMax);
  }

  void SetMC(Bool_t answer){ fIsMC = answer; }

  void SetIsIonColl(Bool_t answer){ isIonColl = answer; }
    
  //Flags for ITS and TOF PID usage
  void UseITSPID(Bool_t answer){ fPIDcutITS = answer; }
  void UseTOFPID(Bool_t answer){ fPIDcutTOF = answer; }

  //Track cut setters. StandardITSTPC2011 cuts used if nothing specified
  void SetTPCminClusters(Int_t number){
    fESDtrackCuts->SetMinNClustersTPC(number);
  }

  void SetTPCminCrossedRows(Int_t number){
    fESDtrackCuts->SetMinNCrossedRowsTPC(number);
  }

  void SetTPCRatio(Double_t number){
    fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(number);
  }

  void SetTPCChi2PerCluster(Float_t number){
    fESDtrackCuts->SetMaxChi2PerClusterTPC(number);
  }

  void SetITSclusterRequirements(AliESDtrackCuts::Detector detector, AliESDtrackCuts::ITSClusterRequirement requirement){
    fESDtrackCuts->SetClusterRequirementITS(detector, requirement);
  }

  void SetITSChi2PerCluster(Float_t number){
    fESDtrackCuts->SetMaxChi2PerClusterITS(number);
  }

  void SetMaxDCAxy(Float_t number){
    fESDtrackCuts->SetMaxDCAToVertexXY(number);
  }

  void SetMaxDCAPtDep(const char* dist){
    fESDtrackCuts->SetMaxDCAToVertexXYPtDep(dist);
  }

  void SetMaxDCAz(Double_t number){
    fESDtrackCuts->SetMaxDCAToVertexZ(number);
  }

  void SetTPCconstrainedGlobalChi2(Double_t number){
    fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(number);
  }
  
 private:
 
  Double_t IsEventAccepted(AliESDEvent *event);
  //Int_t GetAcceptedTracks(AliVEvent *event, Double_t gCentrality);
  
  AliAnalysisTaskSimpleTreeMaker(const AliAnalysisTaskSimpleTreeMaker&); // not implemented

  AliAnalysisTaskSimpleTreeMaker& operator=(const AliAnalysisTaskSimpleTreeMaker&); // not implemented

  AliESDEvent* esdEvent;  //ESD object
  AliMCEvent* mcEvent;    //MC object

  AliESDtrackCuts *fESDtrackCuts; // ESD track cuts object

  AliPIDResponse *fPIDResponse; //! PID response object

  TTreeStream* fStream;
  TTree* fTree;

  TH1F* fQAhist;
  Double_t fCentralityPercentileMin;// minimum centrality threshold (default = 0)
  Double_t fCentralityPercentileMax;// maximum centrality threshold (default = 80)

  Bool_t fIsMC;
  
  Double_t fPtMin;// minimum pT threshold (default = 0)
  Double_t fPtMax;// maximum pT threshold (default = 1000)
  Double_t fEtaMin;// minimum eta threshold (default = -10)
  Double_t fEtaMax;// maximum eta threshold (default = 10)

 
  Double_t fESigTPCMin; 
  Double_t fESigTPCMax; 

    //Values and flags for PID cuts in ITS and TOF
  Double_t fESigITSMin;
  Double_t fESigITSMax;
  Double_t fESigTOFMin;
  Double_t fESigTOFMax;

  Bool_t fPIDcutITS;
  Bool_t fPIDcutTOF;

  //Values and flag for pion PID cuts in TPC
  Double_t fPSigTPCMin;
  Double_t fPSigTPCMax;

  Bool_t fPionPIDcutTPC;
  Bool_t isIonColl;

  ClassDef(AliAnalysisTaskSimpleTreeMaker, 1); //

};



#endif

