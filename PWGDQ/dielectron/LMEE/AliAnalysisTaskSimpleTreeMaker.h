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

  //Int_t	fCutMaxChi2TPCConstrainedVsGlobalVertexType;
 
  Double_t fESigTPCMin; 
  Double_t fESigTPCMax; 

  ClassDef(AliAnalysisTaskSimpleTreeMaker, 1); //

};



#endif

