

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


 private:
     
 
  AliPIDResponse *fPIDResponse;     //! PID response object

  //AliPIDCombined *fPIDCombined;    

     
  vector<Double_t> eta;
  vector<Double_t> phi;
  vector<Double_t> pt;
  vector<Int_t> charge;
  
  vector<Int_t> IsBG;
  vector<Int_t> runn;
  
  Double_t pttemp;
  Double_t etatemp;
  
  Int_t n;
  
  Double_t IsEventAccepted(AliESDEvent *event);
  Int_t GetAcceptedTracks(AliVEvent *event, Double_t gCentrality);
  
 AliAnalysisTaskMLTreeMaker(const AliAnalysisTaskMLTreeMaker&); // not implemented

 AliAnalysisTaskMLTreeMaker& operator=(const AliAnalysisTaskMLTreeMaker&); // not implemented


  TList *fList;//output list for QA histograms


  Double_t fCentralityPercentileMin;// minimum centrality threshold (default = 0)

  Double_t fCentralityPercentileMax;// maximum centrality threshold (default = 80)

  Double_t fPtMin;// minimum pT threshold (default = 0)

  Double_t fPtMax;// maximum pT threshold (default = 1000)

  Double_t fEtaMin;// minimum eta threshold (default = -10)

  Double_t fEtaMax;// maximum eta threshold (default = 10)


  Int_t fFilterBit;// track cut bit from track selection (default = 96)

  AliESDtrackCuts* fESDTrackCuts;
  
  Int_t gMultiplicity;
  
  AliMCEvent* fMcArray; 

  vector<Double_t> EsigTPC;
  vector<Double_t> EsigTOF;
  vector<Double_t> EsigITS;
  
  
 Double_t tempEsigTPC;
 Double_t tempEsigTPC2;
 Double_t tempEsigTOF;
 Double_t tempEsigITS;

  Bool_t hasMC;
  
  vector<Double_t> MCpt;
  vector<Double_t> MCeta;
  vector<Double_t> MCphi;
  
  vector<Float_t> dcar;    //DCA

  vector<Float_t> dcaz;
  
  Float_t tempdca[2];

  Double_t vertx;
  Double_t verty;
  Double_t vertz; 
  
  Double_t vert[3];

  
  vector<Int_t> nITS;
  Int_t tempnits;
  Double_t nitssharedtemp;
  vector<Double_t> nITSshared;
  vector<Double_t> chi2ITS;
  vector<Double_t> chi2TPC;
  vector<Double_t> chi2Global;
  vector<Double_t> chi2GlobalvsTPC;
  Int_t	fCutMaxChi2TPCConstrainedVsGlobalVertexType;
  
  vector<Int_t> pdg;
  vector<Int_t> pdgmother;
  vector<Int_t> hasmother;
  vector<Int_t> motherlabel;
  
  
  
  
  //TBits*            fUsedVars;                // used variables by AliDielectronVarManager

  
//  Double_t probs[AliPID::kSPECIESC];

  TTree* fTree;
  TH1F* fQAHist;
//  TH2D* fHistTrackStats;//QA histogram for track filter bit statistics vs. centrality

//  TH3D* fHistEtaPhiPt;//QA histogram for eta/phi/pt distribution


  ClassDef(AliAnalysisTaskMLTreeMaker, 1); //

};



#endif

