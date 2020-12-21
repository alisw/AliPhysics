#ifndef ALIANALYSISTASKLCTOD0MC_CXX
#define ALIANALYSISTASKLCTOD0MC_CXX

// Analysis task for the  extraction @ the MC level
// Authors: Panos Cristakoglou@cern.ch

class TList;
class TH1F;
class TH2F;

class AliMCEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskLcToD0MC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskLcToD0MC(const char *name = "AliAnalysisTaskLcToD0MC");
  virtual ~AliAnalysisTaskLcToD0MC(); 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);

  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;
    fVyMax = vy;
    fVzMax = vz;
  }

  void SetKinematicsCutsMC(Double_t ptmin, Double_t ptmax, 
			   Double_t etamin, Double_t etamax){
    fPtMin  = ptmin;
    fPtMax  = ptmax;
    fEtaMin = etamin;
    fEtaMax = etamax;

  }

 private:
  TList *fList; //fList object

  TH1F *fHistEventStats; //event stats
  TH1F *fHistTrackStats; //Track filter bit stats
  TH1F *fHistVx; //x coordinate of the primary vertex
  TH1F *fHistVy; //y coordinate of the primary vertex
  TH1F *fHistVz; //z coordinate of the primary vertex

  TH1F *fHistMultiplicity; //multiplicity of accepted tracks
  
  TH2F *fHistPtVsNch; //pt spectrum vs multiplicity

  TH2F *fHistProtonALICE;
  TH2F *fHistPionALICE;
  TH2F *fHistLambdaALICE;
  TH2F *fHistK0sALICE;
  TH2F *fHistLambdacALICE;
  TH2F *fHistDZeroALICE;

  TH2F *fHistProtonCMS;
  TH2F *fHistPionCMS;
  TH2F *fHistLambdaCMS;
  TH2F *fHistK0sCMS;
  TH2F *fHistLambdacCMS;
  TH2F *fHistDZeroCMS;

  TH2F *fHistProtonLHCb;
  TH2F *fHistPionLHCb;
  TH2F *fHistLambdaLHCb;
  TH2F *fHistK0sLHCb;
  TH2F *fHistLambdacLHCb;
  TH2F *fHistDZeroLHCb;

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  Double_t fPtMin;//only used for AODs
  Double_t fPtMax;//only used for AODs
  Double_t fEtaMin;//only used for AODs
  Double_t fEtaMax;//only used for AODs

  Double_t fEtaMinCMS;
  Double_t fEtaMaxCMS;
  Double_t fEtaMinLHCb;
  Double_t fEtaMaxLHCb;

  AliAnalysisTaskLcToD0MC(const AliAnalysisTaskLcToD0MC&); // not implemented
  AliAnalysisTaskLcToD0MC& operator=(const AliAnalysisTaskLcToD0MC&); // not implemented
  
  ClassDef(AliAnalysisTaskLcToD0MC, 1); // example of analysis
};

#endif
