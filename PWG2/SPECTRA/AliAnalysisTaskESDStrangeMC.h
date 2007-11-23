#ifndef AliAnalysisTaskESDStrangeMC_cxx
#define AliAnalysisTaskESDStrangeMC_cxx

// macro to study V0s, with Monte Carlo information access
// Author: H.Ricaud, Helene.Ricaud@IReS.in2p3.fr

class TH1F;
class TH2F;
class TList;
class AliESDEvent;

Double_t myRap(Double_t,Double_t);

#include "AliAnalysisTask.h"

class AliAnalysisTaskESDStrangeMC : public AliAnalysisTask {
 public:
  AliAnalysisTaskESDStrangeMC(const char *name = "AliAnalysisTaskESDStrangeMC");
  virtual ~AliAnalysisTaskESDStrangeMC() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESDEvent *fESD;    //ESD object
  TList       *fListHist;

  // MC histograms
  TH1F        *fHistPtMC; //Pt spectrum

  TH2F        *fHistMCPtVsYK0s;
  TH2F        *fHistMCPtVsYLambda;
  TH2F        *fHistMCPtVsYAntiLambda;

  // ESD histograms
  TH1F        *fHistTrackPerEvent;
  TH1F        *fHistMCDaughterTrack;

  TH2F        *fHistDcaPosToPrimVertex;
  TH2F        *fHistDcaNegToPrimVertex;
  TH2F        *fHistRadiusV0;
  TH2F        *fHistDecayLengthV0;
  TH2F        *fHistDcaV0Daughters;
  TH2F        *fHistChi2;
  TH2F        *fHistPtVsYK0s;
  TH2F        *fHistPtVsYK0sMI;
  TH2F        *fHistPtVsYLambda;
  TH2F        *fHistPtVsYLambdaMI;
  TH2F        *fHistPtVsYAntiLambda;
  TH2F        *fHistPtVsYAntiLambdaMI;

  TH1F        *fHistMassK0;
  TH1F        *fHistMassK0MI;
  TH1F        *fHistMassLambda;
  TH1F        *fHistMassLambdaMI;
  TH1F        *fHistMassAntiLambda;
  TH1F        *fHistMassAntiLambdaMI;


  // Associated particles histograms
  TH1F        *fHistAsMcPtK0;
  TH1F        *fHistAsMcPtK0MI;
  TH1F        *fHistAsMcPtLambda;
  TH1F        *fHistAsMcPtLambdaMI;
  TH1F        *fHistAsMcPtAntiLambda;
  TH1F        *fHistAsMcPtAntiLambdaMI;
  TH1F        *fHistPidMcMassK0;
  TH1F        *fHistPidMcMassK0MI;
  TH1F        *fHistPidMcMassLambda;
  TH1F        *fHistPidMcMassLambdaMI;
  TH1F        *fHistPidMcMassAntiLambda;
  TH1F        *fHistPidMcMassAntiLambdaMI;
  TH1F        *fHistAsMcMassK0;
  TH1F        *fHistAsMcMassK0MI;
  TH1F        *fHistAsMcMassLambda;
  TH1F        *fHistAsMcMassLambdaMI;
  TH1F        *fHistAsMcMassAntiLambda;
  TH1F        *fHistAsMcMassAntiLambdaMI;
   
  AliAnalysisTaskESDStrangeMC(const AliAnalysisTaskESDStrangeMC&); 
  AliAnalysisTaskESDStrangeMC& operator=(const AliAnalysisTaskESDStrangeMC&); 

  ClassDef(AliAnalysisTaskESDStrangeMC, 1); 
};

#endif
