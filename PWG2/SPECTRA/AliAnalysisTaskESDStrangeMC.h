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
  TH1F        *fHistMCMultiplicity;

  TH2F        *fHistMCPtVsYK0s;
  TH2F        *fHistMCPtVsYLambda;
  TH2F        *fHistMCPtVsYAntiLambda;

  // ESD histograms
  TH1F        *fHistTrackPerEvent;
  TH1F        *fHistMCDaughterTrack;

  TH1F        *fHistPrimaryVertexX;
  TH1F        *fHistPrimaryVertexY;
  TH1F        *fHistPrimaryVertexZ;

  TH2F        *fHistDcaPosToPrimVertex;
  TH2F        *fHistDcaNegToPrimVertex;
  TH2F        *fHistDcaPosToPrimVertexZoom;
  TH2F        *fHistDcaNegToPrimVertexZoom;
  TH2F        *fHistRadiusV0;
  TH2F        *fHistDecayLengthV0;
  TH2F        *fHistDcaV0Daughters;
  TH2F        *fHistChi2;
  TH2F        *fHistCosPointAngle;
  TH2F        *fHistCosPointAngleZoom;
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

  TH2F        *fHistMassVsRadiusK0;
  TH2F        *fHistMassVsRadiusK0MI;
  TH2F        *fHistMassVsRadiusLambda;
  TH2F        *fHistMassVsRadiusLambdaMI;
  TH2F        *fHistMassVsRadiusAntiLambda;
  TH2F        *fHistMassVsRadiusAntiLambdaMI;

  TH2F        *fHistArmenterosPodolanski;
  TH2F        *fHistArmenterosPodolanskiMI;


  // Associated particles histograms
  TH1F        *fHistAsMcPtK0;
  TH1F        *fHistAsMcPtK0MI;
  TH1F        *fHistAsMcPtLambda;
  TH1F        *fHistAsMcPtLambdaMI;
  TH1F        *fHistAsMcPtAntiLambda;
  TH1F        *fHistAsMcPtAntiLambdaMI;
  TH1F        *fHistAsMcPtZoomK0;
  TH1F        *fHistAsMcPtZoomK0MI;
  TH1F        *fHistAsMcPtZoomLambda;
  TH1F        *fHistAsMcPtZoomLambdaMI;
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

  TH2F        *fHistAsMcMassVsRadiusK0;
  TH2F        *fHistAsMcMassVsRadiusK0MI;
  TH2F        *fHistAsMcMassVsRadiusLambda;
  TH2F        *fHistAsMcMassVsRadiusLambdaMI;
  TH2F        *fHistAsMcMassVsRadiusAntiLambda;
  TH2F        *fHistAsMcMassVsRadiusAntiLambdaMI;

  TH1F        *fHistAsMcResxK0;
  TH1F        *fHistAsMcResyK0;
  TH1F        *fHistAsMcReszK0;

  TH2F        *fHistAsMcResrVsRadiusK0;
  TH2F        *fHistAsMcReszVsRadiusK0;

  TH1F        *fHistAsMcResxK0MI;
  TH1F        *fHistAsMcResyK0MI;
  TH1F        *fHistAsMcReszK0MI;

  TH2F        *fHistAsMcResrVsRadiusK0MI;
  TH2F        *fHistAsMcReszVsRadiusK0MI;

  TH1F        *fHistAsMcResxLambda;
  TH1F        *fHistAsMcResyLambda;
  TH1F        *fHistAsMcReszLambda;

  TH2F        *fHistAsMcResrVsRadiusLambda;
  TH2F        *fHistAsMcReszVsRadiusLambda;
    
  TH1F        *fHistAsMcResxLambdaMI;
  TH1F        *fHistAsMcResyLambdaMI;
  TH1F        *fHistAsMcReszLambdaMI;

  TH2F        *fHistAsMcResrVsRadiusLambdaMI;
  TH2F        *fHistAsMcReszVsRadiusLambdaMI;

  TH1F        *fHistAsMcResxAntiLambda;
  TH1F        *fHistAsMcResyAntiLambda;
  TH1F        *fHistAsMcReszAntiLambda;

  TH2F        *fHistAsMcResrVsRadiusAntiLambda;
  TH2F        *fHistAsMcReszVsRadiusAntiLambda;
    
  TH1F        *fHistAsMcResxAntiLambdaMI;
  TH1F        *fHistAsMcResyAntiLambdaMI;
  TH1F        *fHistAsMcReszAntiLambdaMI;

  TH2F        *fHistAsMcResrVsRadiusAntiLambdaMI;
  TH2F        *fHistAsMcReszVsRadiusAntiLambdaMI;


  AliAnalysisTaskESDStrangeMC(const AliAnalysisTaskESDStrangeMC&); 
  AliAnalysisTaskESDStrangeMC& operator=(const AliAnalysisTaskESDStrangeMC&); 

  ClassDef(AliAnalysisTaskESDStrangeMC, 1); 
};

#endif
