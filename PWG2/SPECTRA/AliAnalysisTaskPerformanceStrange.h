#ifndef ALIANALYSISTASKPERFORMANCESTRANGE_H
#define ALIANALYSISTASKPERFORMANCESTRANGE_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//	        AliAnalysisTaskPerformanceSrange class
//    This task is for a performance study of V0 identification.
//                It works with MC info and ESD tree.
//                 Author: H.Ricaud, H.Ricaud@gsi.de
//-----------------------------------------------------------------

class TString;
class TList;
class TH1F;
class TH2F;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPerformanceStrange : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPerformanceStrange();
  AliAnalysisTaskPerformanceStrange(const char *name);
  virtual ~AliAnalysisTaskPerformanceStrange() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
  void   SetCollidingSystems(Int_t collidingSystems = 0) {fCollidingSystems = collidingSystems;}
  void   SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  void   SetUsePID(const char* usePID) {fUsePID = usePID;}
  void   SetAnalysisCut(const char* useCut) {fUseCut = useCut;}
  Double_t MyRapidity(Double_t rE, Double_t rPz) const;
 
 private:
  TString      fAnalysisType;                   //  "ESD" or "AOD"
  Int_t        fCollidingSystems;               //  Colliding systems 0/1 for pp/PbPb  
  TString      fUsePID;                         //  "withPID" or "noPID"
  TString      fUseCut;                         //  "yes" or "no"

  TList       *fListHist;			//! Output List

  // MC histograms
  TH1F        *fHistMCMultiplicityPrimary;       //! Histo
  TH1F        *fHistMCMultiplicityTracks;       //! Histo
  
  TH2F        *fHistMCtracksProdRadiusK0s;       //! Histo
  TH2F        *fHistMCtracksProdRadiusLambda;       //! Histo
  TH2F        *fHistMCtracksProdRadiusAntiLambda;       //! Histo

  TH1F        *fHistMCtracksDecayRadiusK0s;       //! Histo
  TH1F        *fHistMCtracksDecayRadiusLambda;       //! Histo
  TH1F        *fHistMCtracksDecayRadiusAntiLambda;       //! Histo

  TH1F        *fHistMCPtAllK0s;       //! Histo
  TH1F        *fHistMCPtAllLambda;       //! Histo
  TH1F        *fHistMCPtAllAntiLambda;       //! Histo

  TH1F        *fHistMCProdRadiusK0s;       //! Histo
  TH1F        *fHistMCProdRadiusLambda;       //! Histo
  TH1F        *fHistMCProdRadiusAntiLambda;       //! Histo

  TH2F        *fHistMCPtVsYK0s;       //! Histo
  TH2F        *fHistMCPtVsYLambda;       //! Histo
  TH2F        *fHistMCPtVsYAntiLambda;       //! Histo

  TH1F        *fHistMCPtLambdaFromSigma;       //! Histo
  TH1F        *fHistMCPtAntiLambdaFromSigma;       //! Histo

  TH1F        *fHistNTimesRecK0s;       //! Histo
  TH1F        *fHistNTimesRecK0sMI;       //! Histo
  TH1F        *fHistNTimesRecLambda;       //! Histo
  TH1F        *fHistNTimesRecLambdaMI;       //! Histo
  TH1F        *fHistNTimesRecAntiLambda;       //! Histo
  TH1F        *fHistNTimesRecAntiLambdaMI;       //! Histo

  TH2F        *fHistNTimesRecK0sVsPt;       //! Histo
  TH2F        *fHistNTimesRecK0sVsPtMI;       //! Histo
  TH2F        *fHistNTimesRecLambdaVsPt;       //! Histo
  TH2F        *fHistNTimesRecLambdaVsPtMI;       //! Histo
  TH2F        *fHistNTimesRecAntiLambdaVsPt;       //! Histo
  TH2F        *fHistNTimesRecAntiLambdaVsPtMI;       //! Histo

  // ESD histograms
  TH1F        *fHistTrackPerEvent;       //! Histo
  TH1F        *fHistMCDaughterTrack;       //! Histo

  TH1F        *fHistPrimaryVertexX;       //! Histo
  TH1F        *fHistPrimaryVertexY;       //! Histo
  TH1F        *fHistPrimaryVertexZ;       //! Histo

  TH2F        *fHistDcaPosToPrimVertex;       //! Histo
  TH2F        *fHistDcaNegToPrimVertex;       //! Histo
  TH2F        *fHistDcaPosToPrimVertexZoom;       //! Histo
  TH2F        *fHistDcaNegToPrimVertexZoom;       //! Histo
  TH2F        *fHistRadiusV0;       //! Histo
  TH2F        *fHistDecayLengthV0;       //! Histo
  TH2F        *fHistDcaV0Daughters;       //! Histo
  TH2F        *fHistChi2;       //! Histo
  TH2F        *fHistCosPointAngle;       //! Histo
  TH2F        *fHistCosPointAngleZoom;       //! Histo
  TH2F        *fHistProdRadius;       //! Histo
  TH2F        *fHistProdRadiusMI;       //! Histo

  TH1F        *fHistV0Multiplicity;  //! Histo
  TH1F        *fHistV0MultiplicityMI; //! Histo

  TH2F        *fHistPtVsYK0s;       //! Histo
  TH2F        *fHistPtVsYK0sMI;       //! Histo
  TH2F        *fHistPtVsYLambda;       //! Histo
  TH2F        *fHistPtVsYLambdaMI;       //! Histo
  TH2F        *fHistPtVsYAntiLambda;       //! Histo
  TH2F        *fHistPtVsYAntiLambdaMI;       //! Histo

  TH1F        *fHistMassK0;       //! Histo
  TH1F        *fHistMassK0MI;       //! Histo
  TH1F        *fHistMassLambda;       //! Histo
  TH1F        *fHistMassLambdaMI;       //! Histo
  TH1F        *fHistMassAntiLambda;       //! Histo
  TH1F        *fHistMassAntiLambdaMI;       //! Histo

  TH2F        *fHistMassVsRadiusK0;       //! Histo
  TH2F        *fHistMassVsRadiusK0MI;       //! Histo
  TH2F        *fHistMassVsRadiusLambda;       //! Histo
  TH2F        *fHistMassVsRadiusLambdaMI;       //! Histo
  TH2F        *fHistMassVsRadiusAntiLambda;       //! Histo
  TH2F        *fHistMassVsRadiusAntiLambdaMI;       //! Histo

  TH2F        *fHistPtVsMassK0;       //! Histo
  TH2F        *fHistPtVsMassK0MI;       //! Histo
  TH2F        *fHistPtVsMassLambda;       //! Histo
  TH2F        *fHistPtVsMassLambdaMI;       //! Histo
  TH2F        *fHistPtVsMassAntiLambda;       //! Histo
  TH2F        *fHistPtVsMassAntiLambdaMI;       //! Histo

  TH2F        *fHistArmenterosPodolanski;       //! Histo
  TH2F        *fHistArmenterosPodolanskiMI;       //! Histo

  //PID
  TH1F        *fHistNsigmaPosPionAntiLambda;    //! Histo
  TH1F        *fHistNsigmaNegProtonAntiLambda;   //! Histo
  TH1F        *fHistNsigmaPosProtonLambda;        //! Histo
  TH1F        *fHistNsigmaNegPionLambda;           //! Histo
  TH1F        *fHistNsigmaPosPionK0;                //! Histo
  TH1F        *fHistNsigmaNegPionK0;                 //! Histo

  // Associated particles histograms
  TH1F        *fHistAsMcPtK0;       //! Histo
  TH1F        *fHistAsMcPtK0MI;       //! Histo
  TH1F        *fHistAsMcPtLambda;       //! Histo
  TH1F        *fHistAsMcPtLambdaMI;       //! Histo
  TH1F        *fHistAsMcPtAntiLambda;       //! Histo
  TH1F        *fHistAsMcPtAntiLambdaMI;       //! Histo
  TH1F        *fHistAsMcPtZoomK0;       //! Histo
  TH1F        *fHistAsMcPtZoomK0MI;       //! Histo
  TH1F        *fHistAsMcPtZoomLambda;       //! Histo
  TH1F        *fHistAsMcPtZoomLambdaMI;       //! Histo

  TH1F        *fHistAsMcProdRadiusK0;       //! Histo
  TH1F        *fHistAsMcProdRadiusK0MI;       //! Histo
  TH1F        *fHistAsMcProdRadiusLambda;       //! Histo
  TH1F        *fHistAsMcProdRadiusLambdaMI;       //! Histo
  TH1F        *fHistAsMcProdRadiusAntiLambda;       //! Histo
  TH1F        *fHistAsMcProdRadiusAntiLambdaMI;       //! Histo

  TH2F        *fHistAsMcProdRadiusXvsYK0s;       //! Histo
  TH2F        *fHistAsMcProdRadiusXvsYK0sMI;       //! Histo
  TH2F        *fHistAsMcProdRadiusXvsYLambda;       //! Histo
  TH2F        *fHistAsMcProdRadiusXvsYLambdaMI;       //! Histo
  TH2F        *fHistAsMcProdRadiusXvsYAntiLambda;       //! Histo
  TH2F        *fHistAsMcProdRadiusXvsYAntiLambdaMI;       //! Histo

  TH1F        *fHistPidMcMassK0;       //! Histo
  TH1F        *fHistPidMcMassK0MI;       //! Histo
  TH1F        *fHistPidMcMassLambda;       //! Histo
  TH1F        *fHistPidMcMassLambdaMI;       //! Histo
  TH1F        *fHistPidMcMassAntiLambda;       //! Histo
  TH1F        *fHistPidMcMassAntiLambdaMI;       //! Histo
  TH1F        *fHistAsMcMassK0;       //! Histo
  TH1F        *fHistAsMcMassK0MI;       //! Histo
  TH1F        *fHistAsMcMassLambda;       //! Histo
  TH1F        *fHistAsMcMassLambdaMI;       //! Histo
  TH1F        *fHistAsMcMassAntiLambda;       //! Histo
  TH1F        *fHistAsMcMassAntiLambdaMI;       //! Histo

  TH2F        *fHistAsMcPtVsMassK0;       //! Histo
  TH2F        *fHistAsMcPtVsMassK0MI;       //! Histo
  TH2F        *fHistAsMcPtVsMassLambda;       //! Histo
  TH2F        *fHistAsMcPtVsMassLambdaMI;       //! Histo
  TH2F        *fHistAsMcPtVsMassAntiLambda;       //! Histo
  TH2F        *fHistAsMcPtVsMassAntiLambdaMI;       //! Histo

  TH2F        *fHistAsMcMassVsRadiusK0;       //! Histo
  TH2F        *fHistAsMcMassVsRadiusK0MI;       //! Histo
  TH2F        *fHistAsMcMassVsRadiusLambda;       //! Histo
  TH2F        *fHistAsMcMassVsRadiusLambdaMI;       //! Histo
  TH2F        *fHistAsMcMassVsRadiusAntiLambda;       //! Histo
  TH2F        *fHistAsMcMassVsRadiusAntiLambdaMI;       //! Histo

  TH1F        *fHistAsMcResxK0;       //! Histo
  TH1F        *fHistAsMcResyK0;       //! Histo
  TH1F        *fHistAsMcReszK0;       //! Histo

  TH2F        *fHistAsMcResrVsRadiusK0;       //! Histo
  TH2F        *fHistAsMcReszVsRadiusK0;       //! Histo

  TH1F        *fHistAsMcResxK0MI;       //! Histo
  TH1F        *fHistAsMcResyK0MI;       //! Histo
  TH1F        *fHistAsMcReszK0MI;       //! Histo

  TH2F        *fHistAsMcResrVsRadiusK0MI;       //! Histo
  TH2F        *fHistAsMcReszVsRadiusK0MI;       //! Histo

  TH1F        *fHistAsMcResxLambda;       //! Histo
  TH1F        *fHistAsMcResyLambda;       //! Histo
  TH1F        *fHistAsMcReszLambda;       //! Histo

  TH2F        *fHistAsMcResrVsRadiusLambda;       //! Histo
  TH2F        *fHistAsMcReszVsRadiusLambda;       //! Histo
    
  TH1F        *fHistAsMcResxLambdaMI;       //! Histo
  TH1F        *fHistAsMcResyLambdaMI;       //! Histo
  TH1F        *fHistAsMcReszLambdaMI;       //! Histo

  TH2F        *fHistAsMcResrVsRadiusLambdaMI;       //! Histo
  TH2F        *fHistAsMcReszVsRadiusLambdaMI;       //! Histo

  TH1F        *fHistAsMcResxAntiLambda;       //! Histo
  TH1F        *fHistAsMcResyAntiLambda;       //! Histo
  TH1F        *fHistAsMcReszAntiLambda;       //! Histo

  TH2F        *fHistAsMcResrVsRadiusAntiLambda;       //! Histo
  TH2F        *fHistAsMcReszVsRadiusAntiLambda;       //! Histo
    
  TH1F        *fHistAsMcResxAntiLambdaMI;       //! Histo
  TH1F        *fHistAsMcResyAntiLambdaMI;       //! Histo
  TH1F        *fHistAsMcReszAntiLambdaMI;       //! Histo

  TH2F        *fHistAsMcResrVsRadiusAntiLambdaMI;       //! Histo
  TH2F        *fHistAsMcReszVsRadiusAntiLambdaMI;       //! Histo

  TH1F        *fHistAsMcResPtK0;       //! Histo
  TH1F        *fHistAsMcResPtK0MI;       //! Histo
  TH1F        *fHistAsMcResPtLambda;       //! Histo
  TH1F        *fHistAsMcResPtLambdaMI;       //! Histo
  TH1F        *fHistAsMcResPtAntiLambda;       //! Histo
  TH1F        *fHistAsMcResPtAntiLambdaMI;       //! Histo

  TH2F        *fHistAsMcResPtVsRapK0;       //! Histo
  TH2F        *fHistAsMcResPtVsRapK0MI;       //! Histo
  TH2F        *fHistAsMcResPtVsRapLambda;       //! Histo
  TH2F        *fHistAsMcResPtVsRapLambdaMI;       //! Histo
  TH2F        *fHistAsMcResPtVsRapAntiLambda;       //! Histo
  TH2F        *fHistAsMcResPtVsRapAntiLambdaMI;       //! Histo
  TH2F        *fHistAsMcResPtVsPtK0;       //! Histo
  TH2F        *fHistAsMcResPtVsPtK0MI;       //! Histo
  TH2F        *fHistAsMcResPtVsPtLambda;       //! Histo
  TH2F        *fHistAsMcResPtVsPtLambdaMI;       //! Histo
  TH2F        *fHistAsMcResPtVsPtAntiLambda;       //! Histo
  TH2F        *fHistAsMcResPtVsPtAntiLambdaMI;       //! Histo
  

  TH1F        *fHistAsMcMotherPdgCodeK0s;       //! Histo
  TH1F        *fHistAsMcMotherPdgCodeK0sMI;       //! Histo
  TH1F        *fHistAsMcMotherPdgCodeLambda;       //! Histo
  TH1F        *fHistAsMcMotherPdgCodeLambdaMI;       //! Histo
  TH1F        *fHistAsMcMotherPdgCodeAntiLambda;       //! Histo
  TH1F        *fHistAsMcMotherPdgCodeAntiLambdaMI;       //! Histo

  TH1F        *fHistAsMcPtLambdaFromSigma;       //! Histo
  TH1F        *fHistAsMcPtLambdaFromSigmaMI;       //! Histo
  TH1F        *fHistAsMcPtAntiLambdaFromSigma;       //! Histo
  TH1F        *fHistAsMcPtAntiLambdaFromSigmaMI;       //! Histo

  // Associated secondary particles:
  TH2F        *fHistAsMcSecondaryPtVsYK0s;       //! Histo
  TH2F        *fHistAsMcSecondaryPtVsYK0sMI;       //! Histo
  TH2F        *fHistAsMcSecondaryPtVsYLambda;       //! Histo
  TH2F        *fHistAsMcSecondaryPtVsYLambdaMI;       //! Histo
  TH2F        *fHistAsMcSecondaryPtVsYAntiLambda;       //! Histo
  TH2F        *fHistAsMcSecondaryPtVsYAntiLambdaMI;       //! Histo

  TH1F        *fHistAsMcSecondaryProdRadiusK0s;       //! Histo
  TH1F        *fHistAsMcSecondaryProdRadiusK0sMI;       //! Histo
  TH1F        *fHistAsMcSecondaryProdRadiusLambda;       //! Histo
  TH1F        *fHistAsMcSecondaryProdRadiusLambdaMI;       //! Histo
  TH1F        *fHistAsMcSecondaryProdRadiusAntiLambda;       //! Histo
  TH1F        *fHistAsMcSecondaryProdRadiusAntiLambdaMI;       //! Histo

  TH2F        *fHistAsMcSecondaryProdRadiusXvsYK0s;       //! Histo
  TH2F        *fHistAsMcSecondaryProdRadiusXvsYK0sMI;       //! Histo
  TH2F        *fHistAsMcSecondaryProdRadiusXvsYLambda;       //! Histo
  TH2F        *fHistAsMcSecondaryProdRadiusXvsYLambdaMI;       //! Histo
  TH2F        *fHistAsMcSecondaryProdRadiusXvsYAntiLambda;       //! Histo
  TH2F        *fHistAsMcSecondaryProdRadiusXvsYAntiLambdaMI;       //! Histo

  TH1F        *fHistAsMcSecondaryMotherPdgCodeK0s;       //! Histo
  TH1F        *fHistAsMcSecondaryMotherPdgCodeK0sMI;       //! Histo
  TH1F        *fHistAsMcSecondaryMotherPdgCodeLambda;       //! Histo
  TH1F        *fHistAsMcSecondaryMotherPdgCodeLambdaMI;       //! Histo
  TH1F        *fHistAsMcSecondaryMotherPdgCodeAntiLambda;       //! Histo
  TH1F        *fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI;       //! Histo

  TH1F        *fHistAsMcSecondaryPtLambdaFromSigma;       //! Histo
  TH1F        *fHistAsMcSecondaryPtLambdaFromSigmaMI;       //! Histo
  TH1F        *fHistAsMcSecondaryPtAntiLambdaFromSigma;       //! Histo
  TH1F        *fHistAsMcSecondaryPtAntiLambdaFromSigmaMI;       //! Histo

  AliAnalysisTaskPerformanceStrange(const AliAnalysisTaskPerformanceStrange&); 
  AliAnalysisTaskPerformanceStrange& operator=(const AliAnalysisTaskPerformanceStrange&); 

  ClassDef(AliAnalysisTaskPerformanceStrange, 1); 
};

#endif
