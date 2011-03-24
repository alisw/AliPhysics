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
class AliAnalysisCentralitySelector;
//class TH3F;
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPerformanceStrange : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPerformanceStrange();
  AliAnalysisTaskPerformanceStrange(const char *name);
  virtual ~AliAnalysisTaskPerformanceStrange() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
  void   SetCollidingSystems(Bool_t collidingSystems = 0) {fCollidingSystems = collidingSystems;}
  void   SetAnalysisMC(Bool_t analysisMC) {fAnalysisMC = analysisMC;}
  void   SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  void   SetUsePID(const char* usePID) {fUsePID = usePID;}
  void   SetAnalysisCut(const char* useCut) {fUseCut = useCut;}
  void   UseOnTheFly(Bool_t useOnTheFly) {fUseOnTheFly = useOnTheFly;}
  void   SetCentralitySelector(AliAnalysisCentralitySelector * centr) { fCentrSelector = centr;}
  Double_t MyRapidity(Double_t rE, Double_t rPz) const;
 
 private:
  Double_t fCuts[7];                              //! V0 finding cuts
  Bool_t       fAnalysisMC;                       //  0->No MC or 1->MC analysis
  TString      fAnalysisType;                     //  "ESD" or "AOD"
  Bool_t       fCollidingSystems;                 //  Colliding systems 0/1 for pp/PbPb  
  TString      fUsePID;                           //  "withPID" or "noPID"
  TString      fUseCut;                           //  "yes" or "no"
  Bool_t       fUseOnTheFly;                      //  0->Offline V0s or 1->Onthefly V0s

  AliESDEvent *fESD;                              //! ESD object

  TList       *fListHist;		          //! Output List

  AliAnalysisCentralitySelector * fCentrSelector; // Centrality selector, used to 

  // MC histograms  ---------------------------------------
  TH1F        *fHistMCPrimaryVertexX;             //! Primary Vertex X
  TH1F        *fHistMCPrimaryVertexY;             //! Primary Vertex Y
  TH1F        *fHistMCPrimaryVertexZ;             //! Primary Vertex Z

  TH1F        *fHistMCMultiplicityPrimary;        //! Multiplicity Primaries
  TH1F        *fHistMCMultiplicityTracks;         //! Multiplicity Tracks
  
  TH2F        *fHistMCtracksProdRadiusK0s;        //! Production Radius vs. tracks
  TH2F        *fHistMCtracksProdRadiusLambda;     //! Production Radius vs. tracks
  TH2F        *fHistMCtracksProdRadiusAntiLambda; //! Production Radius vs. tracks

  TH1F        *fHistMCtracksDecayRadiusK0s;       //! Decay Radius
  TH1F        *fHistMCtracksDecayRadiusLambda;    //! Decay Radius
  TH1F        *fHistMCtracksDecayRadiusAntiLambda;//! Decay Radius

  TH1F        *fHistMCPtAllK0s;                   //! Pt Distribution
  TH1F        *fHistMCPtAllLambda;                //! Pt Distribution
  TH1F        *fHistMCPtAllAntiLambda;            //! Pt Distribution

  TH1F        *fHistMCProdRadiusK0s;              //! Production Radius
  TH1F        *fHistMCProdRadiusLambda;           //! Production Radius
  TH1F        *fHistMCProdRadiusAntiLambda;       //! Production Radius

  TH1F        *fHistMCRapK0s;                     //! Rapidity distribution
  TH1F        *fHistMCRapInPtRangeK0s;            //! Rapidity distribution
  TH1F        *fHistMCRapLambda;                  //! Rapidity distribution
  TH1F        *fHistMCRapInPtRangeLambda;         //! Rapidity distribution
  TH1F        *fHistMCRapAntiLambda;              //! Rapidity distribution
  TH1F        *fHistMCRapInPtRangeAntiLambda;     //! Rapidity distribution
  TH1F        *fHistMCRapXi;                      //! Rapidity distribution
  TH1F        *fHistMCRapInPtRangeXi;             //! Rapidity distribution
  TH1F        *fHistMCRapPhi;                     //! Rapidity distribution
  TH1F        *fHistMCRapInPtRangePhi;            //! Rapidity distribution

  TH1F        *fHistMCPtK0s;                      //! Pt distribution K0s
  TH1F        *fHistMCPtLambda;                   //! Pt distribution Lambda

  TH1F        *fHistMCPtLambdaFromSigma;          //! Pt distribution of Lambda <- Sigma decay
  TH1F        *fHistMCPtAntiLambdaFromSigma;      //! Pt distribution of anti-Lambda <- Sigma decay

  TH1F        *fHistNTimesRecK0s;                 //! Multiple reconstruction studies
  TH1F        *fHistNTimesRecLambda;              //! Multiple reconstruction studies
  TH1F        *fHistNTimesRecAntiLambda;          //! Multiple reconstruction studies

  TH2F        *fHistNTimesRecK0sVsPt;             //! Multiple reconstruction studies
  TH2F        *fHistNTimesRecLambdaVsPt;          //! Multiple reconstruction studies
  TH2F        *fHistNTimesRecAntiLambdaVsPt;      //! Multiple reconstruction studies
  // ------------------------------------------------------

  // Reconstructed particle histograms  -------------------
  TH1F        *fHistNumberEvents;                 //! Number of events
  TH1F        *fHistTrackPerEvent;                //! Multiplicity
  TH1F        *fHistTrackletPerEvent;             //! Multiplicity
  TH1F        *fHistMCDaughterTrack;              //! Multiplicity

  TH1F        *fHistSPDPrimaryVertexZ;            //! Primary vertex

  TH1F        *fHistPrimaryVertexX;               //! Primary vertex
  TH1F        *fHistPrimaryVertexY;               //! Primary vertex
  TH1F        *fHistPrimaryVertexZ;               //! Primary vertex

  TH1F        *fHistPrimaryVertexResX;            //! Primary vertex resolution
  TH1F        *fHistPrimaryVertexResY;            //! Primary vertex resolution
  TH1F        *fHistPrimaryVertexResZ;            //! Primary vertex resolution

  TH1F        *fHistPrimaryVertexPosXV0events;    //! Primary vertex position in X in events with V0 candidates
  TH1F        *fHistPrimaryVertexPosYV0events;    //! Primary vertex position in Y in events with V0 candidates
  TH1F        *fHistPrimaryVertexPosZV0events;    //! Primary vertex position in Z in events with V0 candidates

  TH2F        *fHistDaughterPt;                   //! Daughters Pt

  TH2F        *fHistDcaPosToPrimVertex;           //! Cut checks
  TH2F        *fHistDcaNegToPrimVertex;           //! Cut checks
  TH2F        *fHistDcaPosToPrimVertexZoom;       //! Cut checks
  TH2F        *fHistDcaNegToPrimVertexZoom;       //! Cut checks
  TH2F        *fHistRadiusV0;                     //! Cut checks
  TH2F        *fHistDecayLengthV0;                //! Cut checks
  TH2F        *fHistDcaV0Daughters;               //! Cut checks
  TH2F        *fHistChi2;                         //! Cut checks
  TH2F        *fHistCosPointAngle;                //! Cut checks
  TH2F        *fHistCosPointAngleZoom;            //! Cut checks
  TH2F        *fHistProdRadius;                   //! Cut checks

  TH1F        *fHistV0Multiplicity;               //! V0 Multiplicity

  TH2F        *fHistChi2KFBeforeCutK0s;           //! Kalman filter Chi2
  TH2F        *fHistChi2KFBeforeCutLambda;        //! Kalman filter Chi2
  TH2F        *fHistChi2KFBeforeCutAntiLambda;    //! Kalman filter Chi2
  TH2F        *fHistChi2KFAfterCutK0s;            //! Kalman filter Chi2
  TH2F        *fHistChi2KFAfterCutLambda;         //! Kalman filter Chi2
  TH2F        *fHistChi2KFAfterCutAntiLambda;     //! Kalman filter Chi2

  TH1F        *fHistMassK0;                       //! Invariant mass
  TH1F        *fHistMassLambda;                   //! Invariant mass
  TH1F        *fHistMassAntiLambda;               //! Invariant mass

  TH2F        *fHistMassVsRadiusK0;               //! Invariant mass vs. radius
  TH2F        *fHistMassVsRadiusLambda;           //! Invariant mass vs. radius
  TH2F        *fHistMassVsRadiusAntiLambda;       //! Invariant mass vs. radius

  TH2F        *fHistPtVsMassK0;                   //! Pt vs. mass
  TH2F        *fHistPtVsMassLambda;               //! Pt vs. mass
  TH2F        *fHistArmenterosPodolanski;         //! Armenteros-Podolanski
  // ------------------------------------------------------


  // PID histograms  --------------------------------------
  TH1F        *fHistNsigmaPosPionAntiLambda;      //! Sigma positive pion
  TH1F        *fHistNsigmaNegProtonAntiLambda;    //! Sigma anti-proton
  TH1F        *fHistNsigmaPosProtonLambda;        //! Sigma proton
  TH1F        *fHistNsigmaNegPionLambda;          //! Sigma negative pion
  TH1F        *fHistNsigmaPosPionK0;              //! Sigma positive pion K0s
  TH1F        *fHistNsigmaNegPionK0;              //! Sigma negative pion K0s
  // ------------------------------------------------------

  // Associated particles ---------------------------------
  TH1F        *fHistAsMcRapK0;                    //! Rapidity distribution
  TH1F        *fHistAsMcRapLambda;                //! Rapidity distribution
  TH1F        *fHistAsMcRapAntiLambda;            //! Rapidity distribution

  TH1F        *fHistAsMcPtK0;                     //! Pt distribution
  TH1F        *fHistAsMcPtLambda;                 //! Pt distribution

  TH1F        *fHistAsMcPtZoomK0;                 //! Pt distribution
  TH1F        *fHistAsMcPtZoomLambda;             //! Pt distribution

  TH1F        *fHistAsMcProdRadiusK0;             //! Radius distribution
  TH1F        *fHistAsMcProdRadiusLambda;         //! Radius distribution
  TH1F        *fHistAsMcProdRadiusAntiLambda;     //! Radius distribution

  TH2F        *fHistAsMcProdRadiusXvsYK0s;        //! Radius distribution vs. rapidity
  TH2F        *fHistAsMcProdRadiusXvsYLambda;     //! Radius distribution vs. rapidity
  TH2F        *fHistAsMcProdRadiusXvsYAntiLambda; //! Radius distribution vs. rapidity

  TH1F        *fHistPidMcMassK0;                  //! Invariant mass distribution with PID checked
  TH1F        *fHistPidMcMassLambda;              //! Invariant mass distribution with PID checked
  TH1F        *fHistPidMcMassAntiLambda;          //! Invariant mass distribution with PID checked

  TH1F        *fHistAsMcMassK0;                   //! Invariant mass distribution
  TH1F        *fHistAsMcMassLambda;               //! Invariant mass distribution
  TH1F        *fHistAsMcMassAntiLambda;           //! Invariant mass distribution

  TH2F        *fHistAsMcPtVsMassK0;               //! Pt vs. invariant mass
  TH2F        *fHistAsMcPtVsMassLambda;           //! Pt vs. invariant mass
  TH2F        *fHistAsMcPtVsMassAntiLambda;       //! Pt vs. invariant mass

  TH2F        *fHistAsMcMassVsRadiusK0;           //! Invariant mass vs. radius
  TH2F        *fHistAsMcMassVsRadiusLambda;       //! Invariant mass vs. radius
  TH2F        *fHistAsMcMassVsRadiusAntiLambda;   //! Invariant mass vs. radius

  TH1F        *fHistAsMcResxK0;                   //! Position resolution for K0s
  TH1F        *fHistAsMcResyK0;                   //! Position resolution for K0s
  TH1F        *fHistAsMcReszK0;                   //! Position resolution for K0s

  TH2F        *fHistAsMcResrVsRadiusK0;           //! Position resolution vs. radius for K0s
  TH2F        *fHistAsMcReszVsRadiusK0;           //! Position resolution vs. radius for K0s

  TH1F        *fHistAsMcResxLambda;               //! Position resolution for Lambda
  TH1F        *fHistAsMcResyLambda;               //! Position resolution for Lambda
  TH1F        *fHistAsMcReszLambda;               //! Position resolution for Lambda

  TH2F        *fHistAsMcResrVsRadiusLambda;       //! Position resolution vs. radius for Lambda
  TH2F        *fHistAsMcReszVsRadiusLambda;       //! Position resolution vs. radius for Lambda

  TH1F        *fHistAsMcResxAntiLambda;           //! Position resolution for anti-Lambda
  TH1F        *fHistAsMcResyAntiLambda;           //! Position resolution for anti-Lambda
  TH1F        *fHistAsMcReszAntiLambda;           //! Position resolution for anti-Lambda

  TH2F        *fHistAsMcResrVsRadiusAntiLambda;   //! Position resolution vs. radius for anti-Lambda
  TH2F        *fHistAsMcReszVsRadiusAntiLambda;   //! Position resolution vs. radius for anti-Lambda

  TH1F        *fHistAsMcResPtK0;                  //! Pt resolution for K0s
  TH1F        *fHistAsMcResPtLambda;              //! Pt resolution for Lambda
  TH1F        *fHistAsMcResPtAntiLambda;          //! Pt resolution for anti-Lambda

  TH2F        *fHistAsMcResPtVsRapK0;             //! Pt resolution vs. rapidity for K0s
  TH2F        *fHistAsMcResPtVsRapLambda;         //! Pt resolution vs. rapidity for Lambda
  TH2F        *fHistAsMcResPtVsRapAntiLambda;     //! Pt resolution vs. rapidity for anti-Lambda

  TH2F        *fHistAsMcResPtVsPtK0;              //! Pt resolution vs. Pt for K0s	
  TH2F        *fHistAsMcResPtVsPtLambda;          //! Pt resolution vs. Pt for Lambda	
  TH2F        *fHistAsMcResPtVsPtAntiLambda;      //! Pt resolution vs. Pt for anti-Lambda

  TH1F        *fHistAsMcMotherPdgCodeK0s;         //! Pdg code of mother particle for K0s
  TH1F        *fHistAsMcMotherPdgCodeLambda;      //! Pdg code of mother particle for Lambda
  TH1F        *fHistAsMcMotherPdgCodeAntiLambda;  //! Pdg code of mother particle for anti-Lambda

  TH1F        *fHistAsMcPtLambdaFromSigma;        //! Pt distribution of Lambda <- Sigma decay
  TH1F        *fHistAsMcPtAntiLambdaFromSigma;    //! Pt distribution of anti-Lambda <- Sigma decay
  // ------------------------------------------------------

  // Associated secondary particle histograms -------------
  TH2F        *fHistAsMcSecondaryPtVsRapK0s;       //! Pt vs. rapidity distribution for K0s
  TH2F        *fHistAsMcSecondaryPtVsRapLambda;    //! Pt vs. rapidity distribution for Lambda
  TH2F        *fHistAsMcSecondaryPtVsRapAntiLambda;//! Pt vs. rapidity distribution for anti-Lambda

  TH1F        *fHistAsMcSecondaryProdRadiusK0s;       //! Production radius for K0s
  TH1F        *fHistAsMcSecondaryProdRadiusLambda;    //! Production radius for Lambda
  TH1F        *fHistAsMcSecondaryProdRadiusAntiLambda;//! Production radius for anti-Lambda

  TH2F        *fHistAsMcSecondaryProdRadiusXvsYK0s;       //! Production radius vs. rapidity for K0s
  TH2F        *fHistAsMcSecondaryProdRadiusXvsYLambda;    //! Production radius vs. rapidity for Lambda
  TH2F        *fHistAsMcSecondaryProdRadiusXvsYAntiLambda;//! Production radius vs. rapidity for anti-Lambda

  TH1F        *fHistAsMcSecondaryMotherPdgCodeK0s; //! Pdg code of mother particle for secondary K0s
  TH1F        *fHistAsMcSecondaryMotherPdgCodeLambda;//! Pdg code of mother particle for secondary for Lambda
  TH1F        *fHistAsMcSecondaryMotherPdgCodeAntiLambda;//! Pdg code of mother particle for secondary for anti-Lambda

  TH1F        *fHistAsMcSecondaryPtLambdaFromSigma;//! Pt distribution of secondary Lambda <- Sigma decay
  TH1F        *fHistAsMcSecondaryPtAntiLambdaFromSigma;//! Pt distribution of secondary anti-Lambda <- Sigma decay

  AliAnalysisTaskPerformanceStrange(const AliAnalysisTaskPerformanceStrange&); 
  AliAnalysisTaskPerformanceStrange& operator=(const AliAnalysisTaskPerformanceStrange&); 

  ClassDef(AliAnalysisTaskPerformanceStrange, 1); 
};

#endif
