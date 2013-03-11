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
class TH3F;
class AliAnalysisCentralitySelector;
class AliPIDResponse;
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPerformanceStrange : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPerformanceStrange();
  AliAnalysisTaskPerformanceStrange(const char *name);
  virtual ~AliAnalysisTaskPerformanceStrange(); // Destructor implemented by Kalinak  
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
  void   SetCollidingSystems(Int_t collidingSystems = 0) {fCollidingSystems = collidingSystems;}
  void   SetAnalysisMC(Bool_t analysisMC) {fAnalysisMC = analysisMC;}
  void   SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  void   SetUsePID(const char* usePID) {fUsePID = usePID;}
  void   SetAnalysisCut(const char* useCut) {fUseCut = useCut;}
  void   SetCentralityRange(Int_t down, Int_t up) {fDown=down; fUp = up;}
  void   SetTrackCuts(AliESDtrackCuts * myTracksCuts) { fTracksCuts = myTracksCuts;}
  void   SetCentralitySelector(AliAnalysisCentralitySelector * centr) { fCentrSelector = centr;}
  void   SetQASelector(Bool_t QA = 0) { fQASelector = QA;}
  Double_t MyRapidity(Double_t rE, Double_t rPz) const;
 
 private:
  Double_t fCuts[7];                            //! V0 finding cuts
  Bool_t       fAnalysisMC;                     //  0->No MC or 1->MC analysis
  TString      fAnalysisType;                   //  "ESD" or "AOD"
  Bool_t       fCollidingSystems;               //  Colliding systems 0/1 for pp/PbPb  
  TString      fUsePID;                         //  "withPID" or "noPID"
  TString      fUseCut;                         //  "yes" or "no"
  Int_t		fDown;				//centrality range 
  Int_t		fUp;                            //centrality range
  AliESDEvent *fESD;                            //! ESD object
  TList       *fListHist;		//! Output List


  AliAnalysisCentralitySelector * fCentrSelector; // Centrality selector, used to 
  AliESDtrackCuts * fTracksCuts;		// track cuts
  AliPIDResponse *fPIDResponse;                 // PID response
  Bool_t      fQASelector;                    // Quality Assurenc Histo switch

  // MC histograms
  TH1F        *fHistMCPrimaryVertexX;      //! Histo
  TH1F        *fHistMCPrimaryVertexY;      //! Histo
  TH1F        *fHistMCPrimaryVertexZ;      //! Histo
  TH1F        *fHistPtTracksITSRefit;      //! Histo
  TH1F        *fHistPtTracks;      //! Histo
  TH1F        *fHistPtTracksPITSRefit;      //! Histo
  TH1F        *fHistPtTracksP;      //! Histo

  TH1F        *fHistMCMultiplicityPrimary;       //! Histo
  TH1F        *fHistMCMultiplicityTracks;       //! Histo
  TH1F        *fHistTPCTracks;                  //! Histo
  
  TH2F        *fHistMCtracksProdRadiusK0s;       //! Histo
  TH2F        *fHistMCtracksProdRadiusLambda;       //! Histo
  TH2F        *fHistMCtracksProdRadiusAntiLambda;       //! Histo

  TH1F        *fHistMCtracksDecayRadiusK0s;       //! Histo
  TH1F        *fHistMCtracksDecayRadiusLambda;       //! Histo
  TH1F        *fHistMCtracksDecayRadiusAntiLambda;       //! Histo

  TH1F        *fHistMCPtAllK0s;       //! Histo
  TH1F        *fHistMCPtAllLambda;       //! Histo
  TH1F        *fHistMCPtAllAntiLambda;       //! Histo

  //Rap3
  TH1F        *fHistMCPtAllK0sRap3;       //! Histo
  TH1F        *fHistMCPtAllLambdaRap3;       //! Histo
  TH1F        *fHistMCPtAllAntiLambdaRap3;       //! Histo

  TH1F        *fHistMCProdRadiusK0s;       //! Histo
  TH1F        *fHistMCProdRadiusLambda;       //! Histo
  TH1F        *fHistMCProdRadiusAntiLambda;       //! Histo

  TH1F    *fHistMCPrimDecayRadiusK0s;          //! Histo
  TH1F    *fHistMCPrimDecayRadiusLambda;          //! Histo
  TH1F    *fHistMCPrimDecayRadiusAntiLambda;     //! Histo


  TH1F        *fHistMCRapK0s;                 //! Histo
  TH1F        *fHistMCRapInPtRangeK0s;        //! Histo
  TH1F        *fHistMCRapLambda;              //! Histo
  TH1F        *fHistMCRapInPtRangeLambda;     //! Histo
  TH1F        *fHistMCRapAntiLambda;          //! Histo
  TH1F        *fHistMCRapInPtRangeAntiLambda; //! Histo
  TH1F        *fHistMCRapXi;                  //! Histo
  TH1F        *fHistMCRapInPtRangeXi;         //! Histo
  TH1F        *fHistMCRapPhi;                 //! Histo
  TH1F        *fHistMCRapInPtRangePhi;        //! Histo
//////////////////////////////////////////////////////////

  TH1F        *fHistMCPtK0s;       //! Histo
  TH1F        *fHistMCPtLambda;       //! Histo
  TH1F        *fHistMCPtAntiLambda;       //! Histo

  //Rap3
  TH1F        *fHistMCPtK0sRap3;       //! Histo
  TH1F        *fHistMCPtLambdaRap3;       //! Histo
  TH1F        *fHistMCPtAntiLambdaRap3;       //! Histo

//////////////////////////////////////////////////////////



  TH1F        *fHistMCPtLambdaFromSigma;       //! Histo
  TH1F        *fHistMCPtAntiLambdaFromSigma;       //! Histo

  TH1F        *fHistNTimesRecK0s;       //! Histo
  TH1F        *fHistNTimesRecLambda;       //! Histo
  TH1F        *fHistNTimesRecAntiLambda;       //! Histo

  TH2F        *fHistNTimesRecK0sVsPt;       //! Histo
  TH2F        *fHistNTimesRecLambdaVsPt;       //! Histo
  TH2F        *fHistNTimesRecAntiLambdaVsPt;       //! Histo


  // ESD histograms
  TH1F        *fHistNumberEvents;        //! Histo
  TH1F        *fHistTrackPerEvent;       //! Histo

  TH1F        *fHistTPCMult;             //! Histo

  TH1F        *fHistTrackletPerEvent;   //! Histo
  TH1F        *fHistMCDaughterTrack;       //! Histo

  TH1F        *fHistSPDPrimaryVertexZ;       //! Histo

  TH1F        *fHistPrimaryVertexX;       //! Histo
  TH1F        *fHistPrimaryVertexY;       //! Histo
  TH1F        *fHistPrimaryVertexZ;       //! Histo

  TH1F        *fHistPrimaryVertexResX;       //! Histo
  TH1F        *fHistPrimaryVertexResY;       //! Histo
  TH1F        *fHistPrimaryVertexResZ;       //! Histo

  TH1F        *fHistPrimaryVertexPosXV0events;  //! Primary vertex position in X in events with V0 candidates
  TH1F        *fHistPrimaryVertexPosYV0events;  //! Primary vertex position in Y in events with V0 candidates
  TH1F        *fHistPrimaryVertexPosZV0events;  //! Primary vertex position in Z in events with V0 candidates

  TH2F        *fHistDaughterPt;               //! Histo

///////////////////////////K0s 2D histos: cut vs on fly status/////////////////

  TH2F        *fHistDcaPosToPrimVertexK0;       //! Histo
  TH2F        *fHistDcaNegToPrimVertexK0;       //! Histo
//  TH2F        *fHistDcaPosToPrimVertexZoomK0;       //! Histo
//  TH2F        *fHistDcaNegToPrimVertexZoomK0;       //! Histo
  TH2F        *fHistRadiusV0K0;       //! Histo
  TH2F        *fHistDecayLengthV0K0;       //! Histo
  TH2F        *fHistDcaV0DaughtersK0;       //! Histo
  TH2F        *fHistChi2K0;       //! Histo
  TH2F        *fHistCosPointAngleK0;       //! Histo
//  TH2F        *fHistCosPointAngleZoomK0;       //! Histo
//  TH2F        *fHistProdRadiusK0;       //! Histo


///////////////////////////K0s 2D histos: cut vs mass//////////////
  TH2F        *fHistDcaPosToPrimVertexK0vsMassK0;  //! Histo
  TH2F        *fHistDcaNegToPrimVertexK0vsMassK0;  //! Histo
  TH2F        *fHistRadiusV0K0vsMassK0;            //! Histo
  TH2F        *fHistDecayLengthV0K0vsMassK0;       //! Histo
  TH2F        *fHistDcaV0DaughtersK0vsMassK0;      //! Histo
  TH2F        *fHistCosPointAngleK0vsMassK0;       //! Histo
  
  // pt1
  TH2F        *fHistDcaPosToPrimVertexK0vsMassK0pt1;  //! Histo
  TH2F        *fHistDcaNegToPrimVertexK0vsMassK0pt1;  //! Histo
  TH2F        *fHistRadiusV0K0vsMassK0pt1;            //! Histo
  TH2F        *fHistDecayLengthV0K0vsMassK0pt1;       //! Histo
  TH2F        *fHistDcaV0DaughtersK0vsMassK0pt1;      //! Histo
  TH2F        *fHistCosPointAngleK0vsMassK0pt1;       //! Histo
  
  // pt2
  TH2F        *fHistDcaPosToPrimVertexK0vsMassK0pt2;  //! Histo
  TH2F        *fHistDcaNegToPrimVertexK0vsMassK0pt2;  //! Histo
  TH2F        *fHistRadiusV0K0vsMassK0pt2;             //! Histo
  TH2F        *fHistDecayLengthV0K0vsMassK0pt2;     //! Histo
  TH2F        *fHistDcaV0DaughtersK0vsMassK0pt2;    //! Histo
  TH2F        *fHistCosPointAngleK0vsMassK0pt2;     //! Histo

  // pt3
  TH2F        *fHistDcaPosToPrimVertexK0vsMassK0pt3;    //! Histo
  TH2F        *fHistDcaNegToPrimVertexK0vsMassK0pt3;    //! Histo
  TH2F        *fHistRadiusV0K0vsMassK0pt3;             //! Histo
  TH2F        *fHistDecayLengthV0K0vsMassK0pt3;       //! Histo
  TH2F        *fHistDcaV0DaughtersK0vsMassK0pt3;      //! Histo
  TH2F        *fHistCosPointAngleK0vsMassK0pt3;      //! Histo

//////////////////////////Lambda 2D histos: cut vs on fly status////////////////////

  TH2F        *fHistDcaPosToPrimVertexL;       //! Histo
  TH2F        *fHistDcaNegToPrimVertexL;       //! Histo
//  TH2F        *fHistDcaPosToPrimVertexZoomL;       //! Histo
//  TH2F        *fHistDcaNegToPrimVertexZoomL;       //! Histo
  TH2F        *fHistRadiusV0L;       //! Histo
  TH2F        *fHistDecayLengthV0L;       //! Histo
  TH2F        *fHistDcaV0DaughtersL;       //! Histo
  TH2F        *fHistChi2L;       //! Histo
  TH2F        *fHistCosPointAngleL;       //! Histo

//  TH2F        *fHistCosPointAngleZoomL;       //! Histo
//  TH2F        *fHistProdRadiusL;       //! Histo    

//////////////////////////Lambda 2D histos: cut vs mass////////////////
  TH2F        *fHistDcaPosToPrimVertexLvsMassL;      //! Histo
  TH2F        *fHistDcaNegToPrimVertexLvsMassL;      //! Histo
  TH2F        *fHistRadiusV0LvsMassL;                 //! Histo
  TH2F        *fHistDecayLengthV0LvsMassL;            //! Histo
  TH2F        *fHistDcaV0DaughtersLvsMassL;         //! Histo
  TH2F        *fHistCosPointAngleLvsMassL;            //! Histo
  TH3F        *fHistCosPointAngleLvsMassVsPtsigL;    //! Histo
  TH3F        *fHistCosPointAngleLvsMassVsPtbackL;    //! Histo



  // pt1
  TH2F        *fHistDcaPosToPrimVertexLambdaVsMasspt1;  //! Histo
  TH2F        *fHistDcaNegToPrimVertexLambdaVsMasspt1;  //! Histo
  TH2F        *fHistRadiusV0LambdaVsMasspt1;            //! Histo
  TH2F        *fHistDecayLengthV0LambdaVsMasspt1;       //! Histo
  TH2F        *fHistDcaV0DaughtersLambdaVsMasspt1;      //! Histo
  TH2F        *fHistCosPointAngleLambdaVsMasspt1;       //! Histo
  
  // pt2
  TH2F        *fHistDcaPosToPrimVertexLambdaVsMasspt2;  //! Histo
  TH2F        *fHistDcaNegToPrimVertexLambdaVsMasspt2;  //! Histo
  TH2F        *fHistRadiusV0LambdaVsMasspt2;             //! Histo
  TH2F        *fHistDecayLengthV0LambdaVsMasspt2;     //! Histo
  TH2F        *fHistDcaV0DaughtersLambdaVsMasspt2;    //! Histo
  TH2F        *fHistCosPointAngleLambdaVsMasspt2;     //! Histo

  // pt3
  TH2F        *fHistDcaPosToPrimVertexLambdaVsMasspt3;    //! Histo
  TH2F        *fHistDcaNegToPrimVertexLambdaVsMasspt3;    //! Histo
  TH2F        *fHistRadiusV0LambdaVsMasspt3;             //! Histo
  TH2F        *fHistDecayLengthV0LambdaVsMasspt3;       //! Histo
  TH2F        *fHistDcaV0DaughtersLambdaVsMasspt3;      //! Histo
  TH2F        *fHistCosPointAngleLambdaVsMasspt3;      //! Histo


//////////////////////////Lambda 2D histos: cut vs on fly status////////////////////

  TH2F        *fHistDcaPosToPrimVertexAntiL;       //! Histo
  TH2F        *fHistDcaNegToPrimVertexAntiL;       //! Histo
//  TH2F        *fHistDcaPosToPrimVertexZoomL;       //! Histo
//  TH2F        *fHistDcaNegToPrimVertexZoomL;       //! Histo
  TH2F        *fHistRadiusV0AntiL;       //! Histo
  TH2F        *fHistDecayLengthV0AntiL;       //! Histo
  TH2F        *fHistDcaV0DaughtersAntiL;       //! Histo
  TH2F        *fHistChi2AntiL;       //! Histo
  TH2F        *fHistCosPointAngleAntiL;       //! Histo
//  TH2F        *fHistCosPointAngleZoomL;       //! Histo
//  TH2F        *fHistProdRadiusL;       //! Histo    

//////////////////////////Lambda 2D histos: cut vs mass////////////////
  TH2F        *fHistDcaPosToPrimVertexAntiLvsMass;      //! Histo
  TH2F        *fHistDcaNegToPrimVertexAntiLvsMass;      //! Histo
  TH2F        *fHistRadiusV0AntiLvsMass;                 //! Histo
  TH2F        *fHistDecayLengthV0AntiLvsMass;            //! Histo
  TH2F        *fHistDcaV0DaughtersAntiLvsMass;         //! Histo
  TH2F        *fHistCosPointAngleAntiLvsMass;            //! Histo



  // pt1
  TH2F        *fHistDcaPosToPrimVertexAntiLVsMasspt1;  //! Histo
  TH2F        *fHistDcaNegToPrimVertexAntiLVsMasspt1;  //! Histo
  TH2F        *fHistRadiusV0AntiLVsMasspt1;            //! Histo
  TH2F        *fHistDecayLengthV0AntiLVsMasspt1;       //! Histo
  TH2F        *fHistDcaV0DaughtersAntiLVsMasspt1;      //! Histo
  TH2F        *fHistCosPointAngleAntiLVsMasspt1;       //! Histo
  
  // pt2
  TH2F        *fHistDcaPosToPrimVertexAntiLVsMasspt2;  //! Histo
  TH2F        *fHistDcaNegToPrimVertexAntiLVsMasspt2;  //! Histo
  TH2F        *fHistRadiusV0AntiLVsMasspt2;             //! Histo
  TH2F        *fHistDecayLengthV0AntiLVsMasspt2;     //! Histo
  TH2F        *fHistDcaV0DaughtersAntiLVsMasspt2;    //! Histo
  TH2F        *fHistCosPointAngleAntiLVsMasspt2;     //! Histo

  // pt3
  TH2F        *fHistDcaPosToPrimVertexAntiLVsMasspt3;    //! Histo
  TH2F        *fHistDcaNegToPrimVertexAntiLVsMasspt3;    //! Histo
  TH2F        *fHistRadiusV0AntiLVsMasspt3;             //! Histo
  TH2F        *fHistDecayLengthV0AntiLVsMasspt3;       //! Histo
  TH2F        *fHistDcaV0DaughtersAntiLVsMasspt3;      //! Histo
  TH2F        *fHistCosPointAngleAntiLVsMasspt3;      //! Histo



//////////////////////////////////////////////////////////////////////


  TH1F        *fHistV0Multiplicity;  //! Histo
  TH1F        *fHistMassK0;       //! Histo
  TH1F        *fHistMassLambda;       //! Histo
  TH1F        *fHistMassAntiLambda;       //! Histo
  TH2F        *fHistMassVsRadiusK0;       //! Histo
  TH2F        *fHistMassVsRadiusLambda;       //! Histo
  TH2F        *fHistMassVsRadiusAntiLambda;       //! Histo

////////////////////////////////////////////////////////////////////////////

  TH2F        *fHistPtVsMassK0;       //! Histo
  TH2F        *fHistPtVsMassLambda;       //! Histo
  TH2F        *fHistPtVsMassAntiLambda;       //! Histo

  //Rap3
  TH2F        *fHistPtVsMassK0Rap3;       //! Histo
  TH2F        *fHistPtVsMassLambdaRap3;       //! Histo
  TH2F        *fHistPtVsMassAntiLambdaRap3;       //! Histo

  TH2F        *fHistTranscTauVsMassL;                 //! Histo
  TH2F        *fHistTranscTauVsMassAntiL;                 //! Histo
  TH2F        *fHistTranscTauVsMassK0s;                 //! Histo

  // cTauVsMass Rap3
  TH2F        *fHistTranscTauVsMassLRap3;       //! Histo
  TH2F        *fHistTranscTauVsMassAntiLRap3;       //! Histo
  TH2F        *fHistTranscTauVsMassK0sRap3;       //! Histo

  // cTauVsMass Low pt
  TH2F    *fHistTranscTauVsMassLptLow;        //! Histo
  TH2F    *fHistTranscTauVsMassAntiLptLow;    //! Histo
  TH2F    *fHistTranscTauVsMassK0sptLow;       //! Histo

  //cTauVsMass Low pt Rap3 
  TH2F    *fHistTranscTauVsMassLptLowRap3;        //! Histo
  TH2F    *fHistTranscTauVsMassAntiLptLowRap3;    //! Histo
  TH2F    *fHistTranscTauVsMassK0sptLowRap3;       //! Histo

/////////////////////////////////////////////

  TH2F        *fHistArmenterosPodolanski;       //! Histo
  TH2F        *fHistK0sMassVsLambdaMass;       //! Histo

  //PID check
  TH2F *fHistTPCsigPLambda;               //! Histo
  TH2F *fHistTPCsigPAntiLambda;               //! Histo
  TH1F *fHistNSigmaProton;               //! Histo

  //PID
  TH1F        *fHistNsigmaPosPionAntiLambda;    //! Histo
  TH1F        *fHistNsigmaNegProtonAntiLambda;   //! Histo
  TH1F        *fHistNsigmaPosProtonLambda;        //! Histo
  TH1F        *fHistNsigmaNegPionLambda;           //! Histo
  TH1F        *fHistNsigmaPosProtonAntiLambda;        //! Histo
  TH1F        *fHistNsigmaNegPionAntiLambda;           //! Histo
  TH1F        *fHistNsigmaPosPionK0;                //! Histo
  TH1F        *fHistNsigmaNegPionK0;                 //! Histo

  // Associated particles histograms
  TH1F        *fHistAsMcRapK0;       //! Histo
  TH1F        *fHistAsMcRapLambda;       //! Histo
  TH1F        *fHistAsMcRapAntiLambda;       //! Histo

////////////////////////////////////////////////////////////////////
  TH1F        *fHistAsMcPtK0;       //! Histo
  TH1F        *fHistAsMcPtLambda;       //! Histo
  TH1F        *fHistAsMcPtAntiLambda;       //! Histo

  TH1F        *fHistAsMcTranscTauL;       //! Histo
  TH1F        *fHistAsMcTranscTauAntiL;       //! Histo
  TH1F        *fHistAsMcTranscTauK0s;       //! Histo

  TH1F        *fHistAsMcTranscTauLRap3;       //! Histo
  TH1F        *fHistAsMcTranscTauAntiLRap3;       //! Histo
  TH1F        *fHistAsMcTranscTauK0sRap3;       //! Histo


  TH1F    *fHistAsMcTranscTauLptLow;        //! Histo
  TH1F    *fHistAsMcTranscTauAntiLptLow;    //! Histo
  TH1F    *fHistAsMcTranscTauK0sptLow;       //! Histo

  //Rap3
  TH1F    *fHistAsMcTranscTauLptLowRap3;        //! Histo
  TH1F    *fHistAsMcTranscTauAntiLptLowRap3;    //! Histo
  TH1F    *fHistAsMcTranscTauK0sptLowRap3;       //! Histo


  //Rap3
  TH1F        *fHistAsMcPtK0Rap3;       //! Histo
  TH1F        *fHistAsMcPtLambdaRap3;       //! Histo
  TH1F        *fHistAsMcPtAntiLambdaRap3;       //! Histo


/////////////////////////////////////////////////////////////////////


  TH1F        *fHistAsMcPtZoomK0;       //! Histo
  TH1F        *fHistAsMcPtZoomLambda;       //! Histo
  TH1F        *fHistAsMcPtZoomAntiLambda;       //! Histo

  TH1F        *fHistAsMcProdRadiusK0;       //! Histo
  TH1F        *fHistAsMcProdRadiusLambda;       //! Histo
  TH1F        *fHistAsMcProdRadiusAntiLambda;       //! Histo

  TH2F        *fHistAsMcProdRadiusXvsYK0s;       //! Histo
  TH2F        *fHistAsMcProdRadiusXvsYLambda;       //! Histo
  TH2F        *fHistAsMcProdRadiusXvsYAntiLambda;       //! Histo

  TH1F        *fHistPidMcMassK0;       //! Histo
  TH1F        *fHistPidMcMassLambda;       //! Histo
  TH1F        *fHistPidMcMassAntiLambda;       //! Histo

    //Mass 

  TH1F        *fHistAsMcMassK0;       //! Histo
  TH1F        *fHistAsMcMassLambda;       //! Histo
  TH1F        *fHistAsMcMassAntiLambda;       //! Histo

  //Rap3
  TH1F        *fHistAsMcMassK0Rap3;       //! Histo
  TH1F        *fHistAsMcMassLambdaRap3;       //! Histo
  TH1F        *fHistAsMcMassAntiLambdaRap3;       //! Histo


  //PtVsMass

  TH2F        *fHistAsMcPtVsMassK0;       //! Histo
  TH2F        *fHistAsMcPtVsMassLambda;       //! Histo
  TH2F        *fHistAsMcPtVsMassAntiLambda;       //! Histo

  //Rap3
  TH2F        *fHistAsMcPtVsMassK0Rap3;       //! Histo
  TH2F        *fHistAsMcPtVsMassLambdaRap3;       //! Histo
  TH2F        *fHistAsMcPtVsMassAntiLambdaRap3;       //! Histo

  TH2F        *fHistAsMcMassVsRadiusK0;       //! Histo
  TH2F        *fHistAsMcMassVsRadiusLambda;       //! Histo
  TH2F        *fHistAsMcMassVsRadiusAntiLambda;       //! Histo

  TH1F        *fHistAsMcResxK0;       //! Histo
  TH1F        *fHistAsMcResyK0;       //! Histo
  TH1F        *fHistAsMcReszK0;       //! Histo

  TH2F        *fHistAsMcResrVsRadiusK0;       //! Histo
  TH2F        *fHistAsMcReszVsRadiusK0;       //! Histo


  TH1F        *fHistAsMcResxLambda;       //! Histo
  TH1F        *fHistAsMcResyLambda;       //! Histo
  TH1F        *fHistAsMcReszLambda;       //! Histo

  TH2F        *fHistAsMcResrVsRadiusLambda;       //! Histo
  TH2F        *fHistAsMcReszVsRadiusLambda;       //! Histo
    

  TH1F        *fHistAsMcResxAntiLambda;       //! Histo
  TH1F        *fHistAsMcResyAntiLambda;       //! Histo
  TH1F        *fHistAsMcReszAntiLambda;       //! Histo

  TH2F        *fHistAsMcResrVsRadiusAntiLambda;       //! Histo
  TH2F        *fHistAsMcReszVsRadiusAntiLambda;       //! Histo
    

  TH1F        *fHistAsMcResPtK0;       //! Histo
  TH1F        *fHistAsMcResPtLambda;       //! Histo
  TH1F        *fHistAsMcResPtAntiLambda;       //! Histo


  TH2F        *fHistAsMcResPtVsRapK0;       //! Histo
  TH2F        *fHistAsMcResPtVsRapLambda;       //! Histo
  TH2F        *fHistAsMcResPtVsRapAntiLambda;       //! Histo
  TH2F        *fHistAsMcResPtVsPtK0;       //! Histo
  TH2F        *fHistAsMcResPtVsPtLambda;       //! Histo
  TH2F        *fHistAsMcResPtVsPtAntiLambda;       //! Histo
  

  TH1F        *fHistAsMcMotherPdgCodeK0s;       //! Histo
  TH1F        *fHistAsMcMotherPdgCodeLambda;       //! Histo
  TH1F        *fHistAsMcMotherPdgCodeAntiLambda;       //! Histo
  TH1F        *fHistAsMcPtLambdaFromSigma;       //! Histo
  TH1F        *fHistAsMcPtAntiLambdaFromSigma;       //! Histo


  // Associated secondary particles:
  TH2F        *fHistAsMcSecondaryPtVsRapK0s;       //! Histo
  TH2F        *fHistAsMcSecondaryPtVsRapLambda;       //! Histo
  TH2F        *fHistAsMcSecondaryPtVsRapAntiLambda;       //! Histo

  TH1F        *fHistAsMcSecondaryProdRadiusK0s;       //! Histo
  TH1F        *fHistAsMcSecondaryProdRadiusLambda;       //! Histo
  TH1F        *fHistAsMcSecondaryProdRadiusAntiLambda;       //! Histo

  TH2F        *fHistAsMcSecondaryProdRadiusXvsYK0s;       //! Histo
  TH2F        *fHistAsMcSecondaryProdRadiusXvsYLambda;       //! Histo
  TH2F        *fHistAsMcSecondaryProdRadiusXvsYAntiLambda;       //! Histo

  TH1F        *fHistAsMcSecondaryMotherPdgCodeK0s;       //! Histo
  TH1F        *fHistAsMcSecondaryMotherPdgCodeLambda;       //! Histo
  TH1F        *fHistAsMcSecondaryMotherPdgCodeAntiLambda;       //! Histo

  TH1F        *fHistAsMcSecondaryPtLambdaFromSigma;       //! Histo
  TH1F        *fHistAsMcSecondaryPtAntiLambdaFromSigma;       //! Histo


  AliAnalysisTaskPerformanceStrange(const AliAnalysisTaskPerformanceStrange&); 
  AliAnalysisTaskPerformanceStrange& operator=(const AliAnalysisTaskPerformanceStrange&); 

  ClassDef(AliAnalysisTaskPerformanceStrange, 1); 
};

#endif
