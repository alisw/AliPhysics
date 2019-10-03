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
  TH1F        *fHistPtTracks;      //! Histo

  TH1F        *fHistMCMultiplicityPrimary;       //! Histo
  TH1F        *fHistMCMultiplicityTracks;       //! Histo
  TH1F        *fHistTPCTracks;                  //! Histo
  
  TH1F        *fHistMCPtAllK0s;       //! Histo
  TH1F        *fHistMCPtAllLambda;       //! Histo
  TH1F        *fHistMCPtAllAntiLambda;       //! Histo
  TH1F        *fHistMCPtAllXi;       //! Histo
  TH1F        *fHistMCPtAllAntiXi;       //! Histo
  TH1F        *fHistMCPtAllOmega;       //! Histo
  TH1F        *fHistMCPtAllAntiOmega;       //! Histo

  TH1F        *fHistMCRapK0s;                 //! Histo
  TH1F        *fHistMCRapLambda;              //! Histo
  TH1F        *fHistMCRapAntiLambda;          //! Histo
  TH1F        *fHistMCRapXi;                  //! Histo
//////////////////////////////////////////////////////////

  TH1F        *fHistMCPtK0s;       //! Histo
  TH1F        *fHistMCPtLambda;       //! Histo
  TH1F        *fHistMCPtAntiLambda;       //! Histo
//////////////////////////////////////////////////////////

  // ESD histograms
  TH1F        *fHistNumberEvents;        //! Histo
  TH1F        *fHistTrackPerEvent;       //! Histo
  TH1F        *fHistTPCMult;             //! Histo
  TH1F        *fHistTrackletPerEvent;   //! Histo
  TH1F        *fHistSPDPrimaryVertexZ;       //! Histo
  TH1F        *fHistPrimaryVertexX;       //! Histo
  TH1F        *fHistPrimaryVertexY;       //! Histo
  TH1F        *fHistPrimaryVertexZ;       //! Histo
//////////////////////////////////////////////////////////////////////

  TH1F        *fHistV0Multiplicity;  //! Histo
  TH1F        *fHistMassK0;       //! Histo
  TH1F        *fHistMassLambda;       //! Histo
  TH1F        *fHistMassAntiLambda;       //! Histo
  TH1F        *fHistMassXi;       //! Histo
  TH1F        *fHistMassAntiXi;       //! Histo
  TH1F        *fHistMassOmega;       //! Histo
  TH1F        *fHistMassAntiOmega;       //! Histo

  TH2F        *fHistMassXiVsPID;       //! Histo
////////////////////////////////////////////////////////////////////////////

  TH2F        *fHistPtVsMassK0;       //! Histo
  TH2F        *fHistPtVsMassLambda;       //! Histo
  TH2F        *fHistPtVsMassAntiLambda;       //! Histo
/////////////////////////////////////////////

  TH2F        *fHistArmenterosPodolanski;       //! Histo
  TH2F        *fHistK0sMassVsLambdaMass;       //! Histo

  //PID check
  TH2F *fHistTPCsigPLambda;               //! Histo
  TH2F *fHistTPCsigPAntiLambda;               //! Histo
  TH1F *fHistNSigmaProton;               //! Histo

  // Associated particles histograms
  TH1F        *fHistAsMcRapK0;       //! Histo
  TH1F        *fHistAsMcRapLambda;       //! Histo
  TH1F        *fHistAsMcRapAntiLambda;       //! Histo

////////////////////////////////////////////////////////////////////
  TH1F        *fHistAsMcPtK0;       //! Histo
  TH1F        *fHistAsMcPtLambda;       //! Histo
  TH1F        *fHistAsMcPtAntiLambda;       //! Histo

  TH1F        *fHistPidMcMassK0;       //! Histo
  TH1F        *fHistPidMcMassLambda;       //! Histo
  TH1F        *fHistPidMcMassAntiLambda;       //! Histo

    //Mass 

  TH1F        *fHistAsMcMassK0;       //! Histo
  TH1F        *fHistAsMcMassLambda;       //! Histo
  TH1F        *fHistAsMcMassAntiLambda;       //! Histo

  //PtVsMass

  TH2F        *fHistAsMcPtVsMassK0;       //! Histo
  TH2F        *fHistAsMcPtVsMassLambda;       //! Histo
  TH2F        *fHistAsMcPtVsMassAntiLambda;       //! Histo

  TH1F        *fHistCompositionXi;             //! Histo
  TH1F        *fHistCompositionAntiXi;          //! Histo
  TH1F        *fHistCompositionOmega;           //! Histo
  TH1F        *fHistCompositionAntiOmega;       //! Histo
  TH1I        *fHistMCIndexes;          //! Histo


  AliAnalysisTaskPerformanceStrange(const AliAnalysisTaskPerformanceStrange&); 
  AliAnalysisTaskPerformanceStrange& operator=(const AliAnalysisTaskPerformanceStrange&); 

  ClassDef(AliAnalysisTaskPerformanceStrange, 1); 
};

#endif
