#ifndef ALIANATRANSVERSEEVENTSHAPE_H
#define ALIANATRANSVERSEEVENTSHAPE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */


// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TProfile.h>
#include <TTreeStream.h>
#include <TRandom.h>
#include <TObject.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliMCEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliAODMCParticle.h> 
#include <AliESDtrackCuts.h>
#include <AliTransverseEventShape.h>

class AliPPVsMultUtils;

class AliAnaTransverseEventShapeTask : public AliAnalysisTaskSE {
 public:
  
  
  AliAnaTransverseEventShapeTask();
  AliAnaTransverseEventShapeTask(const char *name);
  virtual ~AliAnaTransverseEventShapeTask();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t*);
  virtual Float_t GetTest();

  Bool_t   GetAnalysisMC() { return fAnalysisMC; }   
  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }     

  //for event shape analysis
  virtual void  SetUseHybridESA(Bool_t usehyb) {fUseHybrid = usehyb;}
  virtual void  SetTrackFilterESAHyb1(AliAnalysisFilter* trackH1F) {fTrackFilterHybrid1 = trackH1F;}
  virtual void  SetTrackFilterESAHyb2(AliAnalysisFilter* trackH2F) {fTrackFilterHybrid2 = trackH2F;}
  virtual void  SetTrackFilterESA(AliAnalysisFilter* trackF) {fTrackFilterESA = trackF;}
  virtual void  SetMinMultForESA(Int_t minnch)     {fMinMultESA = minnch;}
  virtual void  SetStepSizeESA(Float_t sizestep)   {fSizeStepESA = sizestep;}
  virtual void  SetIsEtaAbsESA(Bool_t isabseta){fIsAbsEtaESA = isabseta;}
  virtual void  SetTrackEtaMinESA(Float_t etaminF) {fEtaMinCutESA = etaminF;}
  virtual void  SetTrackEtaMaxESA(Float_t etamaxF) {fEtaMaxCutESA = etamaxF;}
  virtual void  SetTrackPtMinESA(Float_t ptminF) {fPtMinCutESA = ptminF;}
  virtual void  SetTrackPtMaxESA(Float_t ptmaxF) {fPtMaxCutESA = ptmaxF;}



  virtual void  SetTrigger(UInt_t ktriggerInt) {ftrigBit = ktriggerInt;}
  virtual void  SetCentralityEstimator(const char * centEst) {fCentEst = centEst;}
  virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;}   
  virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
  virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
  virtual void  SetStoreMcIn(Bool_t value) {fStoreMcIn = value;}
  virtual void  SetAnalysisPbPb(Bool_t isanaPbPb) {fAnalysisPbPb = isanaPbPb;}
  
 private:
  virtual Float_t GetVertex(const AliVEvent* event) const;
  
  AliESDEvent*      fESD;                  //! ESD object
  AliAODEvent*      fAOD;                  //! AOD object
  AliPPVsMultUtils *fPPVsMultUtils;   //! object for Vzero based multiplicity
  AliMCEvent*       fMC;                   //! MC object
  AliStack*         fMCStack;              //! MC ESD stack
  TClonesArray*     fMCArray;             //! MC array for AOD

  Bool_t            fUseHybrid;
  AliAnalysisFilter *fTrackFilterHybrid1;
  AliAnalysisFilter *fTrackFilterHybrid2;
  AliAnalysisFilter *fTrackFilterESA; // Track filter for Event Shapes
  Int_t             fMinMultESA;
  Float_t           fSizeStepESA;
  Bool_t            fIsAbsEtaESA;
  Float_t           fEtaMaxCutESA;
  Float_t           fEtaMinCutESA;
  Float_t           fPtMaxCutESA;
  Float_t           fPtMinCutESA;
  
  TString           fCentEst;             // V0A , V0M, 
  TString           fAnalysisType;        //  "ESD" or "AOD"
  Bool_t            fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
  Bool_t            fAnalysisPbPb;        //  true you want to analyze PbPb data, false for pp
  UInt_t            ftrigBit;
  TRandom*          fRandom;              //! random number generator
  Bool_t            fPileUpRej;           // kTRUE is pile-up is rejected


  //
  // Cuts and options
  //

  Double_t     fVtxCut;             // Vtx cut on z position in cm
  Double_t     fEtaCut;             // Eta cut used to select particles
  Float_t      fMinCent; //minimum centrality
  Float_t      fMaxCent; //maximum centrality
  Bool_t       fStoreMcIn;          // Store MC input tracks
  //
  // Help variables
  //
  Short_t      fMcProcessType;      // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
  Short_t      fTriggeredEventMB;   // 1 = triggered, 0 = not trigged (MC only)
  Short_t      fVtxStatus;          // -1 = no vtx, 0 = outside cut, 1 = inside cut
  Float_t      fZvtx;               // z vertex
  Float_t      fZvtxMC;             // z vertex MC (truth)
  Int_t        fRun;                // run no
  ULong64_t    fEventId;            // unique event id

  //
  // Output objects
  //
  TList*        fListOfObjects;     //! Output list of objects
  TH1D*         hVtxBeforeCuts;     //! Vertex z dist before cuts
  TH1D*         hVtxAfterCuts;      //! Vertex z dist after cuts
  TH1D*         hn1;
  TH1D*         hso;
  TH1D*         hst;
  TH1D*         HStMultRef0;
  TH1D*         HSoMultRef0;
  TH1D*         HMultRef;

  // protected:
  AliTransverseEventShape* fESASelection; // event selection class


  AliAnaTransverseEventShapeTask(const AliAnaTransverseEventShapeTask&);            // not implemented
  AliAnaTransverseEventShapeTask& operator=(const AliAnaTransverseEventShapeTask&); // not implemented

  ClassDef(AliAnaTransverseEventShapeTask, 1);   
};

#endif
