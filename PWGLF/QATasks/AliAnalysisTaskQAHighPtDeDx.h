#ifndef ALIANALYSISTASKQAHIGHPTDEDX_H
#define ALIANALYSISTASKQAHIGHPTDEDX_H
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



class AliAnalysisTaskQAHighPtDeDx : public AliAnalysisTaskSE {
 public:
 

  AliAnalysisTaskQAHighPtDeDx();
  AliAnalysisTaskQAHighPtDeDx(const char *name);
  virtual ~AliAnalysisTaskQAHighPtDeDx();






  //AliAnalysisTaskQAHighPtDeDx(const char *name="<default name>");
  //virtual ~AliAnalysisTaskQAHighPtDeDx() { /*if (fOutputList) delete fOutputList;*/}//;

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  Bool_t   GetAnalysisMC() { return fAnalysisMC; }   
  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }     
  //Double_t GetMinPt() { return fMinPt; }   
  //Int_t    GetTreeOption() { return fTreeOption; }  

  virtual void  SetTrigger(UInt_t ktriggerInt) {ftrigBit = ktriggerInt;}
  virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
  virtual void  SetTrackFilterTPC(AliAnalysisFilter* trackF) {fTrackFilterTPC = trackF;}
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
  virtual void AnalyzeESD(AliESDEvent* esd); 
  virtual void AnalyzeAOD(AliAODEvent* aod); 
  virtual void ProduceArrayTrksESD(AliESDEvent* event);
  virtual void ProduceArrayV0ESD(AliESDEvent* event);
  virtual void ProduceArrayTrksAOD(AliAODEvent* event);
  virtual void ProduceArrayV0AOD(AliAODEvent* event);
  Short_t   GetPidCode(Int_t pdgCode) const;
  void      ProcessMCTruthESD();
  void      ProcessMCTruthAOD(); 

  Short_t   GetPythiaEventProcessType(Int_t pythiaType);
  Short_t   GetDPMjetEventProcessType(Int_t dpmJetType);
  ULong64_t GetEventIdAsLong(AliVHeader* header) const;

  TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
  Int_t      FindPrimaryMotherLabel(AliStack* stack, Int_t label);

  AliAODMCParticle* FindPrimaryMotherAOD(AliAODMCParticle* startParticle);

  TParticle* FindPrimaryMotherV0(AliStack* stack, Int_t label);
  Int_t      FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps);
  Bool_t PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh);


  AliAODMCParticle* FindPrimaryMotherAODV0(AliAODMCParticle* startParticle, Int_t& nSteps);



  static const Double_t fgkClight;   // Speed of light (cm/ps)

  AliESDEvent* fESD;                  //! ESD object
  AliAODEvent* fAOD;                  //! AOD object
  AliMCEvent*  fMC;                   //! MC object
  AliStack*    fMCStack;              //! MC ESD stack
  TClonesArray* fMCArray;             //! MC array for AOD
  AliAnalysisFilter* fTrackFilterGolden;    //  Track Filter, set 2010 with golden cuts
  AliAnalysisFilter* fTrackFilterTPC; // track filter for TPC only tracks
  TString       fCentEst;             // V0A , V0M, 
  TString       fAnalysisType;        //  "ESD" or "AOD"
  Bool_t        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
  Bool_t        fAnalysisPbPb;        //  true you want to analyze PbPb data, false for pp
  UInt_t       ftrigBit;
  TRandom*      fRandom;              //! random number generator
  Bool_t        fPileUpRej;           // kTRUE is pile-up is rejected
 


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
  TH1I*         fEvents;            //! No of accepted events
  TH1I*         fVtx;               //! Event vertex info
  TH1F*         fVtxMC;             //! Event vertex info for ALL MC events
  TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
  TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts
  TH1F* fn1;
  TH1F* fcent;

  TH2D *hMIPVsEta;
  TProfile *pMIPVsEta;
  TH2D *hMIPVsEtaV0s;
  TProfile *pMIPVsEtaV0s;
  TH2D *hPlateauVsEta;
  TProfile *pPlateauVsEta;
  TH2D *hPhi;

  TH2D     *hMIPVsNch[9];
  TProfile *pMIPVsNch[9];

  TH2D     *hMIPVsPhi[9];
  TProfile *pMIPVsPhi[9];
  TH2D     *hPlateauVsPhi[9];
  TProfile *pPlateauVsPhi[9];

  TH2D* histPiV0[9];
  TH1D* histpPiV0[9];

  TH2D* histPV0[9];
  TH1D* histpPV0[9];

  TH2D* histAllCh[9];

  TH2D* histPiTof[9];
  TH1D* histpPiTof[9];

  TH2D* histEV0[9];

  TH1D* hMcIn[7][9];
  TH1D* hMcOut[7][9];

  AliAnalysisTaskQAHighPtDeDx(const AliAnalysisTaskQAHighPtDeDx&);            // not implemented
  AliAnalysisTaskQAHighPtDeDx& operator=(const AliAnalysisTaskQAHighPtDeDx&); // not implemented

  //TTree*        fTree;              //! Debug tree 

  ClassDef(AliAnalysisTaskQAHighPtDeDx, 1);    //Analysis task for high pt analysis 
};

#endif
