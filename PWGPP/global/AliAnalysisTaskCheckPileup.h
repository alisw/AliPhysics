#ifndef ALIANALYSISTASKPILEUP_H
#define ALIANALYSISTASKPILEUP_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskCheckPileup
// AliAnalysisTask to check pileup tagging performance on ESDs
// Author: F. Prino
//*************************************************************************


class TH1F;

#include "AliAnalysisTaskSE.h"
#include "AliCounterCollection.h"

class AliAnalysisTaskCheckPileup : public AliAnalysisTaskSE 
{
 public:

  AliAnalysisTaskCheckPileup();
  virtual ~AliAnalysisTaskCheckPileup(); 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  Bool_t         GetReadMC() const { return fReadMC; }
  void           SetReadMC(Bool_t flag) {fReadMC=flag;}
  void           SetFillTree(Bool_t flag) {fFillTree=flag;}
  void           SetPhysicsSelectionOption(Int_t opt){fUsePhysSel=opt;}
  void           SetUseAllTriggers(){fTrigOpt=-1;}
  void           SetTriggerMB(){fTrigOpt=0;}
  void           SetTriggerHighMultV0(){fTrigOpt=1;}
  void           SetTriggerHighMultSPD(){fTrigOpt=2;}
  void           SetUsePFProtection(Bool_t opt){fUsePFProtection=opt;}

  AliCounterCollection* GetCounter(){return fCounterPerRun;}

  void SetCutOnContribToSPDPileupVert(Int_t cutc) {fSPDContributorsCut=cutc;}
  void SetCutOnSPDZDiff(Double_t cutz) {fSPDZDiffCut=cutz;}
  void ConfigureMultiTrackVertexPileup(Int_t nc, Double_t wdz, Double_t chi2, Bool_t bcflag){
    fMVContributorsCut=nc;  fMVCChi2Cut=chi2;
    fMVWeiZDiffCut=wdz;  fMVCheckPlpFromDifferentBC=bcflag;
  }
  void SetZDiamondCut(Double_t zd){fZDiamondCut=zd;}

  void SetTriggerMask(Int_t mask) {fTriggerMask = mask;}
  Bool_t ApplyPhysSel(AliESDEvent* esd);

 protected:
  Bool_t fReadMC;           // flag to read Monte Carlo information
  Bool_t fFillTree;      // flag to switch off ntuple
  Int_t  fUsePhysSel;   // phys sel configuration: 0=no, 1= standard, 2=by hand
  Int_t  fTrigOpt;      // trigger selection
  Bool_t fUsePFProtection;    // flag to switch off ntuple
  TList* fOutputPrimV;        //! 1st list of output histos
  TList* fOutputSPDPil;      //! 2nd list of output histos
  TList* fOutputMVPil;       //! 3rd list of output histos
  TList* fOutputMC;          //! 4th list of output histos
  // Primary vertex histos
  TH1F* fHistoXVertSPD;      // histogram
  TH1F* fHistoYVertSPD;      // histogram
  TH1F* fHistoZVertSPD;      // histogram
  TH1F* fHistoXVertTRK;      // histogram
  TH1F* fHistoYVertTRK;      // histogram
  TH1F* fHistoZVertTRK;      // histogram
  // Global histos
  TH2F* fHistoTPCTracksVsTracklets; // histogram
  TH2F* fHistoGloTracksVsTracklets; // histogram
  TH2F* fHistoSPDvertVsFastOr; // histogram
  // SPD Vertex Pileup histos
  TH1F* fHistoNOfPileupVertSPD; // histogram
  TH1F* fHistoNtracklPilSPD;    // histogram
  TH1F* fHistoNtracklNoPilSPD;  // histogram
  TH1F* fHistoNCL1PilSPD;       // histogram
  TH1F* fHistoNCL1NoPilSPD;     // histogram
  TH1F* fHistoContribPrimVertPilSPD;    // histogram
  TH1F* fHistoContribPrimVertNoPilSPD;  // histogram
  TH1F* fHistoContribFirstPilSPD;    // histogram
  TH1F* fHistoZDiffFirstPilSPD;      // histogram
  TH1F* fHistoContribSecondPilSPD;   // histogram
  TH1F* fHistoZDiffSecondPilSPD;     // histogram
  TH1F* fHistoContribTaggingPilSPD;    // histogram
  TH1F* fHistoZDiffTaggingPilSPD;      // histogram
  TH1F* fHistoZDiffTaggingPilZDiamcutSPD;    // histogram
  // Track Vertex Pileup histos
  TH1F* fHistoNOfPileupVertMV; // histogram
  TH1F* fHistoNtracklPilMV;    // histogram
  TH1F* fHistoNtracklNoPilMV;  // histogram
  TH1F* fHistoNCL1PilMV;       // histogram
  TH1F* fHistoNCL1NoPilMV;     // histogram
  TH1F* fHistoContribPrimVertPilMV;    // histogram
  TH1F* fHistoContribPrimVertNoPilMV;  // histogram
  TH1F* fHistoContribFirstPilMV;    // histogram
  TH1F* fHistoZDiffFirstPilMV;      // histogram
  TH1F* fHistoContribSecondPilMV;   // histogram
  TH1F* fHistoZDiffSecondPilMV;     // histogram
  TH1F* fHistoContribTaggingPilMV;    // histogram
  TH1F* fHistoZDiffTaggingPilMV;      // histogram
  TH1F* fHistoZDiffTaggingPilZDiamcutMV;   // histogram
  // MC histos
  TH1F* fHistoNGenCollis; // histogram
  TH1F* fHistoNGenCollisSPDstrobe; // histogram
  TH2F* fHistoCollisTimeOrbit; // histogram
  TH2F* fHistoCollisTimeSPDstrobe; // histogram

  AliCounterCollection* fCounterPerRun; // counters

  TTree* fTrackTree;    // track counters
  UInt_t fTimeStamp;    // tree variables
  UInt_t fNTracksTPC;    // tree variables
  UInt_t fNTracksTPCITS; // tree variables
  UInt_t fNTracklets;    // tree variables
  UInt_t fNContribSPD;   // tree variables
  UInt_t fNFastOr;       // tree variables

  Int_t fSPDContributorsCut;  // cut on cotrtributors to SPD pileup vertex
  Double_t fSPDZDiffCut;      // cut on z diff of SPD pileup vertex

  Int_t  fMVContributorsCut;  // cut on cotrtributors to MV pileup vertex
  Double_t fMVCChi2Cut;       //minimum value of Chi2perNDF of the pileup vertex, multi-vertex
  Float_t  fMVWeiZDiffCut;    //minimum of the sqrt of weighted distance between the primary and the pilup vertex, multi-vertex
  Bool_t  fMVCheckPlpFromDifferentBC; //pileup from different BC (Bunch Crossings)
  Double_t fZDiamondCut; // cut on z of diamond

  Int_t  fTriggerMask;    //flag to select trigger type

 private:    

  AliAnalysisTaskCheckPileup(const AliAnalysisTaskCheckPileup&); // not implemented
  AliAnalysisTaskCheckPileup& operator=(const AliAnalysisTaskCheckPileup&); // not implemented
  
  ClassDef(AliAnalysisTaskCheckPileup,4); // primary vertex analysis
};

#endif
