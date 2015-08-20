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
  void           SetFillNtuple(Bool_t flag) {fFillNtuple=flag;}
  
  AliCounterCollection* GetCounter(){return fCounterPerRun;}

  void SetCutOnContribToSPDPileupVert(Int_t cutc) {fSPDContributorsCut=cutc;}
  void SetCutOnSPDZDiff(Double_t cutz) {fSPDZDiffCut=cutz;}
  void ConfigureMultiTrackVertexPileup(Int_t nc, Double_t wdz, Double_t chi2, Bool_t bcflag){
    fMVContributorsCut=nc;  fMVCChi2Cut=chi2;
    fMVWeiZDiffCut=wdz;  fMVCheckPlpFromDifferentBC=bcflag;
  }
  void SetZDiamondCut(Double_t zd){fZDiamondCut=zd;}

  void SetTriggerMask(Int_t mask) {fTriggerMask = mask;}

 protected:
  Bool_t fReadMC;           // flag to read Monte Carlo information
  Bool_t fFillNtuple;      // flag to switch off ntuple

  TList* fOutputPrimV;        //! 1st list of output histos
  TList* fOutputSPDPil;      //! 2nd list of output histos
  TList* fOutputMVPil;       //! 3rd list of output histos
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

  AliCounterCollection* fCounterPerRun; // counters

  TNtuple* fTrackNtuple; // track counters

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
  
  ClassDef(AliAnalysisTaskCheckPileup,2); // primary vertex analysis
};

#endif
