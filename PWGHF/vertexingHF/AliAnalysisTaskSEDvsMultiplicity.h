#ifndef ALIANALYSISTASKSEDVSMULTIPLICITY_H
#define ALIANALYSISTASKSEDVSMULTIPLICITY_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSEDvsMultiplicity
// AliAnalysisTaskSE for the D meson vs. multiplcity analysis
// Authors: Renu Bala, Zaida Conesa del Valle, Francesco Prino
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TArrayD.h>
#include <TFile.h>
#include <TRandom.h>
#include <TProfile.h>
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliNormalizationCounter.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliVertexingHFUtils.h"
#include "AliVEvent.h"



class AliAnalysisTaskSEDvsMultiplicity : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEDvsMultiplicity();
  AliAnalysisTaskSEDvsMultiplicity(const char *name, Int_t pdgMeson, AliRDHFCuts* cuts);
  virtual ~AliAnalysisTaskSEDvsMultiplicity();


  void SetMassLimits(Double_t lowlimit, Double_t uplimit);
  void SetMassLimits(Int_t pdg, Double_t range);
  Double_t GetUpperMassLimit() const {return fUpmasslimit;}
  Double_t GetLowerMassLimit() const {return fLowmasslimit;}
  void SetNMassBins(Int_t nbins){fNMassBins=nbins;}
  Int_t GetNMassBins() const{return fNMassBins;}

  void SetImpactParameterBinning(Int_t nbins, Double_t dmin, Double_t dmax){
    fNImpParBins=nbins;
    fLowerImpPar=dmin;
    fHigherImpPar=dmax;
  }

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetUseBit(Bool_t use=kTRUE){fUseBit=use;}
  void SetDoImpactParameterHistos(Bool_t doImp=kTRUE){fDoImpPar=doImp;}


  void SetMultiplVsZProfileLHC10b(TProfile* hprof){
    if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
    fMultEstimatorAvg[0]=new TProfile(*hprof);
  }
  void SetMultiplVsZProfileLHC10c(TProfile* hprof){
    if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
    fMultEstimatorAvg[1]=new TProfile(*hprof);
  }
  void SetMultiplVsZProfileLHC10d(TProfile* hprof){
    if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
    fMultEstimatorAvg[2]=new TProfile(*hprof);
  }
  void SetMultiplVsZProfileLHC10e(TProfile* hprof){
    if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
    fMultEstimatorAvg[3]=new TProfile(*hprof);
  }


  Int_t CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  AliAnalysisTaskSEDvsMultiplicity(const AliAnalysisTaskSEDvsMultiplicity &source);
  AliAnalysisTaskSEDvsMultiplicity& operator=(const AliAnalysisTaskSEDvsMultiplicity& source); 

  TProfile* GetEstimatorHistogram(const AliVEvent *event);
  void CreateImpactParameterHistos();


  TList  *fOutput; //! list send on output slot 1
  TList  *fListCuts; //list of cuts
  TList  *fOutputCounters; //! list send on output slot 3
  TList  *fListProfiles; //list of profile histos for z-vtx correction

  TH1F *fHistNEvents;     //!hist. for No. of events

  TH2F* fHistNtrEta16vsNtrEta1; //!hist. for Ntracklets in eta<1.6 vs. eta<1.
  TH2F* fHistNtrCorrEta1vsNtrRawEta1; //!hist. for Ntracklets in eta<1 with and w/o corrections 
  TH2F* fHistNtrVsZvtx; //!  hist of ntracklets vs Zvertex
  TH2F* fHistNtrCorrVsZvtx; //!  hist of ntracklets vs Zvertex

  TH2F* fHistNtrVsNchMC; //!  hist of ntracklets vs Nch (Generated)
  TH2F* fHistNtrCorrVsNchMC; //!  hist of ntracklets vs Nch (Generated)
  TH2F* fHistNtrVsNchMCPrimary; //!  hist of ntracklets vs Nch (Primary)
  TH2F* fHistNtrCorrVsNchMCPrimary; //!  hist of ntracklets vs Nch (Primary)
  TH2F* fHistNtrVsNchMCPhysicalPrimary; //!  hist of ntracklets vs Nch (Physical Primary)
  TH2F* fHistNtrCorrVsNchMCPhysicalPrimary; //!  hist of ntracklets vs Nch (Physical Primary)
  TH1F* fHistGenPrimaryParticlesInelGt0; //!hist. of geenrated multiplcity
  
  TH1F* fHistNtrCorrEvSel; //! hist. of ntracklets for selected events
  TH1F* fHistNtrCorrEvWithCand; //! hist. of ntracklets for evnts with a candidate
  TH1F* fHistNtrCorrEvWithD;//! hist. of ntracklets for evnts with a candidate in D mass peak


  TH3F *fPtVsMassVsMult;  //! hist. of Pt vs Mult vs. mass (
  TH3F *fPtVsMassVsMultNoPid;  //! hist. of Pt vs Mult vs. mass (no pid)
  TH3F *fPtVsMassVsMultUncorr;  //! hist. of Pt vs Mult vs. mass (raw mult)
  TH3F *fPtVsMassVsMultPart;  //! hist. of Pt vs Mult vs. mass (particle)
  TH3F *fPtVsMassVsMultAntiPart;  //! hist. of Pt vs Mult vs. mass (antiparticle)

  THnSparseF *fHistMassPtImpPar[5];//! histograms for impact paramter studies

  Double_t fUpmasslimit;  //upper inv mass limit for histos
  Double_t fLowmasslimit; //lower inv mass limit for histos
  Int_t   fNMassBins;    // nbins for invariant mass histos

  AliRDHFCuts *fRDCutsAnalysis; // Cuts for Analysis
  AliNormalizationCounter *fCounter;  //!Counter for normalization
  AliNormalizationCounter *fCounterU; //!Counter for normalization, uncorr mult

  Bool_t fDoImpPar;  //swicth for D impact parameter THnSparse
  Int_t  fNImpParBins;   // nunber of bins in impact parameter histos
  Double_t fLowerImpPar;  // lower limit in impact parameter (um)
  Double_t fHigherImpPar; // higher limit in impact parameter (um)

  Bool_t fReadMC;    //flag for access to MC
  Int_t  fMCOption;  // 0=keep all cand, 1=keep only signal, 2= keep only back
  Bool_t fUseBit;    // flag to use bitmask
  
  TProfile* fMultEstimatorAvg[4]; // TProfile with mult vs. Z per period
  Double_t fRefMult;   // refrence multiplcity (period b)
  Int_t fPdgMeson;   // pdg code of analyzed meson

  
  ClassDef(AliAnalysisTaskSEDvsMultiplicity,3); // D vs. mult task
};

#endif
