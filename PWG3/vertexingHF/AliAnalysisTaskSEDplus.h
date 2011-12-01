#ifndef ALIANALYSISTASKSEDPLUS_H
#define ALIANALYSISTASKSEDPLUS_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSEDplus
// AliAnalysisTaskSE for the D+ candidates Invariant Mass Histogram and 
//comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// Renu Bala, bala@to.infn.it
// F. Prino, prino@to.infn.it
// G. Ortona, ortona@to.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TArrayD.h>

#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliNormalizationCounter.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

class AliAnalysisTaskSEDplus : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEDplus();
  AliAnalysisTaskSEDplus(const char *name, AliRDHFCutsDplustoKpipi* analysiscuts,AliRDHFCutsDplustoKpipi* productioncuts,Bool_t fillNtuple=kFALSE);
  virtual ~AliAnalysisTaskSEDplus();

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetDoLikeSign(Int_t dols=0){fDoLS=dols;}
  void SetCutsDistr(Bool_t cutsDistr=kTRUE){fCutsDistr=cutsDistr;}
  void SetDoImpactParameterHistos(Bool_t doImp=kTRUE){fDoImpPar=doImp;}
  void SetImpactParameterBinning(Int_t nbins, Float_t dmin, Float_t dmax){
    fNImpParBins=nbins;
    fLowerImpPar=dmin;
    fHigherImpPar=dmax;
  }
  void SetUseStrangeness(Bool_t uses=kTRUE){fUseStrangeness=uses;}
  void SetMassLimits(Float_t range);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetPtBinLimit(Int_t n, Float_t *limitarray);
  void SetBinWidth(Float_t w);
  void SetUseBit(Bool_t dols=kTRUE){fUseBit=dols;}


  
  Float_t GetUpperMassLimit(){return fUpmasslimit;}
  Float_t GetLowerMassLimit(){return fLowmasslimit;}
  Int_t GetNBinsPt(){return fNPtBins;}
  Double_t GetPtBinLimit(Int_t ibin);
  Float_t GetBinWidth(){return fBinWidth;}
  Int_t GetNBinsHistos();
  
  void LSAnalysis(TClonesArray *arrayOppositeSign,TClonesArray *arrayLikeSign,AliAODEvent *aod,AliAODVertex *vtx1, Int_t nDplusOS);
  Int_t CheckOrigin(TClonesArray* arrayMC, const AliAODMCParticle *mcPartCandidate) const;
  void CreateLikeSignHistos();
  void CreateImpactParameterHistos();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  AliAnalysisTaskSEDplus(const AliAnalysisTaskSEDplus &source);
  AliAnalysisTaskSEDplus& operator=(const AliAnalysisTaskSEDplus& source); 
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*3;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*3+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*3+2;}
  Int_t GetLSHistoIndex(Int_t iPtBin)const { return iPtBin*5;}
  Float_t GetTrueImpactParameter(const AliAODMCHeader *mcHeader, TClonesArray* arrayMC, const AliAODMCParticle *partDp) const;

  enum {kMaxPtBins=20};

  TList   *fOutput; //! list send on output slot 0
  TH1F *fHistNEvents; //!hist. for No. of events
  TH1F *fMassHist[3*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fCosPHist[3*kMaxPtBins]; //!hist. for PointingAngle (LC)
  TH1F *fDLenHist[3*kMaxPtBins]; //!hist. for Dec Length (LC)
  TH1F *fSumd02Hist[3*kMaxPtBins]; //!hist. for sum d02 (LC)
  TH1F *fSigVertHist[3*kMaxPtBins]; //!hist. for sigVert (LC)
  TH1F *fPtMaxHist[3*kMaxPtBins]; //!hist. for Pt Max (LC)
  TH1F *fPtKHist[3*kMaxPtBins]; //!hist. for PtK (LC)
  TH1F *fPtpi1Hist[3*kMaxPtBins]; //!hist. for PtPi1 (LC)
  TH1F *fPtpi2Hist[3*kMaxPtBins]; //!hist. for PtPi2 (LC)
  TH1F *fDCAHist[3*kMaxPtBins]; //!hist. for DCA (LC)
  TH1F *fDLxy[3*kMaxPtBins]; //!hist. for DLxy (LC)
  TH1F *fDLxyTC[3*kMaxPtBins]; //!hist. for DLxy (TC)
  TH1F *fCosxy[3*kMaxPtBins]; //!hist. for Cosxy (LC)
  TH1F *fCosxyTC[3*kMaxPtBins]; //!hist. for Cosxy (TC)
  TH1F *fMassHistTC[3*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistTCPlus[3*kMaxPtBins]; //!hist. for D+ inv mass (TC)
  TH1F *fMassHistTCMinus[3*kMaxPtBins]; //!hist. for D- inv mass (TC)
  TH1F *fMassHistLS[5*kMaxPtBins];//!hist. for LS inv mass (LC)
  TH1F *fCosPHistLS[3*kMaxPtBins];//!hist. for LS cuts variable 1 (LC)
  TH1F *fDLenHistLS[3*kMaxPtBins];//!hist. for LS cuts variable 2 (LC)
  TH1F *fSumd02HistLS[3*kMaxPtBins];//!hist. for LS cuts variable 3 (LC)
  TH1F *fSigVertHistLS[3*kMaxPtBins];//!hist. for LS cuts variable 4 (LC)
  TH1F *fPtMaxHistLS[3*kMaxPtBins];//!hist. for LS cuts variable 5 (LC)
  TH1F *fDCAHistLS[3*kMaxPtBins];//!hist. for LS cuts variable 6 (LC)
  TH1F *fMassHistLSTC[5*kMaxPtBins];//!hist. for LS inv mass (TC)
  TH2F *fCorreld0Kd0pi[3]; //!hist. for d0k*d0pi vs. d0k*d0pi (LC)
  TH1F *fHistCentrality[3];//!hist. for cent distr (all,sel ev, )
  THnSparseF *fHistMassPtImpParTC[5];//! histograms for impact paramter studies
  TH2F *fPtVsMass;    //! hist. of pt vs. mass (prod. cuts)
  TH2F *fPtVsMassTC;  //! hist. of pt vs. mass (analysis cuts)
  TH2F *fYVsPt;       //! hist. of Y vs. Pt (prod. cuts)
  TH2F *fYVsPtTC;     //! hist. of Y vs. Pt (analysis cuts)
  TH2F *fYVsPtSig;    //! hist. of Y vs. Pt (MC, only sig, prod. cuts)
  TH2F *fYVsPtSigTC;    //! hist. of Y vs. Pt (MC, only sig, analysis cuts)
  TNtuple *fNtupleDplus; //! output ntuple
  Float_t fUpmasslimit;  //upper inv mass limit for histos
  Float_t fLowmasslimit; //lower inv mass limit for histos
  Int_t fNPtBins; //Number of Pt Bins
  Float_t fBinWidth;//width of one bin in output histos
  TList *fListCuts; //list of cuts
  AliRDHFCutsDplustoKpipi *fRDCutsProduction; //Production D+ Cuts
  AliRDHFCutsDplustoKpipi *fRDCutsAnalysis; //Cuts for Analysis
  AliNormalizationCounter *fCounter;//!Counter for normalization
  Double_t fArrayBinLimits[kMaxPtBins+1]; //limits for the Pt bins
  Bool_t fFillNtuple;   // flag for filling ntuple
  Bool_t fReadMC;    //flag for access to MC
  Bool_t fUseStrangeness;//flag to enhance strangeness in MC to fit to data
  Bool_t fUseBit;      // flag to use bitmask
  Bool_t fCutsDistr;    // flag to activate cuts distr histos
  Bool_t fDoImpPar;    // flag to activate impact paramter histos
  Int_t  fNImpParBins;   // nunber of bins in impact parameter histos
  Float_t fLowerImpPar;  // lower limit in impact parameter (um)
  Float_t fHigherImpPar; // higher limit in impact parameter (um)
  Int_t  fDoLS;        // flag to do LS analysis
  
  ClassDef(AliAnalysisTaskSEDplus,16); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif
