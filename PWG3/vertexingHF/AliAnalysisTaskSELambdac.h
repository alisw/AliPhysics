#ifndef ALIANALYSISTASKLAMBDAC_H
#define ALIANALYSISTASKLAMBDAC_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSELambdac
// AliAnalysisTaskSE for the Lambdac candidates Invariant Mass Histogram and 
//comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// r.romita@gsi.de
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TArrayD.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODMCParticle.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCuts.h"
#include "TClonesArray.h"

class AliAnalysisTaskSELambdac : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSELambdac();
  AliAnalysisTaskSELambdac(const char *name, Bool_t fillNtuple,AliRDHFCutsLctopKpi *lccutsana, AliRDHFCutsLctopKpi *lccutsprod);
  virtual ~AliAnalysisTaskSELambdac();

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetMCPid(){fMCPid=kTRUE;fReadMC=kTRUE;fRealPid=kFALSE;fResPid=kFALSE;return;}
  void SetRealPid(){fRealPid=kTRUE;fMCPid=kFALSE;fResPid=kFALSE;return;}
  void SetResonantPid(){fResPid=kTRUE;fRealPid=kFALSE;fMCPid=kFALSE;return;}
  void SetCutsKF(Float_t cutsKF[10]){for(Int_t i=0;i<10;i++){fCutsKF[i]=cutsKF[i];}return;}
  void SetUseKF(Bool_t useKF=kTRUE){fUseKF=useKF;}
  void SetMassLimits(Float_t range);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetPtBinLimit(Int_t n, Float_t *limitarray);
  
  Float_t GetUpperMassLimit() const {return fUpmasslimit;}
  Float_t GetLowerMassLimit() const {return fLowmasslimit;}
  Int_t GetNBinsPt() const {return fNPtBins;}
  Double_t GetPtBinLimit(Int_t ibin) const ;
  Bool_t IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t IspiKpResonant(AliAODRecoDecayHF3Prong *d,Double_t field) const ;
  Bool_t IspKpiResonant(AliAODRecoDecayHF3Prong *d,Double_t field) const ;
  Bool_t VertexingKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const ;
  Int_t MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t GetLambdacDaugh(AliAODMCParticle *part, TClonesArray *arrayMC) const ;

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  AliAnalysisTaskSELambdac(const AliAnalysisTaskSELambdac &source);
  AliAnalysisTaskSELambdac& operator=(const AliAnalysisTaskSELambdac& source); 
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*3;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*3+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*3+2;}
  Int_t GetLSHistoIndex(Int_t iPtBin)const { return iPtBin*5;}

  Bool_t ReconstructKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const;
 
  enum {kMaxPtBins=10};

  TList   *fOutput; //! list send on output slot 0
  TH1F    *fHistNEvents; //!hist. for No. of events
  TH1F    *fhChi2; //!hist. for No. of events
  TH1F    *fhMassPtGreater3; //!hist. for No. of events
  TH1F    *fhMassPtGreater3TC; //!hist. for No. of events
  TH1F *fMassHist[3*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F*   fCosPHist[3*kMaxPtBins]; //!hist. for PointingAngle (LC)
  TH1F*   fDLenHist[3*kMaxPtBins]; //!hist. for Dec Length (LC)
  TH1F*   fSumd02Hist[3*kMaxPtBins]; //!hist. for sum d02 (LC)
  TH1F*   fSigVertHist[3*kMaxPtBins]; //!hist. for sigVert (LC)
  TH1F*   fPtMaxHist[3*kMaxPtBins]; //!hist. for Pt Max (LC)
  TH1F*   fDCAHist[3*kMaxPtBins]; //!hist. for DCA (LC)
  TH1F *fMassHistTC[3*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistLS[5*kMaxPtBins];//!hist. for LS inv mass (LC)
  TH1F *fCosPHistLS[3*kMaxPtBins];//!hist. for LS cuts variable 1 (LC)
  TH1F *fDLenHistLS[3*kMaxPtBins];//!hist. for LS cuts variable 2 (LC)
  TH1F *fSumd02HistLS[3*kMaxPtBins];//!hist. for LS cuts variable 3 (LC)
  TH1F *fSigVertHistLS[3*kMaxPtBins];//!hist. for LS cuts variable 4 (LC)
  TH1F *fPtMaxHistLS[3*kMaxPtBins];//!hist. for LS cuts variable 5 (LC)
  TH1F *fDCAHistLS[3*kMaxPtBins];//!hist. for LS cuts variable 6 (LC)
  TH1F *fMassHistLSTC[5*kMaxPtBins];//!hist. for LS inv mass (TC)
  TNtuple *fNtupleLambdac; //! output ntuple
  Float_t fUpmasslimit;  //upper inv mass limit for histos
  Float_t fLowmasslimit; //lower inv mass limit for histos
  Float_t fCutsKF[10]; //cuts with KF vertexer
  Int_t fNPtBins; //number of bins in Pt for histograms
  AliRDHFCutsLctopKpi *fRDCutsAnalysis; //Cuts for Analysis
  AliRDHFCutsLctopKpi *fRDCutsProduction; //Production Cuts
  TList *fListCuts; //list of cuts
  Double_t fArrayBinLimits[kMaxPtBins+1]; //limits for the Pt bins
  Bool_t fFillNtuple;   // flag for filling ntuple
  Bool_t fReadMC;    //flag for access to MC
  Bool_t fMCPid;    //flag for access to MC
  Bool_t fRealPid;    //flag for real PID
  Bool_t fResPid;      //flag for PID with resonant channels
  Bool_t fUseKF;      //flag to cut with KF vertexer
  AliAnalysisVertexingHF *fVHF;  // Vertexer heavy flavour (used to pass the cuts)
  
  ClassDef(AliAnalysisTaskSELambdac,3); // AliAnalysisTaskSE for the invariant mass analysis of heavy-flavour decay candidates (Lambdac)
};

#endif

