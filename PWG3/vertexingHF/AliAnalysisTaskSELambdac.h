#ifndef ALIANALYSISTASKLAMBDAC_H
#define ALIANALYSISTASKLAMBDAC_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSELambdac
// AliAnalysisTaskSE for the D+ candidates Invariant Mass Histogram and 
//comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// Renu Bala, bala@to.infn.it
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
  
  Float_t GetUpperMassLimit(){return fUpmasslimit;}
  Float_t GetLowerMassLimit(){return fLowmasslimit;}
  Int_t GetNBinsPt(){return fNPtBins;}
  Double_t GetPtBinLimit(Int_t ibin);
  Bool_t IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC);
  Bool_t IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC);
  Bool_t IspiKpReal(AliAODRecoDecayHF3Prong *d);
  Bool_t IspKpiReal(AliAODRecoDecayHF3Prong *d);
  Bool_t IspiKpResonant(AliAODRecoDecayHF3Prong *d,Double_t field);
  Bool_t IspKpiResonant(AliAODRecoDecayHF3Prong *d,Double_t field);
  Bool_t VertexingKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field);
  Int_t MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC);
  Bool_t GetLambdacDaugh(AliAODMCParticle *part, TClonesArray *arrayMC);

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
 
  enum {kMaxPtBins=10};

  TList   *fOutput; //! list send on output slot 0
  TH1F    *fHistNEvents; //!hist. for No. of events
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
  TH2F *fTPCSignal3Sigma;//!hist. for LS inv mass (TC)
  TH2F *fTPCSignal3SigmaReK;//!hist. for LS inv mass (TC)
  TH2F *fTPCSignal3SigmaRep;//!hist. for LS inv mass (TC)
  TH2F *fTPCSignal2Sigma;//!hist. for LS inv mass (TC)
  TH2F *fTPCSignal2SigmaReK;//!hist. for LS inv mass (TC)
  TH2F *fTPCSignal2SigmaRep;//!hist. for LS inv mass (TC)
  TNtuple *fNtupleLambdac; //! output ntuple
  Float_t fUpmasslimit;  //upper inv mass limit for histos
  Float_t fLowmasslimit; //lower inv mass limit for histos
  Float_t fCutsKF[10];
  Int_t fNPtBins; //number of bins in Pt for histograms
  AliRDHFCutsLctopKpi *fRDCutsAnalysis; //Cuts for Analysis
  AliRDHFCutsLctopKpi *fRDCutsProduction; //Production Cuts
  TList *fListCuts; //list of cuts
  Double_t fArrayBinLimits[kMaxPtBins+1]; //limits for the Pt bins
  Bool_t fFillNtuple;   // flag for filling ntuple
  Bool_t fReadMC;    //flag for access to MC
  Bool_t fMCPid;    //flag for access to MC
  Bool_t fRealPid;    //flag for access to MC
  Bool_t fResPid;      //flag to do LS analysis
  Bool_t fUseKF;      //flag to do LS analysis
  AliAnalysisVertexingHF *fVHF;  // Vertexer heavy flavour (used to pass the cuts)
  
  ClassDef(AliAnalysisTaskSELambdac,1); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif

