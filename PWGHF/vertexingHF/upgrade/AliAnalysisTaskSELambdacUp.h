#ifndef ALIANALYSISTASKSELAMBDACUP_H
#define ALIANALYSISTASKSELAMBDACUP_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//*************************************************************************
// Class AliAnalysisTaskSELambdac
// AliAnalysisTaskSE for the Lambdac candidates Invariant Mass Histogram and 
//comparison of heavy-flavour decay candidates
// Specific for the ITS upgrade 
// r.romita@liv.ac.uk
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
//#include "AliAODpidUtil.h"
#include "AliPIDResponse.h"
#include "AliNormalizationCounter.h"

class AliAnalysisTaskSELambdacUp : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSELambdacUp();
  AliAnalysisTaskSELambdacUp(const char *name, Bool_t fillNtuple,AliRDHFCutsLctopKpi *lccutsana, AliRDHFCutsLctopKpi *lccutsprod);
  virtual ~AliAnalysisTaskSELambdacUp();

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetMCPid(){fMCPid=kTRUE;fReadMC=kTRUE;fRealPid=kFALSE;fResPid=kFALSE;return;}
  void SetRealPid(){fRealPid=kTRUE;fMCPid=kFALSE;return;}
  void SetResonantPid(){fResPid=kTRUE;fRealPid=kTRUE;fMCPid=kFALSE;return;}
  void SetCutsKF(Float_t cutsKF[2]){for(Int_t i=0;i<2;i++){fCutsKF[i]=cutsKF[i];}return;}
  void SetUseKF(Bool_t useKF=kTRUE){fUseKF=useKF;}
  void SetAnalysis(Bool_t analysis=kTRUE){fAnalysis=analysis;}
  void SetMassLimits(Float_t range);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetPtBinLimit(Int_t n, Float_t *limitarray);
  void SetFillVarHists(Bool_t setter) {fFillVarHists=setter;return;}
  
  Float_t GetUpperMassLimit() const {return fUpmasslimit;}
  Float_t GetLowerMassLimit() const {return fLowmasslimit;}
  Int_t GetNBinsPt() const {return fNPtBins;}
  Double_t GetPtBinLimit(Int_t ibin) const ;
  Bool_t IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  void IspiKpResonant(AliAODRecoDecayHF3Prong *d,Double_t field,Int_t *resNumber) const ;
  void IspKpiResonant(AliAODRecoDecayHF3Prong *d,Double_t field,Int_t *resNumber) const ;
  Bool_t VertexingKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const ;
  Int_t MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t GetLambdacDaugh(AliAODMCParticle *part, TClonesArray *arrayMC) const ;
  Bool_t CheckOrigin(AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC);

  void FillMassHists(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC, AliRDHFCutsLctopKpi *cuts,Int_t *nSelectedloose,Int_t *nSelectedtight);
  void FillVarHists(AliAODRecoDecayHF3Prong *part, TClonesArray *arrMC, AliRDHFCutsLctopKpi *cuts, /*TList *listout,*/ AliAODEvent *aod);
 
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  AliAnalysisTaskSELambdacUp(const AliAnalysisTaskSELambdacUp &source);
  AliAnalysisTaskSELambdacUp& operator=(const AliAnalysisTaskSELambdacUp& source); 
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*4;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*4+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*4+2;}
  Int_t GetPYTHIAHistoIndex(Int_t iPtBin)const { return iPtBin*4+3;}

  Bool_t ReconstructKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const;
 
  enum {kMaxPtBins=10};

  TList   *fOutput; //! list send on output slot 0
  TH1F    *fHistNEvents; //!hist. for No. of events
  TH1F    *fhChi2; //!hist. for No. of events
  TH1F *fMassHist[4*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fMassHistTC[4*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistLpi[4*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fMassHistLpiTC[4*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistKp[4*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fMassHistKpTC[4*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistDk[4*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fMassHistDkTC[4*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHist3Pr[4*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fMassHist3PrTC[4*kMaxPtBins]; //!hist. for inv mass (TC)
  TNtuple *fNtupleLambdac; //! output ntuple
  Float_t fUpmasslimit;  //upper inv mass limit for histos
  Float_t fLowmasslimit; //lower inv mass limit for histos
  Float_t fCutsKF[2]; //cuts with KF vertexer
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
  Bool_t fAnalysis;      //apply analysis cuts
  AliAnalysisVertexingHF *fVHF;  // Vertexer heavy flavour (used to pass the cuts)
  Bool_t fFillVarHists;  // flag for creation and fill of histograms with vars
  TH1F *fNentries;      // histo with number of entries
  TList *fOutputMC;     // output1
  //AliAODpidUtil* fUtilPid;
  AliPIDResponse *fPIDResponse;     //! PID response object
  AliNormalizationCounter *fCounter;//!AliNormalizationCounter on output slot 7

  ClassDef(AliAnalysisTaskSELambdacUp,1); // AliAnalysisTaskSE for the invariant mass analysis of heavy-flavour decay candidates (Lambdac)
};

#endif

