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
//#include "AliAODpidUtil.h"
#include "AliPIDResponse.h"

class AliAnalysisTaskSELambdac : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSELambdac();
  AliAnalysisTaskSELambdac(const char *name, Bool_t fillNtuple,AliRDHFCutsLctopKpi *lccutsana, AliRDHFCutsLctopKpi *lccutsprod);
  virtual ~AliAnalysisTaskSELambdac();

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetMCPid(){fMCPid=kTRUE;fReadMC=kTRUE;fRealPid=kFALSE;fResPid=kFALSE;return;}
  void SetRealPid(){fRealPid=kTRUE;fMCPid=kFALSE;return;}
  void SetResonantPid(){fResPid=kTRUE;fRealPid=kTRUE;fMCPid=kFALSE;return;}
  void SetCutsKF(Float_t cutsKF[10]){for(Int_t i=0;i<10;i++){fCutsKF[i]=cutsKF[i];}return;}
  void SetUseKF(Bool_t useKF=kTRUE){fUseKF=useKF;}
  void SetAnalysis(Bool_t analysis=kTRUE){fAnalysis=analysis;}
  void SetMassLimits(Float_t range);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetPtBinLimit(Int_t n, Float_t *limitarray);
  void SetFillVarHists(Bool_t setter) {fFillVarHists=setter;return;}
  void SetMultiplicityHists(Bool_t setter) {fMultiplicityHists=setter;return;}
  void SetPriorsHists(Bool_t setter) {fPriorsHists=setter;return;}
  
  Float_t GetUpperMassLimit() const {return fUpmasslimit;}
  Float_t GetLowerMassLimit() const {return fLowmasslimit;}
  Int_t GetNBinsPt() const {return fNPtBins;}
  Double_t GetPtBinLimit(Int_t ibin) const ;
  Bool_t IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Int_t IspiKpResonant(AliAODRecoDecayHF3Prong *d,Double_t field) const ;
  Int_t IspKpiResonant(AliAODRecoDecayHF3Prong *d,Double_t field) const ;
  Bool_t VertexingKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const ;
  Int_t MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t GetLambdacDaugh(AliAODMCParticle *part, TClonesArray *arrayMC) const ;

  void FillMassHists(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC, AliRDHFCutsLctopKpi *cuts);
  void FillVarHists(AliAODRecoDecayHF3Prong *part, TClonesArray *arrMC, AliRDHFCutsLctopKpi *cuts, /*TList *listout,*/ AliAODEvent *aod);
  Bool_t Is3ProngFromPDG(AliAODRecoDecayHF3Prong *part, TClonesArray *arrMC, Int_t pdgToBeCompared=4);
  Bool_t IsTrackFromPDG(AliAODTrack *daugh, TClonesArray *arrayMC, Int_t pdgToBeCompared);
  Bool_t IsThereAGeneratedLc(TClonesArray *arrayMC);
  Int_t NumberPrimaries(AliAODEvent *aods);
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
  void FillAPrioriConcentrations(AliAODRecoDecayHF3Prong *part, AliRDHFCutsLctopKpi *cuts,
				 AliAODEvent* aod, TClonesArray *arrMC);
  void MultiplicityStudies(AliAODRecoDecayHF3Prong *part, AliRDHFCutsLctopKpi *cuts,
			   AliAODEvent* aod, TClonesArray *arrMC,
			   Bool_t &flag1,Bool_t &flag2,Bool_t &flag3,
			   Bool_t &flag4, Bool_t &flag5, Bool_t &flag6); 
  enum {kMaxPtBins=10};

  TList   *fOutput; //! list send on output slot 0
  TH1F    *fHistNEvents; //!hist. for No. of events
  TH1F    *fhChi2; //!hist. for No. of events
  TH1F    *fhMassPtGreater3; //!hist. for No. of events
  TH1F    *fhMassPtGreater3TC; //!hist. for No. of events
  TH1F    *fhMassPtGreater3Kp; //!hist. for No. of events
  TH1F    *fhMassPtGreater3KpTC; //!hist. for No. of events
  TH1F    *fhMassPtGreater3Lpi; //!hist. for No. of events
  TH1F    *fhMassPtGreater3LpiTC; //!hist. for No. of events
  TH1F    *fhMassPtGreater2; //!hist. for No. of events
  TH1F    *fhMassPtGreater2TC; //!hist. for No. of events
  TH1F    *fhMassPtGreater2Kp; //!hist. for No. of events
  TH1F    *fhMassPtGreater2KpTC; //!hist. for No. of events
  TH1F    *fhMassPtGreater2Lpi; //!hist. for No. of events
  TH1F    *fhMassPtGreater2LpiTC; //!hist. for No. of events
  TH1F *fMassHist[3*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fMassHistTC[3*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistLpi[3*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fMassHistLpiTC[3*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistKp[3*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fMassHistKpTC[3*kMaxPtBins]; //!hist. for inv mass (TC)
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
  Bool_t fAnalysis;      //apply analysis cuts
  AliAnalysisVertexingHF *fVHF;  // Vertexer heavy flavour (used to pass the cuts)
  Bool_t fFillVarHists;
  Bool_t fMultiplicityHists;
  Bool_t fPriorsHists;
  TH1F *fNentries;
  TList *fOutputMC;
  TList *fAPriori;
  TList *fMultiplicity;
  //AliAODpidUtil* fUtilPid;
  AliPIDResponse *fPIDResponse;     //! PID response object

  ClassDef(AliAnalysisTaskSELambdac,5); // AliAnalysisTaskSE for the invariant mass analysis of heavy-flavour decay candidates (Lambdac)
};

#endif

