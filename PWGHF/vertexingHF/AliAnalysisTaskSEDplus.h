#ifndef ALIANALYSISTASKSEDPLUS_H
#define ALIANALYSISTASKSEDPLUS_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
/// \class Class AliAnalysisTaskSEDplus
/// \brief AliAnalysisTaskSE for the D+ candidates Invariant Mass Histogram and
/// comparison of heavy-flavour decay candidates
/// to MC truth (kinematics stored in the AOD)
/// \author Renu Bala, bala@to.infn.it
/// \author F. Prino, prino@to.infn.it
/// \author G. Ortona, ortona@to.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
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
  AliAnalysisTaskSEDplus(const char *name, AliRDHFCutsDplustoKpipi* analysiscuts,Int_t fillNtuple=0);
  virtual ~AliAnalysisTaskSEDplus();

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetDoLikeSign(Int_t dols=0){fDoLS=dols;}
  void SetSystem(Int_t system=0){fSystem=system;}
  void SetCutsDistr(Bool_t cutsDistr=kTRUE){fCutsDistr=cutsDistr;}
  void SetDoImpactParameterHistos(Bool_t doImp=kTRUE){fDoImpPar=doImp;}
  void SetDoTrackVarHistos(Bool_t doTrackHist=kTRUE){fDoTrackVarHist=doTrackHist;}
  void SetDoMCAcceptanceHistos(Bool_t doMCAcc=kTRUE){fStepMCAcc=doMCAcc;}
  void SetImpactParameterBinning(Int_t nbins, Float_t dmin, Float_t dmax){
    fNImpParBins=nbins;
    fLowerImpPar=dmin;
    fHigherImpPar=dmax;
  }
  void SetUseStrangeness(Bool_t uses=kTRUE){fUseStrangeness=uses;}
  void SetMassLimits(Float_t range);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetBinWidth(Float_t w);
  void SetUseBit(Bool_t dols=kTRUE){fUseBit=dols;}

  void SetUseOnlyPositiveEta(){fEtaSelection=1;}
  void SetUseOnlyNegativeEta(){fEtaSelection=-1;}
  void SetUseFullEta(){fEtaSelection=0;}
  void SetUseQuarkLevelTag(Bool_t opt){fUseQuarkTagInKine=opt;}

  Float_t GetUpperMassLimit(){return fUpmasslimit;}
  Float_t GetLowerMassLimit(){return fLowmasslimit;}
  Int_t GetNBinsPt(){return fNPtBins;}
  Float_t GetBinWidth(){return fBinWidth;}
  Int_t GetNBinsHistos();
  
  void LSAnalysis(TClonesArray *arrayOppositeSign,TClonesArray *arrayLikeSign,AliAODEvent *aod,AliAODVertex *vtx1, Int_t nDplusOS);

  void CreateLikeSignHistos();
  void CreateImpactParameterHistos();
  void CreateTrackVarHistos();
  void CreateMCAcceptanceHistos();

  Bool_t CheckAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau);
  void FillMCAcceptanceHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader);

  /// Implementation of interface methods
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
  Float_t GetStrangenessWeights(const AliAODRecoDecayHF3Prong* d, TClonesArray* arrayMC, Float_t factor[3]) const;

  enum {kMaxPtBins=20};
  enum {kVarForSparse=12,kVarForSparseFD=13,kVarForTrackSparse=7};

  TList   *fOutput; //!<! list send on output slot 0
  TH1F *fHistNEvents; //!<!hist. for No. of events
  TH1F *fMassHistNoPid[3*kMaxPtBins]; //!<!hist. for inv mass (w/o PID)
  TH1F *fCosPHist[3*kMaxPtBins]; //!<!hist. for PointingAngle (topol+PID)
  TH1F *fDLenHist[3*kMaxPtBins]; //!<!hist. for Dec Length (topol+PID)
  TH1F *fSumd02Hist[3*kMaxPtBins]; //!<!hist. for sum d02 (topol+PID)
  TH1F *fSigVertHist[3*kMaxPtBins]; //!<!hist. for sigVert (topol+PID)
  TH1F *fPtMaxHist[3*kMaxPtBins]; //!<!hist. for Pt Max (topol+PID)
  TH1F *fPtKHist[3*kMaxPtBins]; //!<!hist. for PtK (topol+PID)
  TH1F *fPtpi1Hist[3*kMaxPtBins]; //!<!hist. for PtPi1 (topol+PID)
  TH1F *fPtpi2Hist[3*kMaxPtBins]; //!<!hist. for PtPi2 (topol+PID)
  TH1F *fDCAHist[3*kMaxPtBins]; //!<!hist. for DCA (topol+PID)
  TH1F *fDLxy[3*kMaxPtBins]; //!<!hist. for DLxy (topol+PID)
  TH1F *fCosxy[3*kMaxPtBins]; //!<!hist. for Cosxy (topol+PID)
  TH1F *fMassHist[3*kMaxPtBins]; //!<!hist. for inv mass (topol+PID cuts)
  TH1F *fMassHistPlus[3*kMaxPtBins]; //!<!hist. for D+ inv mass (topol+PID cuts)
  TH1F *fMassHistMinus[3*kMaxPtBins]; //!<!hist. for D- inv mass (topol+PID cuts)
  TH1F *fMassHistLS[5*kMaxPtBins];//!<!hist. for LS inv mass (topol+PID)
  TH1F *fCosPHistLS[3*kMaxPtBins];//!<!hist. for LS cuts variable 1 (topol+PID)
  TH1F *fDLenHistLS[3*kMaxPtBins];//!<!hist. for LS cuts variable 2 (topol+PID)
  TH1F *fSumd02HistLS[3*kMaxPtBins];//!<!hist. for LS cuts variable 3 (topol+PID)
  TH1F *fSigVertHistLS[3*kMaxPtBins];//!<!hist. for LS cuts variable 4 (topol+PID)
  TH1F *fPtMaxHistLS[3*kMaxPtBins];//!<!hist. for LS cuts variable 5 (topol+PID)
  TH1F *fDCAHistLS[3*kMaxPtBins];//!<!hist. for LS cuts variable 6 (topol+PID)
  TH2F *fCorreld0Kd0pi[3]; //!<!hist. for d0k*d0pi vs. d0k*d0pi (topol+PID)
  TH2F *fHistCentrality[3];//!<!hist. for cent distr (all,sel ev, )
  THnSparseF *fHistMassPtImpPar[5];//!<! histograms for impact parameter and cut variation study
  THnSparseF *fHistTrackVar; //!<! histograms for track cuts study
  THnSparseF *fMCAccPrompt; //!<!histo for StepMCAcc for Dplus prompt (pt,y,ptB)
  THnSparseF *fMCAccBFeed; //!<!histo for StepMCAcc for Dplus FD (pt,y,ptB)
  TH2F *fPtVsMassNoPid;    //!<! hist. of pt vs. mass (w/o PID)
  TH2F *fPtVsMass;  //!<! hist. of pt vs. mass (topol+PID cuts)
  TH3F *fYVsPtNoPid;       //!<! hist. of Y vs. Pt vs. Mass(w/o PID)
  TH3F *fYVsPt;     //!<! hist. of Y vs. Pt vs. Mass (topol+PID cuts)
  TH2F *fYVsPtSigNoPid;    //!<! hist. of Y vs. Pt (MC, only sig, w/o PID)
  TH2F *fYVsPtSig;    //!<! hist. of Y vs. Pt (MC, only sig, topol+PID cuts)
  TH2F *fPhiEtaCand;      //!<! hist. with eta/phi distribution of candidates
  TH2F *fPhiEtaCandSigReg;//!<! hist. eta/phi of candidates in D+ mass region
  TH1F *fSPDMult;    //!<! hist. of spd mult
  TNtuple *fNtupleDplus; //!<! output ntuple
  Float_t fUpmasslimit;  /// upper inv mass limit for histos
  Float_t fLowmasslimit; /// lower inv mass limit for histos
  Int_t fNPtBins; /// Number of Pt Bins
  Float_t fBinWidth;/// width of one bin in output histos
  TList *fListCuts; /// list of cuts
  AliRDHFCutsDplustoKpipi *fRDCutsAnalysis; /// Cuts for Analysis
  AliNormalizationCounter *fCounter;//!<!Counter for normalization
  Double_t fArrayBinLimits[kMaxPtBins+1]; /// limits for the Pt bins
  Int_t fFillNtuple;   /// flag for filling ntuple 0 no NTuple 1 big Ntuple 2 small NTuple
  Bool_t fReadMC;    /// flag for access to MC
  Bool_t fUseStrangeness;/// flag to enhance strangeness in MC to fit to data
  Bool_t fUseBit;      /// flag to use bitmask
  Bool_t fCutsDistr;    /// flag to activate cuts distr histos
  Bool_t fDoImpPar;    /// flag to activate impact paramter histos
  Bool_t fDoTrackVarHist; ///flag to activate track variable cut studies
  Bool_t fStepMCAcc;   /// flag to activate histos for StepMCAcc
  Bool_t fUseQuarkTagInKine; /// flag for quark/hadron level identification of prompt and feeddown
  Int_t  fNImpParBins;   /// nunber of bins in impact parameter histos
  Float_t fLowerImpPar;  /// lower limit in impact parameter (um)
  Float_t fHigherImpPar; /// higher limit in impact parameter (um)
  Int_t  fDoLS;        /// flag to do LS analysis
  Int_t fEtaSelection; /// eta region to accept D+ 0=all, -1 = negative, 1 = positive
  Int_t fSystem;   /// 0=pp,1=PbPb
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEDplus,26); /// AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
  /// \endcond
};

#endif
