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
  void SetDoCutVarsSparses(Bool_t doSparse=kTRUE){fDoSparse=doSparse;}
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
  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}

  void SetUseOnlyPositiveEta(){fEtaSelection=1;}
  void SetUseOnlyNegativeEta(){fEtaSelection=-1;}
  void SetUseFullEta(){fEtaSelection=0;}
  void SetUseQuarkLevelTag(Bool_t opt){fUseQuarkTagInKine=opt;}
  
  void SetCutOnNtracklets(Bool_t applycut=kTRUE, Int_t Ntrckmin=0, Int_t Ntrckmax=100) {
    fCutOnTrckl=applycut;
    fNtrcklMin=Ntrckmin;
    fNtrcklMax=Ntrckmax;
  }
  
  Float_t GetUpperMassLimit(){return fUpmasslimit;}
  Float_t GetLowerMassLimit(){return fLowmasslimit;}
  Int_t GetNBinsPt(){return fNPtBins;}
  Float_t GetBinWidth(){return fBinWidth;}
  Int_t GetNBinsHistos();

  void LSAnalysis(TClonesArray *arrayOppositeSign,TClonesArray *arrayLikeSign,AliAODEvent *aod,AliAODVertex *vtx1, Int_t nDplusOS);

  void CreateLikeSignHistos();
  void CreateImpactParameterHistos();
  void CreateCutVarsSparses();
  void CreateTrackVarHistos();
  void CreateMCAcceptanceHistos();

  Bool_t CheckAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau);
  void FillMCAcceptanceHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader, Int_t tracklets);

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

  enum {kVarForSparse=13,kVarForSparseFD=14,kVarForTrackSparse=7,kVarForImpPar=3};

  TList* fOutput;            //!<! list send on output slot 0
  TH1F* fHistNEvents;        //!<! hist. for No. of events
  TH1F* fHistNCandidates;    //!<! hist. for No. of candidates
  TH1F** fMassHist;          //!<! hist. for inv mass (topol+PID cuts)
  TH1F** fMassHistPlus;      //!<! hist. for D+ inv mass (topol+PID cuts)
  TH1F** fMassHistMinus;     //!<! hist. for D- inv mass (topol+PID cuts)
  TH1F** fMassHistNoPid;     //!<! hist. for inv mass (w/o PID)
  TH1F** fCosPHist;          //!<!hist. for PointingAngle (topol+PID)
  TH1F** fDLenHist;          //!<!hist. for Dec Length (topol+PID)
  TH1F** fSumd02Hist;        //!<!hist. for sum d02 (topol+PID)
  TH1F** fSigVertHist;       //!<!hist. for sigVert (topol+PID)
  TH1F** fPtMaxHist;         //!<!hist. for Pt Max (topol+PID)
  TH1F** fPtKHist;           //!<!hist. for PtK (topol+PID)
  TH1F** fPtpi1Hist;         //!<!hist. for PtPi1 (topol+PID)
  TH1F** fPtpi2Hist;         //!<!hist. for PtPi2 (topol+PID)
  TH1F** fDCAHist;           //!<!hist. for DCA (topol+PID)
  TH1F** fDLxy;              //!<!hist. for DLxy (topol+PID)
  TH1F** fCosxy;             //!<!hist. for Cosxy (topol+PID)
  TH2F *fCorreld0Kd0pi[3];   //!<!hist. for d0k*d0pi vs. d0k*d0pi (topol+PID)
  TH2F *fHistCentrality[3];  //!<!hist. for cent distr (all,sel ev, )
  THnSparseF *fHistMassPtImpPar[5];  //!<! histograms for impact parameter
  THnSparseF *fSparseCutVars[3];     //!<! histograms for cut variation study
  THnSparseF *fHistTrackVar;         //!<! histograms for track cuts study
  THnSparseF *fMCAccPrompt;          //!<!histo for StepMCAcc for Dplus prompt (pt,y,ptB)
  THnSparseF *fMCAccBFeed;           //!<!histo for StepMCAcc for Dplus FD (pt,y,ptB)
  TH2F *fPtVsMassNoPid;      //!<! hist. of pt vs. mass (w/o PID)
  TH2F *fPtVsMass;           //!<! hist. of pt vs. mass (topol+PID cuts)
  TH2F *fPtVsMassPlus;       //!<! hist. of pt vs. mass, D+ candidates (topol+PID cuts)
  TH2F *fPtVsMassMinus;      //!<! hist. of pt vs. mass, D- candidates (topol+PID cuts)
  TH2F *fPtVsMassBadDaus;    //!<! hist. of pt vs. mass (topol+PID cuts)
  TH2F *fPtVsMassGoodDaus;   //!<! hist. of pt vs. mass (topol+PID cuts)
  TH3F *fYVsPtNoPid;         //!<! hist. of Y vs. Pt vs. Mass(w/o PID)
  TH3F *fYVsPt;              //!<! hist. of Y vs. Pt vs. Mass (topol+PID cuts)
  TH2F *fYVsPtSigNoPid;      //!<! hist. of Y vs. Pt (MC, only sig, w/o PID)
  TH2F *fYVsPtSig;           //!<! hist. of Y vs. Pt (MC, only sig, topol+PID cuts)
  TH2F *fPhiEtaCand;         //!<! hist. with eta/phi distribution of candidates
  TH2F *fPhiEtaCandSigReg;   //!<! hist. eta/phi of candidates in D+ mass region
  TH1F *fSPDMult;            //!<! hist. of spd mult
  TH1F* fDaughterClass;      //!<! hist
  TH1F* fDeltaID;            //!<! hist
  TH2F* fIDDauVsIDTra;       //!<! hist

  TH1F** fMassHistLS;        //!<!hist. for LS inv mass (topol+PID)
  TH1F** fCosPHistLS;        //!<!hist. for LS cuts variable 1 (topol+PID)
  TH1F** fDLenHistLS;        //!<!hist. for LS cuts variable 2 (topol+PID)
  TH1F** fSumd02HistLS;      //!<!hist. for LS cuts variable 3 (topol+PID)
  TH1F** fSigVertHistLS;     //!<!hist. for LS cuts variable 4 (topol+PID)
  TH1F** fPtMaxHistLS;       //!<!hist. for LS cuts variable 5 (topol+PID)
  TH1F** fDCAHistLS;         //!<!hist. for LS cuts variable 6 (topol+PID)

  TNtuple *fNtupleDplus; //!<! output ntuple
  Float_t fUpmasslimit;  /// upper inv mass limit for histos
  Float_t fLowmasslimit; /// lower inv mass limit for histos
  Int_t fNPtBins; /// Number of Pt Bins
  Float_t fBinWidth;/// width of one bin in output histos
  TList *fListCuts; /// list of cuts
  AliRDHFCutsDplustoKpipi *fRDCutsAnalysis; /// Cuts for Analysis
  AliNormalizationCounter *fCounter;//!<!Counter for normalization
  Int_t fFillNtuple;   /// flag for filling ntuple 0 no NTuple 1 big Ntuple 2 small NTuple
  Int_t fAODProtection;  /// flag to activate protection against AOD-dAOD mismatch.
                         /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
  Bool_t fReadMC;    /// flag for access to MC
  Bool_t fUseStrangeness;/// flag to enhance strangeness in MC to fit to data
  Bool_t fUseBit;      /// flag to use bitmask
  Bool_t fCutsDistr;    /// flag to activate cuts distr histos
  Bool_t fDoImpPar;    /// flag to activate impact paramter histos
  Bool_t fDoSparse;    /// flag to activate sparses for cut variation study
  Bool_t fDoTrackVarHist; ///flag to activate track variable cut studies
  Bool_t fStepMCAcc;   /// flag to activate histos for StepMCAcc
  Bool_t fUseQuarkTagInKine; /// flag for quark/hadron level identification of prompt and feeddown
  Int_t  fNImpParBins;   /// nunber of bins in impact parameter histos
  Float_t fLowerImpPar;  /// lower limit in impact parameter (um)
  Float_t fHigherImpPar; /// higher limit in impact parameter (um)
  Int_t  fDoLS;        /// flag to do LS analysis
  Int_t fEtaSelection; /// eta region to accept D+ 0=all, -1 = negative, 1 = positive
  Int_t fSystem;   /// 0=pp,1=PbPb
  Int_t fNtrcklMin;   ///minimum number of tracklets
  Int_t fNtrcklMax;   ///maximum number of tracklets
  Bool_t fCutOnTrckl;  ///flag to activate the cut on the number of tracklets
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEDplus,31); /// AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
  /// \endcond
};

#endif
