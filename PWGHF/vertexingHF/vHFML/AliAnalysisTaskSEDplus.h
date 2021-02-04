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
/// \author F. Grosa, fabrizio.grosa@cern.ch
/// \author F. Catalano, fabio.catalano@cern.ch
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TArrayD.h>

#include "AliRDHFCutsDplustoKpipi.h"
#include "AliHFMLResponseDplustoKpipi.h"
#include "AliHFMLVarHandlerDplustoKpipi.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliNormalizationCounter.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

class AliAnalysisTaskSEDplus : public AliAnalysisTaskSE
{
 public:

  enum {kpp, kPbPb};

  AliAnalysisTaskSEDplus();
  AliAnalysisTaskSEDplus(const char *name, AliRDHFCutsDplustoKpipi* analysiscuts, Bool_t createMLtree);
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
  void SetAODMismatchProtection(Int_t opt=0) {fAODProtection=opt;}

  void SetUseOnlyPositiveEta(){fEtaSelection=1;}
  void SetUseOnlyNegativeEta(){fEtaSelection=-1;}
  void SetUseFullEta(){fEtaSelection=0;}
  void SetUseQuarkLevelTag(Bool_t opt){fUseQuarkTagInKine=opt;}

  void SetCutOnNtracklets(Bool_t applycut=kTRUE, Int_t Ntrckmin=0, Int_t Ntrckmax=100) {
    fCutOnTrckl=applycut;
    fNtrcklMin=Ntrckmin;
    fNtrcklMax=Ntrckmax;
  }

  void SetFillOnlySignalSparses(Bool_t fillonlysig=kTRUE) {fFillOnlySignalSparses=fillonlysig;}
  void SetUseFinePtBinsForSparse(Bool_t usefinebins=kTRUE) {fUseFinPtBinsForSparse=usefinebins;} //use only in case of few candidates (e.g. MC signal only)
  void SetKeepOnlyBkgFromHIJING(Bool_t keeponlyhijing=kTRUE) {fKeepOnlyBkgFromHIJING = keeponlyhijing;}

  /// methods for ML application
  void SetDoMLApplication(Bool_t flag = kTRUE, Bool_t isMultiClass = kFALSE) {fApplyML = flag; fMultiClassML = isMultiClass;}
  void SetMLConfigFile(TString path = "") {fConfigPath = path;}
  void SetMLBinsForSparse(Int_t nbins = 300, Double_t min = 0.85, Double_t max = 1.){ fNMLBins[0] = nbins; fMLOutputMin[0] = min; fMLOutputMax[0] = max;}
  void SetMultiClassMLBinsForSparse(Int_t nbinsBkg = 100,
                                    Int_t nbinsPrompt = 100,
                                    Int_t nbinsFD = 100,
                                    Double_t minBkg = 0., Double_t maxBkg = 1.,
                                    Double_t minPrompt = 0., Double_t maxPrompt = 1.,
                                    Double_t minFD = 0., Double_t maxFD = 1.) 
  {
    fNMLBins[0] = nbinsBkg; fNMLBins[1] = nbinsPrompt; fNMLBins[2] = nbinsFD;
    fMLOutputMin[0] = minBkg; fMLOutputMin[1] = minPrompt; fMLOutputMin[2] = minFD;
    fMLOutputMax[0] = maxBkg; fMLOutputMax[1] = maxPrompt; fMLOutputMax[2] = maxFD;
  }
  void SetMinimalVarForMLSparse(Bool_t flag = kTRUE) {fUseMinimalVarForSparse = flag;}
  /// methods for ML tree creation
  void SetCreateMLTree(Bool_t flag = kTRUE) {fCreateMLtree = flag;}
  void SetMLTreePIDopt(int opt) {fPIDopt = opt;} // default AliHFMLVarHandlerDplustoKpipi::kNsigmaDetAndCombPID
  void SetMLTreeAddTrackVar(Bool_t flag = kTRUE) {fAddSingleTrackVar = flag;}
  void SetMLTreeAddImpParProd(Bool_t flag = kTRUE) {fAddImpParProdProngs = flag;}
  void SetFillOnlySignalInMLtree(Bool_t opt = kTRUE) {
    if(fReadMC) fFillOnlySignal = opt;
    else {
      if(opt)
        AliError("fReadMC has to be kTRUE");
    }
  }

  void EnableMLTreeEvtSampling(Float_t fractokeep, ULong_t seed) {
    fEnableEvtSampling = kTRUE;
    fFracEvtToKeep = fractokeep;
    fSeedSampling = seed;
  }
  void EnableMLTreeCandSampling(Float_t fractokeep, Float_t maxptsampling) {
    fEnableCandSampling = kTRUE;
    fFracCandToKeep = fractokeep;
    fMaxCandPtSampling = maxptsampling;
  }

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
  void FillMCAcceptanceHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader);

  enum
  {
    kVarForSparse = 13,
    kVarForSparseFD = 14,
    knVarForSparseMLMinimal = 3,
    kVarForSparseAcc = 2,
    kVarForSparseAccFD = 3,
    kVarForTrackSparse = 7,
    kVarForImpPar = 3
  };

  TList* fOutput = nullptr;                                           //!<! list send on output slot 0
  TH1F* fHistNEvents = nullptr;                                       //!<! hist. for No. of events
  TH1F* fHistNCandidates = nullptr;                                   //!<! hist. for No. of candidates
  TH1F** fMassHist = nullptr;                                         //!<! hist. for inv mass (topol+PID cuts)
  TH1F** fMassHistPlus = nullptr;                                     //!<! hist. for D+ inv mass (topol+PID cuts)
  TH1F** fMassHistMinus = nullptr;                                    //!<! hist. for D- inv mass (topol+PID cuts)
  TH1F** fMassHistNoPid = nullptr;                                    //!<! hist. for inv mass (w/o PID)
  TH1F** fCosPHist = nullptr;                                         //!<! hist. for PointingAngle (topol+PID)
  TH1F** fDLenHist = nullptr;                                         //!<! hist. for Dec Length (topol+PID)
  TH1F** fSumd02Hist = nullptr;                                       //!<! hist. for sum d02 (topol+PID)
  TH1F** fSigVertHist = nullptr;                                      //!<! hist. for sigVert (topol+PID)
  TH1F** fPtMaxHist = nullptr;                                        //!<! hist. for Pt Max (topol+PID)
  TH1F** fPtKHist = nullptr;                                          //!<! hist. for PtK (topol+PID)
  TH1F** fPtpi1Hist = nullptr;                                        //!<! hist. for PtPi1 (topol+PID)
  TH1F** fPtpi2Hist = nullptr;                                        //!<! hist. for PtPi2 (topol+PID)
  TH1F** fDCAHist = nullptr;                                          //!<! hist. for DCA (topol+PID)
  TH1F** fDLxy = nullptr;                                             //!<! hist. for DLxy (topol+PID)
  TH1F** fCosxy = nullptr;                                            //!<! hist. for Cosxy (topol+PID)
  TH2F *fCorreld0Kd0pi[3] = {};                                       //!<! hist. for d0k*d0pi vs. d0k*d0pi (topol+PID)
  TH2F *fHistCentrality[3] = {};                                      //!<! hist. for cent distr (all,sel ev, )
  THnSparseF *fHistMassPtImpPar[5] = {};                              //!<! histograms for impact parameter
  THnSparseF *fSparseCutVars[3] = {};                                 //!<! histograms for cut variation study
  THnSparseF *fHistTrackVar = nullptr;                                //!<! histograms for track cuts study
  THnSparseF *fMCAccPrompt = nullptr;                                 //!<! histo for StepMCAcc for Dplus prompt (pt,y,ptB)
  THnSparseF *fMCAccBFeed = nullptr;                                  //!<! histo for StepMCAcc for Dplus FD (pt,y,ptB)
  TH2F *fPtVsMassNoPid = nullptr;                                     //!<! hist. of pt vs. mass (w/o PID)
  TH2F *fPtVsMass = nullptr;                                          //!<! hist. of pt vs. mass (topol+PID cuts)
  TH2F *fPtVsMassPlus = nullptr;                                      //!<! hist. of pt vs. mass, D+ candidates (topol+PID cuts)
  TH2F *fPtVsMassMinus = nullptr;                                     //!<! hist. of pt vs. mass, D- candidates (topol+PID cuts)
  TH2F *fPtVsMassBadDaus = nullptr;                                   //!<! hist. of pt vs. mass (topol+PID cuts)
  TH2F *fPtVsMassGoodDaus = nullptr;                                  //!<! hist. of pt vs. mass (topol+PID cuts)
  TH3F *fYVsPtNoPid = nullptr;                                        //!<! hist. of Y vs. Pt vs. Mass(w/o PID)
  TH3F *fYVsPt = nullptr;                                             //!<! hist. of Y vs. Pt vs. Mass (topol+PID cuts)
  TH2F *fYVsPtSigNoPid = nullptr;                                     //!<! hist. of Y vs. Pt (MC, only sig, w/o PID)
  TH2F *fYVsPtSig = nullptr;                                          //!<! hist. of Y vs. Pt (MC, only sig, topol+PID cuts)
  TH2F *fPhiEtaCand = nullptr;                                        //!<! hist. with eta/phi distribution of candidates
  TH2F *fPhiEtaCandSigReg = nullptr;                                  //!<! hist. eta/phi of candidates in D+ mass region
  TH1F *fSPDMult = nullptr;                                           //!<! hist. of spd mult
  TH1F* fDaughterClass = nullptr;                                     //!<! hist
  TH1F* fDeltaID = nullptr;                                           //!<! hist
  TH2F* fIDDauVsIDTra = nullptr;                                      //!<! hist

  TH1F** fMassHistLS = nullptr;                                       //!<!hist. for LS inv mass (topol+PID)
  TH1F** fCosPHistLS = nullptr;                                       //!<!hist. for LS cuts variable 1 (topol+PID)
  TH1F** fDLenHistLS = nullptr;                                       //!<!hist. for LS cuts variable 2 (topol+PID)
  TH1F** fSumd02HistLS = nullptr;                                     //!<!hist. for LS cuts variable 3 (topol+PID)
  TH1F** fSigVertHistLS = nullptr;                                    //!<!hist. for LS cuts variable 4 (topol+PID)
  TH1F** fPtMaxHistLS = nullptr;                                      //!<!hist. for LS cuts variable 5 (topol+PID)
  TH1F** fDCAHistLS = nullptr;                                        //!<!hist. for LS cuts variable 6 (topol+PID)

  TH2F* fHistChi2OvernClsVsPtD[3] = {};                               //!<! hist for chi2/ndf vs ptD before and after D-meson selection

  Float_t fUpmasslimit = 1.965;                                       /// upper inv mass limit for histos
  Float_t fLowmasslimit = 1.765;                                      /// lower inv mass limit for histos
  Int_t fNPtBins = 1;                                                 /// Number of Pt Bins
  Float_t fBinWidth = 0.002;                                          /// width of one bin in output histos
  TList *fListCuts = nullptr;                                         /// list of cuts
  AliRDHFCutsDplustoKpipi *fRDCutsAnalysis = nullptr;                 /// Cuts for Analysis
  AliNormalizationCounter *fCounter = nullptr;                        //!<! Counter for normalization
  Int_t fAODProtection = 0;                                           /// flag to activate protection against AOD-dAOD mismatch.
                                                                      /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
  Bool_t fReadMC = kFALSE;                                            /// flag for access to MC
  Bool_t fUseStrangeness = kFALSE;                                    /// flag to enhance strangeness in MC to fit to data
  Bool_t fUseBit = kTRUE;                                             /// flag to use bitmask
  Bool_t fCutsDistr = kFALSE;                                         /// flag to activate cuts distr histos
  Bool_t fDoImpPar = kFALSE;                                          /// flag to activate impact paramter histos
  Bool_t fDoSparse = kFALSE;                                          /// flag to activate sparses for cut variation study
  Bool_t fDoTrackVarHist = kFALSE;                                    /// flag to activate track variable cut studies
  Bool_t fStepMCAcc = kFALSE;                                         /// flag to activate histos for StepMCAcc
  Bool_t fUseQuarkTagInKine = kTRUE;                                  /// flag for quark/hadron level identification of prompt and feeddown
  Int_t  fNImpParBins = 400;                                          /// nunber of bins in impact parameter histos
  Float_t fLowerImpPar = -1000.;                                      /// lower limit in impact parameter (um)
  Float_t fHigherImpPar = 1000.;                                      /// higher limit in impact parameter (um)
  Int_t  fDoLS = 0;                                                   /// flag to do LS analysis
  Int_t fEtaSelection = 0;                                            /// eta region to accept D+ 0=all, -1 = negative, 1 = positive
  Int_t fSystem = kpp;                                                /// 0=pp,1=PbPb
  Int_t fNtrcklMin = 0;                                               /// minimum number of tracklets
  Int_t fNtrcklMax = 1000000000;                                      /// maximum number of tracklets
  Bool_t fCutOnTrckl = kFALSE;                                        /// flag to activate the cut on the number of tracklets
  Bool_t fFillOnlySignalSparses = kFALSE;                             /// flag to fill only signal sparses in case of MC
  Bool_t fUseFinPtBinsForSparse = kFALSE;                             /// flag to fill pt axis of sparse with 0.1 GeV/c wide bins
  Bool_t fKeepOnlyBkgFromHIJING = kFALSE;                             /// flag to keep background from HIJING only

  /// variables for ML application
  Bool_t fApplyML = kFALSE;                                           /// flag to enable ML application
  Bool_t fMultiClassML = kFALSE;                                      /// flag to enable multi-class models (bkg, prompt, FD)
  TString fConfigPath = "";                                           /// path to ML config file
  AliHFMLResponseDplustoKpipi* fMLResponse = nullptr;                 //!<! object to handle ML response
  Int_t fNMLBins[3] = {100, 100, 100};                                /// number of bins for ML output axis in THnSparse
  Double_t fMLOutputMin[3] = {0.0, 0.0, 0.0};                         /// min for ML output axis in THnSparse
  Double_t fMLOutputMax[3] = {1.0, 1.0, 1.0};                         /// max for ML output axis in THnSparse

  /// variables for tree creation
  Bool_t fCreateMLtree = kFALSE;
  AliHFMLVarHandlerDplustoKpipi* fMLhandler = nullptr;                //!<! object to handle ML tree creation and filling
  TTree* fMLtree = nullptr;                                           //!<! tree with candidates for ML
  int fPIDopt = AliHFMLVarHandlerDplustoKpipi::kNsigmaDetAndCombPID;  /// option for PID variables
  Bool_t fAddSingleTrackVar = kFALSE;                                 /// option to store single track variables
  Bool_t fAddImpParProdProngs = kFALSE;                               /// option to store d0K*d0pi1 and d0K*d0pi2 variables
  Bool_t fFillOnlySignal = kFALSE;                                    /// option to store only signal when using MC
  Bool_t fUseMinimalVarForSparse = kFALSE;                            /// flag to keep only inv. mass, pt and prob. in the sparse

  Bool_t fEnableEvtSampling = kFALSE;                                 /// flag to apply event sampling
  Float_t fFracEvtToKeep = 1.1;                                       /// fraction of events to be kept by event sampling
  ULong_t fSeedSampling = 0;                                          /// seed for sampling
  Bool_t fEnableCandSampling = kFALSE;                                /// flag to apply candidate sampling
  Float_t fFracCandToKeep = 1.1;                                      /// fraction of candidates to be kept by sampling
  Float_t fMaxCandPtSampling = 0.;                                    /// maximun candidate pt to apply sampling

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEDplus,38); /// AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
  /// \endcond
};

#endif
