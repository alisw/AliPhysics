#ifndef ALIANALYSISTASKSEITSSASPECTRA_H
#define ALIANALYSISTASKSEITSSASPECTRA_H

///////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskSE for the extraction of the various histograms to
// study the pt spectra of identified hadrons:
// - log(dEdx)-log(dEdxBB) distributions for pions, kaons and protons in pt bins
// - Pt distributions of pions, kaons and protons with nSigma PID
// Authors:
// E. Biolcati, biolcati@to.infn.it
// L. Milano, milano@to.infn.it
// F. Prino, prino@to.infn.it
// Y. Corrales, corrales@to.infn.it
// N. Jacazio, jacazio@to.infn.it
///////////////////////////////////////////////////////////////////////////

/* $Id$ */
class TF1;
class TH1I;
class TH1F;
class TH2F;
class TList;
class TNtuple;
class TRandom3;

class AliEventCuts;
class AliESDEvent;
class AliESDVertex;
class AliESDtrack;
class AliITSPidParams;
class AliITSPIDResponse;
class AliMCEvent;
class AliVTrack;

#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSEITSsaSpectra : public AliAnalysisTaskSE {

public:
  // 3 track-by-track PID approaches
  // kNSigCut : nSigma with symmetric band and ITSsaHybrid BB parameterization
  // kMeanCut : nSigma with asymmetric band and ITSsaHybrid BB parameterization
  // kLanGaus : Landau-Gauss Bayessian approach
  // in addition, statistic pid analysis available setting fFillIntDistHist
  enum EPID_Type       {kNSigCut, kMeanCut, kLanGaus};
  enum EPileup_Type    {kNoPileup, kPileupSPD, kPileupMV};
  enum EEvtCut_Type {
    kIsReadable = 1,
    kIsNotIncDAQ,
    kPassTrig,
    kIsNotPileup,
		kCorrelation,
    kHasRecVtx,
    kHasGoodVtxZ,
    kIsSDDIn,
		kPassSPDclsVsTCut,
		kPassMultSel,
    kNEvtCuts
  };// event selection criteria
  enum ETrkCut_Type {
    kHasNoSelection = 1,
    kIsITSsa,
    kIsITSrefit,
    kIsNotNeutralParticle,
    kPassSPD,
    kPassPIDcls,
    kPassChi2Ncls,
    kIsInEta,
    kPassdEdx,
    kPassPtCut,
    kPassDCAzcut,
    kPassDCAxycut
  };// track selection criteria
  enum {
    kNchg  =  2, // pos = 0, neg = 1;
    kNspc  =  3, // pion = 0, kaon, 1, proton = 2
    kNbins = 22
  };

  AliAnalysisTaskSEITSsaSpectra();
  virtual ~AliAnalysisTaskSEITSsaSpectra();

  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   LocalInit() {Init();}
  virtual void   UserExec(Option_t*);
  virtual void   Terminate(Option_t*);

  //Setters for event selection settings
  void SetITSPidParams(AliITSPidParams* pidParams)     { fITSPidParams = pidParams; }
  void SetTriggerSel   (UInt_t   tg = AliVEvent::kMB)  { fTriggerSel   = tg;   }
  void SetVtxQACheck   (Bool_t chkSPDres = kTRUE, Bool_t chkZsep = kTRUE, Bool_t reqBoth = kFALSE)
  {fChkVtxSPDRes = chkSPDres; fChkVtxZSep = chkZsep; fReqBothVtx = reqBoth;}
  void SetMaxVtxZCut               (Double_t vz = 10)  { fMaxVtxZCut   = vz;   }
  void SetDoCheckSDDIn          (Bool_t flag = kTRUE)  { fChkIsSDDIn   = flag; }
  void SetDoRejIncDAQ           (Bool_t flag = kTRUE)  { fRejIncDAQ    = flag; }
  void SetDoSPDclsVsTrackletsCut(Bool_t flag = kTRUE)  { fDoSPDCvsTCut = flag; }
  void SetUseExtEventCuts       (Bool_t flag = kTRUE)  { fExtEventCuts = flag; }

  //Setters for track cut settings
  void SetMinSPDPoints   (Int_t     np =   1) { fMinSPDPts       =  np; }
  void SetMinNdEdxSamples(Int_t     np =   3) { fMinNdEdxSamples =  np; }
  void SetAbsEtaCut      (Double_t eta =  .8) { fAbsEtaCut       = eta; }
  void SetMinRapCut      (Double_t   y = -.5) { fMinRapCut       =   y; }
  void SetMaxRapCut      (Double_t   y =  .5) { fMaxRapCut       =   y; }
  void SetCMSRapFct      (Double_t  dy =  .0) { fCMSRapFct       =  dy; }
  void SetMindEdx        (Double_t  md =  .0) { fMindEdx         =  md; }
  void SetMinNSigma      (Double_t  ns = 1.5) { fMinNSigma       =  ns; }
  void SetMaxChi2Clu     (Double_t chi = 2.5) { fMaxChi2Clu      = chi; }

  void SetDCACuts(Double_t nsxy = 7., Double_t nsz = 7.) { fNSigmaDCAxy = nsxy; fNSigmaDCAz  = nsz; }

  //Setters for mult sel.
	void SetMultMethod(unsigned int meth=0) {fMultMethod = meth;}
  void SetMultEstimator(TString centEst = "V0M") {fMultEstimator = centEst;}
  void SetMultiplicityCut(Float_t low, Float_t up)
  {
    if ((up > low) && (!(low < 0.)) && (!(up > 100.))) {
      fLowMult = low; fUpMult = up;
    }
  }

  //Setters for global settings
  void SetPidTech(Int_t tech) {fPidMethod = static_cast<EPID_Type>(tech);}
  void SetUseDefaultPriors(Bool_t flag = kTRUE) { fUseDefaultPriors = flag; }
  void SetIsMC            (Bool_t flag = kTRUE) { fIsMC       = flag; }
  void SetFillNtuple      (Bool_t flag = kTRUE) { fFillNtuple = flag;}
  void SetFillIntDistHist () { fFillIntDistHist = kTRUE; }

  //Setters for pileup settings
  void SetSPDPileupSelection(Int_t cont = 3, Float_t distance = 0.8);
  void SetMVPileUpSelection(Int_t cont = 5, Float_t chi2 = 5, Float_t wgthZdiff = 15., Bool_t chkDiffBC = kFALSE);

	AliEventCuts* GetAliEventCuts() { return &fEventCuts; }

  //Setters for smearing settings
  void SetSmearMC(Double_t smearp, Double_t smeardedx)
  {
    fSmearMC = kTRUE;
    fSmearP = smearp;
    fSmeardEdx = smeardedx;
  }

  // Default event selection cuts
  void   SetupStandardEventCutsForRun1();
  void   SetupEventCutsForRun1pPb();
  void   SetupStandardEventCutsForRun2();

protected:
  Bool_t   IsEventAccepted(EEvtCut_Type& evtSel);
	Bool_t 	 CheckExtraEvtSelStep(EEvtCut_Type& evtSel);
  Bool_t   IsMultSelected();
  Bool_t   IsPileUp();
  Bool_t   IsGoodVtxZ();
  Bool_t   IsGoodSPDvtxRes(const AliESDVertex* spdVtx);
  Bool_t   IsVtxReconstructed();

  Bool_t   IsRapIn(Double_t y) { return (y > fMinRapCut && y < fMaxRapCut) ? kTRUE : kFALSE;}
  Bool_t   DCAcut(Double_t impactXY, Double_t impactZ, Double_t pt) const;
  Bool_t   DCAcutXY(Double_t impactXY, Double_t pt) const;
  Bool_t   DCAcutZ(Double_t impactZ, Double_t pt) const;

  Double_t CookdEdx(Double_t* s) const;
  Double_t Eta2y(Double_t pt, Double_t m, Double_t eta) const;
	
  void     AnalyseMCParticles(AliMCEvent* lMCevent, EEvtCut_Type step, bool lHasGoodVtxGen);
  void     CreateDCAcutFunctions();
  void     PostAllData();

  Int_t    GetTrackPid(AliESDtrack* track, Double_t* logdiff) const;
  Int_t    GetMostProbable(const Double_t* pDens, const Double_t* priors) const;
  void     GetPriors(const AliVTrack* track, Double_t* priors) const;
  void     ComputeBayesProbabilities(Double_t* probs, const Double_t* pDens, const Double_t* prior);

private:
  AliAnalysisTaskSEITSsaSpectra(const AliAnalysisTaskSEITSsaSpectra& source);
  AliAnalysisTaskSEITSsaSpectra& operator=(const AliAnalysisTaskSEITSsaSpectra& source);

  AliESDEvent* fESD; //ESD object
  AliITSPidParams*   fITSPidParams;
  AliITSPIDResponse* fITSPIDResponse; //! class with BB parameterizations

  //Standard event cut
  AliEventCuts fEventCuts;      //!<! basic cut variables for events

  /////////////////////////////////////////////////////
  // List
  /////////////////////////////////////////////////////
  TList* fOutput;      //! list with output
  TList* fListCuts;    //! list with functions storing DCA cut
  TList* fListTree;    //! list with trees
  TList* fListPriors;  //! list with priors in case of bayesian pid approach

  TF1* fDCAxyCutFunc;  // function with DCAz cut vs. pt
  TF1* fDCAzCutFunc;   // function with DCAxy cut vs. pt

  TNtuple* fNtupleData; //! output ntuple
  TNtuple* fNtupleMC;   //! output MC ntuple

  /////////////////////////////////////////////////////
  // General Histos
  /////////////////////////////////////////////////////
  TH1I* fHistNEvents;    //! histo with number of events
  TH1I* fHistMCEvents;
  TH1F* fHistMultBefEvtSel;  //! histo with multiplicity of the events before event selection
  TH1F* fHistMultAftEvtSel;  //! histo with multiplicity of the events after all event selection
  TH1F* fHistVtxZ;       //! histo with the distribution of the primary vertex Z coordinate

  TH2F* fHistNTracks[kNchg];      //! histo with number of tracks vs Pt
  TH2F* fHistDEDX;                //! histo with dedx versus momentum
  TH2F* fHistDEDXdouble;          //! histo with dedx versus signed momentum
  TH2F* fHistNSigmaSep[kNchg * kNspc]; //! histo nsigma separation vs momentum

  // MC histograms with spectra of primaries from the MC truth
  TH2F* fHistPrimMCGenVtxZall[kNchg * kNspc]; //! histo from events with gen Zvtx cut
  TH2F* fHistPrimMCGenVtxZcut[kNchg * kNspc]; //! histo from events with rec Zvtx cut

  //Reconstructed
  TH1F* fHistReco      [kNchg * kNspc]; //! NSigma histos for 6 species
  TH2F* fHistMCReco    [kNchg * kNspc]; //! NSigma histos for 6 species
  TH2F* fHistMCPrimReco[kNchg * kNspc]; //! NSigma histos for 6 species

  // MC histograms using reco values
  TH1F* fHistPrimMCReco[kNchg * kNspc]; //! histo with spectra of primaries from the MC truth
  TH1F* fHistSstrMCReco[kNchg * kNspc]; //! histo with spectra of strange decays from the MC truth
  TH1F* fHistSmatMCReco[kNchg * kNspc]; //! histo with spectra of sec. from material from the MC truth

  // DCAxy distributions
  TH1F* fHistRecoDCA[kNchg * kNspc][kNbins]; //! histo with DCA distibution in the pion hypotesis (positive)

  //DCA Templates
  TH1F* fHistPrimDCA[kNchg * kNspc][kNbins]; //! histo with DCA distibution and dedx PID
  TH1F* fHistSstrDCA[kNchg * kNspc][kNbins]; //! histo with DCA distibution and dedx PID
  TH1F* fHistSmatDCA[kNchg * kNspc][kNbins]; //! histo with DCA distibution and dedx PID

  TH1F* fHistMCtruthPrimDCA[kNchg * kNspc][kNbins]; //! histo with DCA distibution and MC truth PID
  TH1F* fHistMCtruthSstrDCA[kNchg * kNspc][kNbins]; //! histo with DCA distibution and MC truth PID
  TH1F* fHistMCtruthSmatDCA[kNchg * kNspc][kNbins]; //! histo with DCA distibution and MC truth PID

  TH1F* fHistCharge[4]; //! histo with charge distribution to check the calibration

  //dEdx distributions
  TH1F* fHistPosHypPi[kNbins]; //! histo with DCA distibution in the kaons hypotesis (positive)
  TH1F* fHistPosHypKa[kNbins]; //! histo with DCA distibution in the kaons hypotesis (positive)
  TH1F* fHistPosHypPr[kNbins]; //! histo with DCA distibution in the protons hypotesis (positive)
  TH1F* fHistNegHypPi[kNbins]; //! histo with DCA distibution in the pions hypotesis (negative)
  TH1F* fHistNegHypKa[kNbins]; //! histo with DCA distibution in the kaons hypotesis (negative)
  TH1F* fHistNegHypPr[kNbins]; //! histo with DCA distibution in the protons hypotesis (negative)

  //dEdx distributions for MC
  TH1F* fHistMCPosOtherHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosOtherHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosOtherHypProt[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosElHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosElHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosElHypProt[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosPiHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosPiHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosPiHypProt[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosKaHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosKaHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosKaHypProt[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosPrHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosPrHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCPosPrHypProt[kNbins]; //! histo with dedx using the MC truth

  TH1F* fHistMCNegOtherHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegOtherHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegOtherHypProt[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegElHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegElHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegElHypProt[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegPiHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegPiHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegPiHypProt[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegKaHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegKaHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegKaHypProt[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegPrHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegPrHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F* fHistMCNegPrHypProt[kNbins]; //! histo with dedx using the MC truth

  Double_t fPtBinLimits[kNbins + 1]; // limits of Pt Bins

  //evt sel.
  UInt_t   fTriggerSel;
  Double_t fMaxVtxZCut;
  Bool_t   fChkIsSDDIn;    // flag for counrint ev. in p-p 2.76
  Bool_t   fRejIncDAQ;     // flag for reject events with incomplete DAQ
  Bool_t   fDoSPDCvsTCut;  // flag for check compatibility between SPD clusters and tracklets
  Bool_t   fChkVtxSPDRes;  // enable check on spd vtx resolution
  Bool_t   fChkVtxZSep;    // enable check on proximity of the z coordinate between both vertexer
  Bool_t   fReqBothVtx;    // ask for both trk and SPD vertex
  Bool_t   fExtEventCuts;  // enable use of AliEventCuts for event selection

  //trk sel.
  Int_t    fMinSPDPts;       // minimum number of SPD Points
  Int_t    fMinNdEdxSamples; // minimum number of SDD+SSD points
  Double_t fAbsEtaCut;       // limits in pseudorap
  Double_t fMinRapCut;
  Double_t fMaxRapCut;
  Double_t fCMSRapFct;
  Double_t fMindEdx;         // minimum dE/dx value in a layer (to cut noise)
  Double_t fMinNSigma;       // minimum number of sigmas
  Double_t fMaxChi2Clu;      // maximum cluster
  Double_t fNSigmaDCAxy;     // DCA cut in bend. plane
  Double_t fNSigmaDCAz;      // DCA cut along z

  //mult sel.
  unsigned int fMultMethod;    // method for cent/mult values: 0=skip mult sel, 1=new cent framework, 2=old cent framework, 3=tracks+tracklets, 4=tracklets, 5=cluster on SPD 	
  TString      fMultEstimator; // centrality/multiplicity framework estimator name
  Float_t      fLowMult;       // low Centrality cut
  Float_t      fUpMult;        // up  Centrality cut
  Float_t      fEvtMult;       // event multiplicity -0.5 by default

  //Global setting
  Int_t  fYear;           // FIXME Year (2009, 2010)
  EPID_Type fPidMethod;   // track-by-track pid approach

  Bool_t fUseDefaultPriors; // flag to use default(equal) priors
  Bool_t fIsMC;             // flag to switch on the MC analysis for the efficiency estimation
  Bool_t fFillNtuple;       // flag to fill ntuples
  Bool_t fFillIntDistHist;  // flag to fill histogram with information for statistic pid analysis

  //Pileup selection setting
  EPileup_Type  fPlpType;

  Int_t    fMinPlpContribMV; //minimum contributors to the pilup vertices, multi-vertex
  Float_t  fMaxPlpChi2MV;    //minimum value of Chi2perNDF of the pileup vertex, multi-vertex
  Float_t  fMinWDistMV;      //minimum of the sqrt of weighted distance between the primary and the pilup vertex, multi-vertex
  Bool_t   fCheckPlpFromDifferentBCMV; //pileup from different BC (Bunch Crossings)

  Int_t    fMinPlpContribSPD; //minimum contributors to the pilup vertices, SPD
  Float_t  fMinPlpZdistSPD;   //minimum distance for the SPD pileup vertex

  //smearing
  TRandom3* fRandGener; // generator for smearing
  Bool_t    fSmearMC;   // flag to apply extra smearing on MC
  Double_t  fSmearP;    // extra relative smearing on simulated momentum
  Double_t  fSmeardEdx; // extra relative smearing on simulated dE/dx

  ClassDef(AliAnalysisTaskSEITSsaSpectra, 11);
};

#endif
