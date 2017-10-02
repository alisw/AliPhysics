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
class TH3F;
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
  enum EPileup_Type    {
		kNoPileup,
		kPileupSPD,
		kPileupInMultBins,
		kPileupMV
	};
  enum EEvtCut_Type {
    kIsReadable=1,
		kPassMultSel,
		kIsSDDIn,
    kIsNotIncDAQ,
    kPassTrig,
		kPassINELgtZERO,
		kCorrelations,
		kPassSPDclsVsTCut,
    kIsPileupSPD,
		kIsPileupSPDinMultBins,
		kIsPileupMV,
    kHasRecVtx,
    kHasGoodVtxZ,
    kNEvtCuts
  };///< event selection criteria
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
  };///< track selection criteria
  enum {
    kNchg  =  2, ///< pos = 0, neg = 1;
    kNspc  =  3, ///< pion = 0, kaon, 1, proton = 2
    kNbins = 22  ///< deprecated parameter
  };

  AliAnalysisTaskSEITSsaSpectra();
  virtual ~AliAnalysisTaskSEITSsaSpectra();

  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   LocalInit() {Init();}
  virtual void   UserExec(Option_t*);
  virtual void   Terminate(Option_t*);

  //Setters for histo bins
  void SetBins     (int nbins, float min, float max, float* bins);
  void SetCentBins (int nbins, float *bins);
  void SetDCABins  (int nbins, float min, float max);
  void SetDCABins  (int nbins, float* bins);
  void SetPtBins   (int nbins, float *bins);

  //Setters for event selection settings
  void SetTriggerSel   (UInt_t   tg = AliVEvent::kMB)  { fTriggerSel   = tg;   }
  void SetMaxVtxZCut               (Double_t vz = 10)  { fMaxVtxZCut   = vz;   }
	void SetChkIsEventINELgtZERO (Bool_t chk=kTRUE)      { fChkIsEventINELgtZERO = chk; }
  void SetDoCheckSDDIn          (Bool_t flag = kTRUE)  { fChkIsSDDIn   = flag; }
  void SetDoRejIncDAQ           (Bool_t flag = kTRUE)  { fRejIncDAQ    = flag; }
  void SetDoSPDclsVsTrackletsCut(Bool_t flag = kTRUE)  { fDoSPDCvsTCut = flag; }
	void SetUseSelectVertex2015pp (Bool_t flag = kTRUE)  { fUseSelectVertex2015pp = flag; }
	void SetVtxQACheck   (Bool_t chkSPDres = kTRUE, Bool_t chkZsep = kTRUE, Bool_t reqBoth = kFALSE)
  {fChkVtxSPDRes = chkSPDres; fChkVtxZSep = chkZsep; fReqBothVtx = reqBoth;}
  void SetUseExtEventCuts       (Bool_t flag = kTRUE)  { fExtEventCuts = flag; }

  //Setters for mult sel.
	void SetMultMethod(unsigned int meth=0) {fMultMethod = meth;}
	void SetMultEvSel(Bool_t flag=kFALSE){ fMultEvSel = flag;}
  void SetMultEstimator(TString centEst = "V0M") {fMultEstimator = centEst;}
  void SetMultiplicityCut(Float_t low, Float_t up)
  {
    if ((up > low) && (!(low < 0.)) && (!(up > 100.))) {
      fLowMult = low; fUpMult = up;
    }
  }

  //Setters for pileup settings
	void SetPileupSelection(unsigned long plp=1ULL) { fPlpType = plp; }
  void SetSPDPileupSelection(Int_t cont = 3, Float_t distance = 0.8)
	{ fMinPlpContribSPD = cont; fMinPlpZdistSPD = distance; }
  void SetMVPileUpSelection(Int_t cont = 5, Float_t chi2 = 5, Float_t wgthZdiff = 15., Bool_t chkDiffBC = kFALSE)
	{ fMinPlpContribMV = cont; fMaxPlpChi2MV = chi2; fMinWDistMV = wgthZdiff; fCheckPlpFromDifferentBCMV = chkDiffBC; }

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

  //Setters for global settings
  void SetPidTech(Int_t tech) {fPidMethod = static_cast<EPID_Type>(tech);}
  void SetUseDefaultPriors(Bool_t flag = kTRUE) { fUseDefaultPriors = flag; }
	void SetITSPidParams(AliITSPidParams* pidParams)     { fITSPidParams = pidParams; }
  void SetIsMC            (Bool_t flag = kTRUE) { fIsMC       = flag; }
  void SetFillNtuple      (Bool_t flag = kTRUE) { fFillNtuple = flag;}
  void SetFillIntDistHist () { fFillIntDistHist = kTRUE; }

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
  Bool_t   IsMultSelected();
  UInt_t   IsPileup();
  Bool_t   IsGoodVtxZ();
  Bool_t   IsGoodSPDvtxRes(const AliESDVertex* spdVtx);
  Bool_t   SelectVertex2015pp();

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

  AliESDEvent*       fESD;            ///<  ESD object
  AliITSPidParams*   fITSPidParams;   //!<! parameterization used for PID
  AliITSPIDResponse* fITSPIDResponse; //!<! class with BB parameterizations

  //Standard event cut
  AliEventCuts fEventCuts;            //!<! basic cut variables for events

  /////////////////////////////////////////////////////
  // List
  /////////////////////////////////////////////////////
  TList*   fOutput;                    //!<! list with output
  TList*   fListCuts;                  //!<! list with functions storing DCA cut
  TList*   fListTree;                  //!<! list with trees
  TList*   fListPriors;                //!<! list with priors in case of bayesian pid approach

  TF1*     fDCAxyCutFunc;              ///< function with DCAz cut vs. pt
  TF1*     fDCAzCutFunc;               ///< function with DCAxy cut vs. pt

  TNtuple* fNtupleData;                //!<! output ntuple
  TNtuple* fNtupleMC;                  //!<! output MC ntuple

  /////////////////////////////////////////////////////
  // General Histos
  /////////////////////////////////////////////////////
  TH2I* fHistNEvents;                           //!<! histo with number of events / mult bin
  TH2I* fHistMCEvents;                          //!<! histo with MC number of events / mult bin
  TH1F* fHistMultBefEvtSel;                     //!<! histo with multiplicity of the events before event selection
  TH1F* fHistMultAftEvtSel;                     //!<! histo with multiplicity of the events after all event selection
  TH2F* fHistVtxZ;                              //!<! histo with the distribution of the primary vertex Z coordinate

  TH3F* fHistNTracks[kNchg];                    //!<! histo with number of tracks vs Pt
  TH2F* fHistDEDX;                              //!<! histo with dedx versus momentum
  TH2F* fHistDEDXdouble;                        //!<! histo with dedx versus signed momentum
  TH2F* fHistNSigmaSep[kNchg * kNspc];          //!<! histo nsigma separation vs momentum

  // MC histograms with spectra of primaries from the MC truth
  TH3F* fHistPrimMCGenVtxZall[kNchg * kNspc];   //!<! histo from events with gen Zvtx cut
  TH3F* fHistPrimMCGenVtxZcut[kNchg * kNspc];   //!<! histo from events with rec Zvtx cut

  //Reconstructed
  TH2F* fHistReco      [kNchg * kNspc];         //!<! NSigma histos for 6 species
  TH3F* fHistMCReco    [kNchg * kNspc];         //!<! NSigma histos for 6 species
  TH3F* fHistMCPrimReco[kNchg * kNspc];         //!<! NSigma histos for 6 species

  // MC histograms using reco values
  TH2F* fHistPrimMCReco[kNchg * kNspc];         //!<! histo with spectra of primaries from the MC truth
  TH2F* fHistSstrMCReco[kNchg * kNspc];         //!<! histo with spectra of strange decays from the MC truth
  TH2F* fHistSmatMCReco[kNchg * kNspc];         //!<! histo with spectra of sec. from material from the MC truth

  // DCAxy distributions
  TH3F* fHistRecoDCA[kNchg * kNspc];            //!<! histo with DCA distibution in the pion hypotesis (positive)

  //DCA Templates
  TH3F* fHistPrimDCA[kNchg * kNspc];            //!<! histo with DCA distibution and dedx PID
  TH3F* fHistSstrDCA[kNchg * kNspc];            //!<! histo with DCA distibution and dedx PID
  TH3F* fHistSmatDCA[kNchg * kNspc];            //!<! histo with DCA distibution and dedx PID

  TH3F* fHistMCtruthPrimDCA[kNchg * kNspc];     //!<! histo with DCA distibution and MC truth PID
  TH3F* fHistMCtruthSstrDCA[kNchg * kNspc];     //!<! histo with DCA distibution and MC truth PID
  TH3F* fHistMCtruthSmatDCA[kNchg * kNspc];     //!<! histo with DCA distibution and MC truth PID

  TH1F* fHistCharge[4];                         //!<! histo with charge distribution to check the calibration

  //dEdx distributions
  TH2F* fHistPosHypPi; //! histo with DCA distibution in the kaons hypotesis (positive)
  TH2F* fHistPosHypKa; //! histo with DCA distibution in the kaons hypotesis (positive)
  TH2F* fHistPosHypPr; //! histo with DCA distibution in the protons hypotesis (positive)
  TH2F* fHistNegHypPi; //! histo with DCA distibution in the pions hypotesis (negative)
  TH2F* fHistNegHypKa; //! histo with DCA distibution in the kaons hypotesis (negative)
  TH2F* fHistNegHypPr; //! histo with DCA distibution in the protons hypotesis (negative)

  //dEdx distributions for MC
  TH2F* fHistMCPosOtherHypPion; //! histo with dedx using the MC truth
  TH2F* fHistMCPosOtherHypKaon; //! histo with dedx using the MC truth
  TH2F* fHistMCPosOtherHypProt; //! histo with dedx using the MC truth
  TH2F* fHistMCPosElHypPion; //! histo with dedx using the MC truth
  TH2F* fHistMCPosElHypKaon; //! histo with dedx using the MC truth
  TH2F* fHistMCPosElHypProt; //! histo with dedx using the MC truth
  TH2F* fHistMCPosPiHypPion; //! histo with dedx using the MC truth
  TH2F* fHistMCPosPiHypKaon; //! histo with dedx using the MC truth
  TH2F* fHistMCPosPiHypProt; //! histo with dedx using the MC truth
  TH2F* fHistMCPosKaHypPion; //! histo with dedx using the MC truth
  TH2F* fHistMCPosKaHypKaon; //! histo with dedx using the MC truth
  TH2F* fHistMCPosKaHypProt; //! histo with dedx using the MC truth
  TH2F* fHistMCPosPrHypPion; //! histo with dedx using the MC truth
  TH2F* fHistMCPosPrHypKaon; //! histo with dedx using the MC truth
  TH2F* fHistMCPosPrHypProt; //! histo with dedx using the MC truth

  TH2F* fHistMCNegOtherHypPion; //! histo with dedx using the MC truth
  TH2F* fHistMCNegOtherHypKaon; //! histo with dedx using the MC truth
  TH2F* fHistMCNegOtherHypProt; //! histo with dedx using the MC truth
  TH2F* fHistMCNegElHypPion; //! histo with dedx using the MC truth
  TH2F* fHistMCNegElHypKaon; //! histo with dedx using the MC truth
  TH2F* fHistMCNegElHypProt; //! histo with dedx using the MC truth
  TH2F* fHistMCNegPiHypPion; //! histo with dedx using the MC truth
  TH2F* fHistMCNegPiHypKaon; //! histo with dedx using the MC truth
  TH2F* fHistMCNegPiHypProt; //! histo with dedx using the MC truth
  TH2F* fHistMCNegKaHypPion; //! histo with dedx using the MC truth
  TH2F* fHistMCNegKaHypKaon; //! histo with dedx using the MC truth
  TH2F* fHistMCNegKaHypProt; //! histo with dedx using the MC truth
  TH2F* fHistMCNegPrHypPion; //! histo with dedx using the MC truth
  TH2F* fHistMCNegPrHypKaon; //! histo with dedx using the MC truth
  TH2F* fHistMCNegPrHypProt; //! histo with dedx using the MC truth

  TArrayF fCentBins;
  TArrayF fDCABins;
  TArrayF fPtBins;

  //evt sel.
  UInt_t   fTriggerSel;
  Double_t fMaxVtxZCut;
	Bool_t   fChkIsEventINELgtZERO;   // flag for select INEL>0 events
  Bool_t   fChkIsSDDIn;             // flag for counrint ev. in p-p 2.76
  Bool_t   fRejIncDAQ;              // flag for reject events with incomplete DAQ
  Bool_t   fDoSPDCvsTCut;           // flag for check compatibility between SPD clusters and tracklets
	Bool_t   fUseSelectVertex2015pp;  // flag to select vertex based on criteria for 2015 pp data
  Bool_t   fChkVtxSPDRes;           // enable check on spd vtx resolution
  Bool_t   fChkVtxZSep;             // enable check on proximity of the z coordinate between both vertexer
  Bool_t   fReqBothVtx;             // ask for both trk and SPD vertex
  Bool_t   fExtEventCuts;           // enable use of AliEventCuts for event selection
  //mult sel.
  unsigned int fMultMethod;    // method for cent/mult values: 0=skip mult sel, 1=new cent framework, 2=old cent framework, 3=tracks+tracklets, 4=tracklets, 5=cluster on SPD 	
  TString      fMultEstimator; // centrality/multiplicity framework estimator name
	Bool_t       fMultEvSel;
  Float_t      fLowMult;       // low Centrality cut
  Float_t      fUpMult;        // up  Centrality cut
  Float_t      fEvtMult;       // event multiplicity -0.5 by default
  //Pileup selection setting
  unsigned long  fPlpType;
	//PileupSPD settings
  Int_t    fMinPlpContribSPD; //minimum contributors to the pilup vertices, SPD
  Float_t  fMinPlpZdistSPD;   //minimum distance for the SPD pileup vertex
  //PileupMV settings
  Int_t    fMinPlpContribMV; //minimum contributors to the pilup vertices, multi-vertex
  Float_t  fMaxPlpChi2MV;    //minimum value of Chi2perNDF of the pileup vertex, multi-vertex
  Float_t  fMinWDistMV;      //minimum of the sqrt of weighted distance between the primary and the pilup vertex, multi-vertex
  Bool_t   fCheckPlpFromDifferentBCMV; //pileup from different BC (Bunch Crossings)

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

  //Global setting
  Int_t  fYear;           // FIXME Year (2009, 2010)
  EPID_Type fPidMethod;   // track-by-track pid approach

  Bool_t fUseDefaultPriors; // flag to use default(equal) priors
  Bool_t fIsMC;             // flag to switch on the MC analysis for the efficiency estimation
  Bool_t fFillNtuple;       // flag to fill ntuples
  Bool_t fFillIntDistHist;  // flag to fill histogram with information for statistic pid analysis

  //smearing
  TRandom3* fRandGener; // generator for smearing
  Bool_t    fSmearMC;   // flag to apply extra smearing on MC
  Double_t  fSmearP;    // extra relative smearing on simulated momentum
  Double_t  fSmeardEdx; // extra relative smearing on simulated dE/dx

  ClassDef(AliAnalysisTaskSEITSsaSpectra, 12);
};

#endif
