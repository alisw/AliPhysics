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

#include <THnSparse.h>
#include <TFile.h>
#include <TAxis.h>
#include <math.h>

#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSEITSsaSpectra : public AliAnalysisTaskSE
{

 public:
  // 3 track-by-track PID approaches
  // kNSigCut : nSigma with symmetric band and ITSsaHybrid BB parameterization
  // kMeanCut : nSigma with asymmetric band and ITSsaHybrid BB parameterization
  // kLanGaus : Landau-Gauss Bayessian approach
  // in addition, statistic pid analysis available by setting fFillIntDistHist
  enum EPID_Type { kNSigCut, kMeanCut, kLanGaus };
  enum EPileup_Type { kNoPileup, kPileupSPD, kPileupInMultBins, kPileupMV };
  enum EEvtCut_Type {
    kIsReadable = 1,
    kPassMultSel,
    kIsSDDIn,
    kIsNotIncDAQ,
    kPassTrig,
    kPassINELgtZERO,
    kCorrelations,
    kPassSPDclsVsTCut,
    kIsNotPileupSPD,
    kIsNotPileupSPDinMultBins,
    kIsNotPileupMV,
    kIsNotPileup,
    kHasRecVtx,
    kHasGoodVtxZ,
    kNEvtCuts
  }; ///< event selection criteria
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
  }; ///< track selection criteria
  enum {
    kNchg = 2,  ///< pos = 0, neg = 1;
    kNspc = 3,  ///< pion = 0, kaon, 1, proton = 2
    kNbins = 22 ///< deprecated parameter
  };

  AliAnalysisTaskSEITSsaSpectra(bool __def_prior = true, bool __fill_ntuple = false);
  virtual ~AliAnalysisTaskSEITSsaSpectra();

  virtual void UserCreateOutputObjects();
  virtual void Initialization();
  virtual void LocalInit() { Initialization(); }
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  // Setters for histo bins
  void SetBins(const int nbins, double min, double max, double *bins);
  void SetCentBins(int nbins, double *bins);
  void SetDCABins(int nbins, double *bins);
  void SetPtBins(int nbins, double *bins);

  //Setter for unfolded probability matrices
  void SetUnfoldingProb(const char *filepath);

  // Setters for event selection settings
  void SetTriggerSel(UInt_t tg = AliVEvent::kMB) { fTriggerSel = tg; }
  void SetMaxVtxZCut(double vz = 10) { fMaxVtxZCut = vz; }
  void SetChkIsEventINELgtZERO(bool chk = kTRUE) { fChkIsEventINELgtZERO = chk; }
  void SetDoCheckSDDIn(bool flag = kTRUE) { fChkIsSDDIn = flag; }
  void SetDoRejIncDAQ(bool flag = kTRUE) { fRejIncDAQ = flag; }
  void SetDoSPDclsVsTrackletsCut(bool flag = kTRUE) { fDoSPDCvsTCut = flag; }
  void SetUseSelectVertex2015pp(bool flag = kTRUE) { fUseSelectVertex2015pp = flag; }
  void SetVtxQACheck(bool chkSPDres = kTRUE, bool chkZsep = kTRUE, bool reqBoth = kFALSE)
  {
    fChkVtxSPDRes = chkSPDres;
    fChkVtxZSep = chkZsep;
    fReqBothVtx = reqBoth;
  }
  void SetUseExtEventCuts(bool flag = kTRUE) { fExtEventCuts = flag; }

  // Setters for mult sel.
  void SetMultMethod(unsigned int meth = 0) { fMultMethod = meth; }
  void SetMultEvSel(bool flag = kFALSE) { fMultEvSel = flag; }
  void SetMultEstimator(TString centEst = "V0M") { fMultEstimator = centEst; }
  void SetMultiplicityCut(float low, float up)
  {
    if ((up > low) && (!(low < 0.)) && (!(up > 100.))) {
      fLowMult = low;
      fUpMult = up;
    }
  }

  // Setters for pileup settings
  void SetPileupSelection(unsigned long plp = 1ULL) { fPlpType = plp; }
  void SetSPDPileupSelection(int cont = 3, float distance = 0.8)
  {
    fMinPlpContribSPD = cont;
    fMinPlpZdistSPD = distance;
  }
  void SetMVPileUpSelection(int cont = 5, float chi2 = 5, float wgthZdiff = 15., bool chkDiffBC = kFALSE)
  {
    fMinPlpContribMV = cont;
    fMaxPlpChi2MV = chi2;
    fMinWDistMV = wgthZdiff;
    fCheckPlpFromDifferentBCMV = chkDiffBC;
  }

  // Setters for track cut settings
  void SetMinSPDPoints(int np = 1) { fMinSPDPts = np; }
  void SetMinNdEdxSamples(int np = 3) { fMinNdEdxSamples = np; }
  void SetAbsEtaCut(double eta = .8) { fAbsEtaCut = eta; }
  void SetMinRapCut(double y = -.5) { fMinRapCut = y; }
  void SetMaxRapCut(double y = .5) { fMaxRapCut = y; }
  void SetCMSRapFct(double dy = .0) { fCMSRapFct = dy; }
  void SetMindEdx(double md = .0) { fMindEdx = md; }
  void SetMinNSigma(double ns = 1.5) { fMinNSigma = ns; }
  void SetMaxChi2Clu(double chi = 2.5) { fMaxChi2Clu = chi; }

  void SetDCACuts(double nsxy = 7., double nsz = 7.)
  {
    fNSigmaDCAxy = nsxy;
    fNSigmaDCAz = nsz;
  }

  // Setters for global settings
  void SetPidTech(int tech) { fPidMethod = static_cast<EPID_Type>(tech); }
  void SetITSPidParams(AliITSPidParams *pidParams) { fITSPidParams = pidParams; }
  void SetIsNominalBfield(bool flag = kTRUE) {fIsNominalBfield = flag;}
  void SetIsMC(bool flag = kTRUE) { fIsMC = flag; }
  void SetIsDCAUnfold(bool flag = kTRUE) { fIsDCAUnfoldHistoEnabled = flag; }
  void SetFillIntDistHist() { fFillIntDistHist = kTRUE; }
  void SetUseUnfolding(bool useUnfolding) {fUseUnfolding = useUnfolding; }

  AliEventCuts *GetAliEventCuts() { return &fEventCuts; }

  // Setters for smearing settings
  void SetSmearMC(double smearp, double smeardedx)
  {
    fSmearMC = kTRUE;
    fSmearP = smearp;
    fSmeardEdx = smeardedx;
  }

  // Default event selection cuts
  void SetupStandardEventCutsForRun1();
  void SetupEventCutsForRun1pPb();
  void SetupStandardEventCutsForRun2();

 protected:
  bool IsEventAccepted(EEvtCut_Type &evtSel);
  bool IsMultSelected();
  UInt_t IsPileup();
  bool IsGoodVtxZ();
  bool IsGoodSPDvtxRes(const AliESDVertex *spdVtx);
  bool SelectVertex2015pp();

  bool IsRapIn(double y) { return (y > fMinRapCut && y < fMaxRapCut) ? kTRUE : kFALSE; }
  bool DCAcut(double impactXY, double impactZ, double pt) const;
  bool DCAcutXY(double impactXY, double pt) const;
  bool DCAcutZ(double impactZ, double pt) const;

  double CookdEdx(double *s) const;
  double Eta2y(double pt, double m, double eta) const;

  void AnalyseMCParticles(AliMCEvent *lMCevent, EEvtCut_Type step, bool lHasGoodVtxGen);
  void CreateDCAcutFunctions();
  void PostAllData();

  double BetheITSsaHybrid(double p, double mass) const;
  int GetTrackPid(AliESDtrack *track, double *logdiff) const;
  int GetMostProbable(const double *pDens, const double *priors) const;
  void GetPriors(const AliVTrack *track, double *priors) const;
  float GetUnfoldedP(double dedx, float p) const;
  void ComputeBayesProbabilities(double *probs, const double *pDens, const double *prior);

 private:
  AliAnalysisTaskSEITSsaSpectra(const AliAnalysisTaskSEITSsaSpectra &source);
  AliAnalysisTaskSEITSsaSpectra &operator=(const AliAnalysisTaskSEITSsaSpectra &source);

  AliESDEvent *fESD;                  ///<  ESD object
  AliITSPidParams *fITSPidParams;     //-> parameterization used for PID
  AliITSPIDResponse *fITSPIDResponse; //!<! class with BB parameterizations

  // Standard event cut
  AliEventCuts fEventCuts; // basic cut variables for events

  /////////////////////////////////////////////////////
  // List
  /////////////////////////////////////////////////////
  TList *fOutput;     //!<! list with output
  TList *fListCuts;   //!<! list with functions storing DCA cut
  TList *fListTree;   //!<! list with trees
  TList *fListPriors; //!<! list with priors in case of bayesian pid approach

  TF1 *fDCAxyCutFunc; ///< function with DCAz cut vs. pt
  TF1 *fDCAzCutFunc;  ///< function with DCAxy cut vs. pt

  TNtuple *fNtupleData; //!<! output ntuple
  TNtuple *fNtupleMC;   //!<! output MC ntuple

  /////////////////////////////////////////////////////
  // General Histos
  /////////////////////////////////////////////////////
  TH2I *fHistNEvents;       //!<! histo with number of events / mult bin
  TH2I *fHistMCEvents;      //!<! histo with MC number of events / mult bin
  TH1F *fHistMultBefEvtSel; //!<! histo with multiplicity of the events before event selection
  TH1F *fHistMultAftEvtSel; //!<! histo with multiplicity of the events after all event selection
  TH2F *fHistVtxZ;          //!<! histo with the distribution of the primary vertex Z coordinate

  TH3F *fHistNTracks[kNchg];           //!<! histo with number of tracks vs Pt
  TH2F *fHistDEDXGen;                  //!<! histo with dedx versus momentum (generated, before track selection)
  TH2F *fHistDEDXGenposlabel;          //!<! histo with dedx versus momentum  with pos label (generated, before track selection)
  TH2F *fHistDEDXGenneglabel;          //!<! histo with dedx versus momentum  with neg label (generated, before track selection)
  TH2F *fHistDEDX;                     //!<! histo with dedx versus momentum
  TH2F *fHistDEDXdouble;               //!<! histo with dedx versus signed momentum
  TH2F *fHistDEDXposlabel;             //!<! histo with dedx versus momentum with positive label
  TH2F *fHistDEDXneglabel;             //!<! histo with dedx versus momentum with negative label
  TH2F *fHistNSigmaSep[kNchg * kNspc]; //!<! histo nsigma separation vs momentum
  TH2F *fHistSepPowerReco[kNchg * kNspc];  //!<!
  TH2F *fHistSepPowerTrue[kNchg * kNspc];  //!<!

  // MC histograms with spectra of primaries from the MC truth
  TH3F *fHistMCPart[kNchg * kNspc];            //!<! histo from events w/o gen Zvtx cut
  TH3F *fHistMCPartGoodGenVtxZ[kNchg * kNspc]; //!<! histo from events w/  gen Zvtx cut (<10cm)
  TH3F *fHistMCGenCharged;                     //!<! histo from events w/o gen Zvtx cut (with p instead of pt)

  // Reconstructed
  TH2F *fHistReco[kNchg * kNspc];         //!<! NSigma histos for 6 species
  THnSparseF *fHistRecoMC[kNchg * kNspc]; //!<! transverse momentum correlation with nsigma PID for 6 species
  THnSparseF *fHistRecoTrueMC[kNchg * kNspc]; //!<! transverse momentum correlation with true PID for 6 species
  THnSparseF *fHistRecoChargedMC; //!<! momentum correlation with true PID for 6 species
  THnSparseF *fHistMCDCA[kNchg * kNspc]; //!<! transverse momentum correlation for DCAxy unfolding

  //Rapidity distributions of identified particles
  TH2F *fHistYdist[kNchg * kNspc]; //!<! y distribution of identified (reco) particles
  TH2F *fHistYdistTruth[kNchg * kNspc]; //!<! y distribution of identified (reco) particle with MC truth for PID

  // MC histograms using reco values
  TH3F *fHistTruePIDMCReco[kNchg * kNspc]; //!<! histo with spectra of primaries from the MC truth (with pt reco)

  //MC histograms using gen values
  TH3F *fHistTruePIDMCGen[kNchg * kNspc]; //!<! histo with spectra of primaries from the MC truth (with pt generated)

  // DCAxy distributions
  TH3F *fHistDCAReco[kNchg * kNspc]; //!<! histo with DCA distibution

  // DCA Templates
  TH3F *fHistDCARecoPID_prim[kNchg * kNspc]; //!<! histo with DCA distibution and dedx PID
  TH3F *fHistDCARecoPID_sstr[kNchg * kNspc]; //!<! histo with DCA distibution and dedx PID
  TH3F *fHistDCARecoPID_smat[kNchg * kNspc]; //!<! histo with DCA distibution and dedx PID

  TH3F *fHistDCATruePID_prim[kNchg * kNspc]; //!<! histo with DCA distibution and MC truth PID
  TH3F *fHistDCATruePID_sstr[kNchg * kNspc]; //!<! histo with DCA distibution and MC truth PID
  TH3F *fHistDCATruePID_smat[kNchg * kNspc]; //!<! histo with DCA distibution and MC truth PID

  TH1F *fHistCharge[4]; //!<! histo with charge distribution to check the calibration

  // dEdx distributions
  TH2F *fHistPosHypPi; //! histo with DCA distibution in the kaons hypotesis (positive)
  TH2F *fHistPosHypKa; //! histo with DCA distibution in the kaons hypotesis (positive)
  TH2F *fHistPosHypPr; //! histo with DCA distibution in the protons hypotesis (positive)
  TH2F *fHistNegHypPi; //! histo with DCA distibution in the pions hypotesis (negative)
  TH2F *fHistNegHypKa; //! histo with DCA distibution in the kaons hypotesis (negative)
  TH2F *fHistNegHypPr; //! histo with DCA distibution in the protons hypotesis (negative)

  // dEdx distributions for MC
  TH2F *fHistMCPosOtherHypPion; //! histo with dedx using the MC truth
  TH2F *fHistMCPosOtherHypKaon; //! histo with dedx using the MC truth
  TH2F *fHistMCPosOtherHypProt; //! histo with dedx using the MC truth
  TH2F *fHistMCPosElHypPion;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosElHypKaon;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosElHypProt;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosPiHypPion;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosPiHypKaon;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosPiHypProt;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosKaHypPion;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosKaHypKaon;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosKaHypProt;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosPrHypPion;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosPrHypKaon;    //! histo with dedx using the MC truth
  TH2F *fHistMCPosPrHypProt;    //! histo with dedx using the MC truth

  TH2F *fHistMCNegOtherHypPion; //! histo with dedx using the MC truth
  TH2F *fHistMCNegOtherHypKaon; //! histo with dedx using the MC truth
  TH2F *fHistMCNegOtherHypProt; //! histo with dedx using the MC truth
  TH2F *fHistMCNegElHypPion;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegElHypKaon;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegElHypProt;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegPiHypPion;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegPiHypKaon;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegPiHypProt;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegKaHypPion;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegKaHypKaon;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegKaHypProt;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegPrHypPion;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegPrHypKaon;    //! histo with dedx using the MC truth
  TH2F *fHistMCNegPrHypProt;    //! histo with dedx using the MC truth

  TArrayD fCentBins;
  TArrayD fDCABins;
  TArrayD fPtBins;

  // evt sel.
  UInt_t fTriggerSel;
  double fMaxVtxZCut;
  bool fChkIsEventINELgtZERO;  // flag for select INEL>0 events
  bool fChkIsSDDIn;            // flag for counrint ev. in p-p 2.76
  bool fRejIncDAQ;             // flag for reject events with incomplete DAQ
  bool fDoSPDCvsTCut;          // flag for check compatibility between SPD clusters and tracklets
  bool fUseSelectVertex2015pp; // flag to select vertex based on criteria for 2015 pp data
  bool fChkVtxSPDRes;          // enable check on spd vtx resolution
  bool fChkVtxZSep;            // enable check on proximity of the z coordinate between both vertexer
  bool fReqBothVtx;            // ask for both trk and SPD vertex
  bool fExtEventCuts;          // enable use of AliEventCuts for event selection
  bool fUseUnfolding;          // enable if you want to use unfolding for PID
  // mult sel.
  unsigned int fMultMethod; // method for cent/mult values: 0=skip mult sel, 1=new cent framework, 2=old cent framework,
                            // 3=tracks+tracklets, 4=tracklets, 5=cluster on SPD
  TString fMultEstimator;   // centrality/multiplicity framework estimator name
  bool fMultEvSel;
  float fLowMult; // low Centrality cut
  float fUpMult;  // up  Centrality cut
  float fEvtMult; // event multiplicity -0.5 by default
  // Pileup selection setting
  unsigned long fPlpType;
  // PileupSPD settings
  int fMinPlpContribSPD; // minimum contributors to the pilup vertices, SPD
  float fMinPlpZdistSPD; // minimum distance for the SPD pileup vertex
  // PileupMV settings
  int fMinPlpContribMV;            // minimum contributors to the pilup vertices, multi-vertex
  float fMaxPlpChi2MV;             // minimum value of Chi2perNDF of the pileup vertex, multi-vertex
  float fMinWDistMV;               // minimum of the sqrt of weighted distance between the primary and the pilup vertex,
                                   // multi-vertex
  bool fCheckPlpFromDifferentBCMV; // pileup from different BC (Bunch Crossings)

  // trk sel.
  int fMinSPDPts;       // minimum number of SPD Points
  int fMinNdEdxSamples; // minimum number of SDD+SSD points
  double fAbsEtaCut;    // limits in pseudorap
  double fMinRapCut;
  double fMaxRapCut;
  double fCMSRapFct;
  double fMindEdx;     // minimum dE/dx value in a layer (to cut noise)
  double fMinNSigma;   // minimum number of sigmas
  double fMaxChi2Clu;  // maximum cluster
  double fNSigmaDCAxy; // DCA cut in bend. plane
  double fNSigmaDCAz;  // DCA cut along z

  // Global setting
  int fYear;            // FIXME Year (2009, 2010)
  EPID_Type fPidMethod; // track-by-track pid approach

  bool fUseDefaultPriors; // flag to use default(equal) priors
  bool fFillNtuple;       // flag to fill ntuples
  bool fIsMC;             // flag to switch on the MC analysis for the efficiency estimation
  bool fIsDCAUnfoldHistoEnabled; //flag to enable the filling of DCA histos used for unfolding
  bool fIsNominalBfield;  // flag to select the magnetic field (nominal = 0.5 T)
  bool fFillIntDistHist;  // flag to fill histogram with information for statistic pid analysis

  // smearing
  TRandom3 *fRandGener; // generator for smearing
  bool fSmearMC;        // flag to apply extra smearing on MC
  double fSmearP;       // extra relative smearing on simulated momentum
  double fSmeardEdx;    // extra relative smearing on simulated dE/dx

  //unfolding
  TH2F* fUnfProb[1170]; //-> histogram with unfolded matrices (probability)

  ClassDef(AliAnalysisTaskSEITSsaSpectra, 12);
};

#endif
