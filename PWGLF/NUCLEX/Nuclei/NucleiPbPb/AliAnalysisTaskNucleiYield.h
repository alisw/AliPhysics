/// \class AliAnalysisTaskNucleiYield
/// \brief This task fills histograms required to perform the analysis on the light nuclei yield.
///
/// The histograms filled in this tasks are used by the analysis macro
/// ## Monte Carlo
/// There are mainly three items studied here:
/// * the acceptance x efficiency;
/// * the \f$DCA_{xy}\f$ distributions with respect to the primary vertex of primary and secondary;
/// * the difference of the reconstructed \f$p_{T}\f$ with respect to the MC truth \f$p_{T}\f$
///
/// ## Data
/// The histograms filled in this case are:
/// * TPC counts after PID cuts as a function of \f$p_{T}\f$ and centrality;
/// * TOF signal for selected tracks as a function of \f$p_{T}\f$ and centrality;
/// * the \f$DCA_{xy}\f$ distributions for selected tracks as a function of \f$p_{T}\f$ and centrality
///
/// \author Maximiliano Puccio <maximiliano.puccio@cern.ch>, University and INFN Torino
/// \date Feb 18, 2015

#ifndef __AliAnalysisTaskNucleiYield__
#define __AliAnalysisTaskNucleiYield__

#include "AliAnalysisTaskSE.h"
#include <Rtypes.h>
#include <TString.h>
#include <TArrayD.h>
#include <TArrayF.h>
#include <AliPIDResponse.h>
#include <AliPID.h>
#include <TLorentzVector.h>
#include "AliESDtrack.h"
#include "AliEventCuts.h"
#include <AliMCEvent.h>
#include <AliAODMCParticle.h>

#include <TH3F.h>
#include <vector>

#define LIGHT_SPEED 2.99792457999999984e-02 // in the units that TOF likes
#define EPS 1.e-16

class TF1;
class TH2F;
class TH3F;
class AliFlowTrackCuts;
class AliAODTrack;
class AliESDtrack;
class AliVVertex;
class AliPIDResponse;
class TList;
class TTree;
class AliPWGFunc;
class AliNanoAODTrack;

struct SLightNucleus {
  enum {
    kPrimary = BIT(0),
    kSecondaryMaterial = BIT(1),
    kSecondaryWeakDecay = BIT(2)
  };
  float pt;
  float eta;
  int   pdg;
  int   flag;
  char  centrality;
};

struct RLightNucleus {
  enum { 
    kT0fill = BIT(0),
    kPrimary = BIT(1),
    kSecondaryMaterial = BIT(2),
    kSecondaryWeakDecay = BIT(3),
    kHasTOF = BIT(4),
    kActiveLengthStatus = BIT(5),
    kCentral = BIT(6),
    kSemiCentral = BIT(7)
  };
  float pt;
  float eta;
  Double32_t dcaxy;          //[-2,2,10]
  Double32_t dcaz;           //[-2,2,10]
  Double32_t tofNsigma;      //[-12.8,12.8,12]
  Double32_t tpcNsigma;      //[-6.4,6.4,8]
  char       centrality;
  char       trackingPID;
  unsigned char tpcPIDcls;
  unsigned char flag;       //
};

class AliAnalysisTaskNucleiYield : public AliAnalysisTaskSE {
public:
  enum {
    kNoPtShape,
    kBlastWaveShape,
    kTsallisShape
  };

  AliAnalysisTaskNucleiYield(TString taskname = "NucleiYieldTask");
  virtual ~AliAnalysisTaskNucleiYield();


  void SetParticleType(AliPID::EParticleType part);
  void SetPDG (int pdg) { fPDG = (pdg >= 0) ? pdg : -pdg; }
  void SetIsMC (bool isMc) { fIsMC = isMc; }
  void SetFillOnlyEventHistos (bool onlyEventHistos) { fFillOnlyEventHistos = onlyEventHistos; }

  void SetRequireITSrecPoints (int rec = 4) { fRequireITSrecPoints = rec; }
  void SetRequireTPCrecPoints (int rec = 70) { fRequireTPCrecPoints = rec; }
  void SetRequireTPCfoundFraction (float rec = 0.8) { fRequireTPCfoundFraction = rec; }
  void SetRequireITSsignal (int sig = 3) { fRequireITSsignal = sig; }
  void SetRequireTPCsignal (int sig = 70) { fRequireTPCsignal = sig; }
  void SetRequireSDDrecPoints (int rec = 1) { fRequireSDDrecPoints = rec; }
  void SetRequireSPDrecPoints (int rec = 1) { fRequireSPDrecPoints = rec; }
  void SetEtaRange (float emin, float emax) { fRequireEtaMin = emin; fRequireEtaMax = emax; }
  void SetYRange (float ymin, float ymax) { fRequireYmin = ymin; fRequireYmax = ymax; }
  void SetRequireMaxChi2 (float maxChi2 = 4.f) { fRequireMaxChi2 = maxChi2; }
  void SetRequireMaxDCAxy (float maxDCA) { fRequireMaxDCAxy = maxDCA; }
  void SetRequireMaxDCAz (float maxDCA) { fRequireMaxDCAz = maxDCA; }
  void SetRequireTPCpidSigmas (float sig) { fRequireTPCpidSigmas = (sig > 0) ? sig : -sig; }
  void SetRequireITSpidSigmas (float sig) { fRequireITSpidSigmas = sig; }
  void SetRequireTOFpidSigmas (float sig) { fRequireTOFpidSigmas = sig; }
  void SetRequireMinEnergyLoss (float ecut) { fRequireMinEnergyLoss = ecut; }
  void SetTPCActiveLengthCut (bool apply = true, float dzone = 3.0, float length = 130, float gcut1 = 1.5, float gcut2 = 0.85, float gcut3 = 0.7) { 
    SetApplyTPCActiveLengthCut(apply);
    SetRequireDeadZoneWidth(dzone);
    SetRequireCutGeoNcrNclLength(length);
    SetRequireCutGeoNcrNclGeom1Pt(gcut1);
    SetRequireCutGeoNcrNclFractionNcr(gcut2);
    SetRequireCutGeoNcrNclFractionNcl(gcut3);
  }
  void SetApplyTPCActiveLengthCut (bool apply) { fApplyTPCLengthCut = apply; }
  void SetRequireDeadZoneWidth (float dzone) { fRequireDeadZoneWidth = dzone; }
  void SetRequireCutGeoNcrNclLength (float length) { fRequireCutGeoNcrNclLength = length; }
  void SetRequireCutGeoNcrNclGeom1Pt (float gcut) { fRequireCutGeoNcrNclGeom1Pt = gcut; }
  void SetRequireCutGeoNcrNclFractionNcr (float gcut) { fCutGeoNcrNclFractionNcr = gcut; }
  void SetRequireCutGeoNcrNclFractionNcl (float gcut) { fCutGeoNcrNclFractionNcl = gcut; }
  void SetRequireVetoSPD (bool veto) { fRequireVetoSPD = veto; }
  void SetRequireMaxMomentum (float p) { fRequireMaxMomentum = p; }
  void SetEnablePtCorrection (bool cut) { fEnablePtCorrection = cut; }
  void SetDisableITSatHighPt (float pt) { fDisableITSatHighPt = pt; }
  void SetDisableTPCpidAtHighPt (float pt) { fDisableTPCpidAtHighPt = pt; }
  void SetFixForLHC14a6 (bool fix) { fFixForLHC14a6 = fix; }
  void SetForceMassAndZ(float mass, float z = 1) { fPDGMass = mass; fPDGMassOverZ = mass / z; }
  void SetITSelectronRejection(float nsigma = 2.) { fITSelectronRejectionSigma = nsigma; }

  void SetCentBins (Int_t nbins, Float_t *bins);
  void SetDCABins (Int_t nbins, Float_t min, Float_t max);
  void SetDCABins (Int_t nbins, Float_t* bins);
  void SetPtBins (Int_t nbins, Float_t min, Float_t max);
  void SetPtBins (Int_t nbins, Float_t *bins);
  void SetCustomTPCpid (Float_t *par, Float_t sigma);
  void SetTOFBins (Int_t nbins, Float_t min, Float_t max);
  void SetDCAzBins (Int_t nbins, Float_t limit);
  void SetSigmaBins (Int_t nbins, Float_t limit);
  void SetTOFSigmaBins (Int_t nbins, Float_t limit);
  void SetFlatteningProbabilities (Int_t n, Float_t *probs) { fFlatteningProbs.Set(n,probs); }
  void SetUseFlattening (bool useIt) { fEnableFlattening = useIt; }
  void SetPtWeightingFunction (int functionID, int nPars, float* pars) {
    fPtShapeFunction = functionID;
    fPtShapeParams.Set(nPars,pars);
  }
  void SetBeamRapidity(float rap) { fBeamRapidity = rap; }
  void SetCentralityEstimator(int est) { fEstimator = est; }
  void SetupTRDstudies(int vintage, bool trdin) { fTRDvintage = vintage; fTRDin = trdin; }

  void SaveTrees(bool save=true) { fSaveTrees = save; }
  void SetTOFminPtTrees(float pt) { fTOFminPtTrees = pt; }

  static int    GetNumberOfITSclustersPerLayer(AliVTrack *track, int &nSPD, int &nSDD, int &nSSD);
  static float  HasTOF(AliVTrack *t, AliPIDResponse* pid);
  static float  HasTOF(AliNanoAODTrack *t, AliPIDResponse* pid);

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  AliEventCuts  fEventCut;
  TArrayD       fTOFfunctionPars;

  UInt_t              fFilterBit;       /// AOD filter bit for the tracks used in this analysis (set to 0 to skip the cut)
  bool                fPropagateTracks; /// Workaround for troublesome productions
  TArrayF             fPtCorrectionA;         ///<  Array containing the parametrisation of the \f$p_{T}\$ correction for anti-matter
  TArrayF             fPtCorrectionM;         ///<  Array containing the parametrisation of the \f$p_{T}\$ correction for matter
  double              fOptionalTOFcleanup;   ///<

  std::vector<float>  fINT7intervals;        ///< Array containing the centrality interval where we select only kINT7 triggers  
private:
  AliAnalysisTaskNucleiYield (const AliAnalysisTaskNucleiYield &source);
  AliAnalysisTaskNucleiYield &operator=(const AliAnalysisTaskNucleiYield &source);

  template<class track> void TrackLoop(track* t, bool nano);
  void SetSLightNucleus(AliAODMCParticle* part, SLightNucleus& snucl);

  bool IsInTRD(float pt, float phi, float sign);

  template <class track_t>
  bool   AcceptTrack(track_t *t, Float_t dca[2]);
  int   PassesPIDSelection(AliAODTrack *t);
  int   PassesPIDSelection(AliNanoAODTrack *t);
  float  GetTPCsigmas(AliVTrack *t);
  float  GetTOFsigmas(AliVTrack* t);

  Bool_t Flatten(float cent);
  void PtCorrection(float &pt, bool positiveCharge);
  Bool_t IsSelectedTPCGeoCut(AliAODTrack *track);
  Bool_t IsSelectedTPCGeoCut(AliNanoAODTrack *track);


  TString               fCurrentFileName;       ///<  Currently analysed file name
  TF1                  *fTOFfunction;           //!<! TOF signal function

  TList                *fList;                  ///<  Output list
  TTree                *fRTree;                 ///<  Output reconstructed ttree
  TTree                *fSTree;                 ///<  Output simulated ttree
  TLorentzVector        fCutVec;                ///<  Vector used to perform some cuts
  Int_t                 fPDG;                   ///<  PDG code of the particle of interest
  Float_t               fPDGMass;               ///<  PDG mass
  Float_t               fPDGMassOverZ;          ///<  PDG mass over z
  float                 fCharge;                ///<  Charge of the particle of interest (absolute value)
  Bool_t                fIsMC;                  ///<  Switch between MC and data
  Bool_t                fFillOnlyEventHistos;   ///<  Set treu to fill only event related histograms

  AliPIDResponse       *fPID;                   //!<! PID response class
  ULong64_t             fTriggerMask;           //!<  Trigger Mask of the Event
  Float_t               fMagField;              //!<! Magnetic field value for the current event
  Float_t               fCentrality;            //!<! Centrality for the current event
  std::vector<int>      fRejectedParticles;     //!<! List of rejected particles for adapting the MC pt shape

  Float_t               fDCAzLimit;             ///<  Limits of the \f$DCA_{z}\f$ histograms
  Int_t                 fDCAzNbins;             ///<  Number of bins used for \f$DCA_{z}\f$ distributions
  Float_t               fSigmaLimit;            ///<  Limits of the \f$n_{sigma}\f$ histograms
  Int_t                 fSigmaNbins;            ///<  Number of bins used for \f$n_{sigma}\f$ distributions
  Float_t               fTOFSigmaLimit;         ///<  Limits of the \f$n_{sigma_{TOF}}\f$ histograms
  Int_t                 fTOFSigmaNbins;         ///<  Number of bins used for \f$n_{sigma_{TOF}}\f$ distributions

  Float_t               fTOFlowBoundary;        ///<  Lower limit for the TOF mass spectra histograms
  Float_t               fTOFhighBoundary;       ///<  Upper limit for the TOF mass spectra histograms
  Int_t                 fTOFnBins;              ///<  Number of bins used for the TOF mass spectra
  Float_t               fDisableITSatHighPt;    ///<  \f$p_{T}\f$ threshold for ITS cuts
  Float_t               fDisableTPCpidAtHighPt; ///<  \f$p_{T}\f$ threshold for TPC pid cut
  Bool_t                fEnablePtCorrection;    ///<  If true enables the MC based \f$p_{T}\f$ correction
  UShort_t              fRequireITSrecPoints;   ///<  Cut on tracks: minimum number of required ITS recpoints
  UShort_t              fRequireTPCrecPoints;   ///<  Cut on tracks: minimum number of required ITS recpoints
  UShort_t              fRequireITSsignal;      ///<  Cut on tracks: minimum number of required ITS PID recpoints
  UShort_t              fRequireSDDrecPoints;   ///<  Cut on tracks: minimum number of required SDD recpoints
  UShort_t              fRequireSPDrecPoints;   ///<  Cut on tracks: minimum number of required SPD recpoints
  UShort_t              fRequireTPCsignal;      ///<  Cut on tracks: minimum number of required TPC PID recpoints
  Float_t               fRequireEtaMin;         ///<  Cut on tracks: minimum eta for the track
  Float_t               fRequireEtaMax;         ///<  Cut on tracks: maximum eta for the track
  Float_t               fRequireYmin;           ///<  Cut on tracks: mimimum y for the track (using PDG mass)
  Float_t               fRequireYmax;           ///<  Cut on tracks: maximum y for the track (using PDG mass)
  Float_t               fRequireMaxChi2;        ///<  Cut on tracks: maximum TPC \f$\chi^{2}/NDF\f$
  Float_t               fRequireMaxDCAxy;       ///<  Cut on tracks: maximum \f$DCA_{xy}\f$ for the track
  Float_t               fRequireMaxDCAz;        ///<  Cut on tracks: maximum \f$DCA_{z}\f$ for the track
  Float_t               fRequireTPCpidSigmas;   ///<  Cut on TPC PID number of sigmas
  Float_t               fRequireITSpidSigmas;   ///<  Cut on ITS PID number of sigmas
  Float_t               fRequireTOFpidSigmas;   ///<  Cut on ITS PID number of sigmas
  Float_t               fRequireMinEnergyLoss;  ///<  Cut on the minimum energy loss counts in TPC
  Bool_t                fApplyTPCLengthCut;     ///<  Apply TPC Active Length cut in Accept track check.
  Float_t               fRequireDeadZoneWidth;  ///<  Cut on on TPC Geometrical Selection Deadzone width
  Float_t               fRequireCutGeoNcrNclLength; ///<  Cut on TPC Geometrical Selection Length
  Float_t               fRequireCutGeoNcrNclGeom1Pt; ///<  Cut on TPC Geometrical Selection 1 Pt
  Float_t               fCutGeoNcrNclFractionNcr; ///<  Cut on TPC Geometrical Selection Fraction
  Float_t               fCutGeoNcrNclFractionNcl; ///<  Cut on TPC Geometrical Selection NCluster
  Bool_t                fRequireVetoSPD;        ///<  Cut away all the tracks with at least 1 SPD cluster
  Float_t               fRequireMaxMomentum;    ///<  Cut in momentum for TPC only spectrum
  Bool_t                fFixForLHC14a6;         ///<  Switch on/off the fix for the MC centrality distribution
  Float_t               fRequireTPCfoundFraction; ///< Found over findable clusters
  Int_t                 fPtShapeFunction;       ///<  Id of the function used to weight the MC input pt shape (see the enum)
  Float_t               fPtShapeMaximum;        ///<  Maximum of the pt shape used
  Float_t               fITSelectronRejectionSigma; ///< nSigma rejection band in ITS response around the electron band for TPC only analysis
  Float_t               fBeamRapidity;          ///< Beam rapidity in case of asymmetric colliding systems
  Int_t                 fEstimator;             ///< Choose the centrality estimator from AliEventCuts

  Bool_t                fEnableFlattening;      ///<  Switch on/off the flattening

  Bool_t                fSaveTrees;             ///<  Switch on/off the output TTrees
  Float_t               fTOFminPtTrees;         ///<  Pt after which the TOF pid is required to save the tree
  RLightNucleus         fRecNucleus;            ///<  Reconstructed nucleus
  SLightNucleus         fSimNucleus;            ///<  Simulated nucleus


  AliPID::EParticleType fParticle;              ///<  Particle specie
  TArrayF               fCentBins;              ///<  Centrality bins
  TArrayF               fDCABins;               ///<  DCA bins
  TArrayF               fPtBins;                ///<  Transverse momentum bins
  TArrayF               fCustomTPCpid;          ///<  Custom parametrisation of the Bethe-Bloch
  TArrayF               fFlatteningProbs;       ///<  Flattening probabilities
  TArrayF               fPtShapeParams;         ///<  Params used by the pt shape function

  AliPWGFunc*           fFunctCollection;       //!<! Collection of functions
  TF1*                  fPtShape;               //!<! Function used to model the pt shape in MC

  // Event related histograms
  TH2F                 *fNormalisationHist;     //!<! Normalisation per centrality classes

  // MC only histograms
  TH1F                 *fProduction;             //!<! *(MC only)* Total number of produced particles
  TH2F                 *fReconstructed[2][2];    //!<! *(MC only)* Positive and negative tracks reconstructed in the acceptance (ITS-TPC,ITS-TPC-TOF)
  TH2F                 *fTotal[2];               //!<! *(MC only)* Positively and negatively charged particles in acceptance
  TH2F                 *fPtCorrection[2];        //!<! *(MC only)* \f$p_{T}^{rec}-p_{T}^{MC}\f$ as a function of \f$p_{T}^{rec}\f$ for positive and negative tracks
  TH2F                 *fPcorrectionTPC[2];      //!<! *(MC only)* \f$p_{T}^{rec}-p_{T}^{MC}\f$ as a function of \f$p_{T}^{rec}\f$ for positive and negative tracks in TPC
  TH3F                 *fDCAPrimary[2][2];       //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of primaries
  TH3F                 *fDCASecondary[2][2];     //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of secondaries from material
  TH3F                 *fDCASecondaryWeak[2][2]; //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of secondaries from Weak Decay

  // Data histograms
  TH3F                 *fMultDistributionTPC;    //!<! *(Data only)* Multiplicity distribution for antimatter
  TH3F                 *fMultDistributionTOF;    //!<! *(Data only)* Multiplicity distribution for antimatter
  TH3F                 *fTOFnSigma[2];           //!<! *(Data only)* TOF nSigma counts for (anti-)matter
  TH3F                 *fTOFT0FillNsigma[2];     //!<! *(Data only)* TOF nSigma counts for (anti-)matter
  TH3F                 *fTOFNoT0FillNsigma[2];   //!<! *(Data only)* TOF nSigma counts for (anti-)matter
  TH3F                 *fTOFsignal[2];           //!<! *(Data only)* TOF signal for (anti-)matter
  TH3F                 *fTOFT0FillSignal[2];     //!<! *(Data only)* TOF signal for (anti-)matter
  TH3F                 *fTOFNoT0FillSignal[2];   //!<! *(Data only)* TOF signal for (anti-)matter
  TH3F                 *fTPCcounts[2];           //!<! *(Data only)* TPC counts for (anti-)matter
  TH3F                 *fTPCsignalTpl[2];        //!<! *(Data only)* TPC counts for (anti-)matter
  TH3F                 *fTPCbackgroundTpl[2];    //!<! *(Data only)* TPC counts for (anti-)matter
  TH3F                 *fDCAxy[2][2];            //!<! *(Data only)* \f$DCA_{xy}\f$ distribution for ITS+TPC tracks
  TH3F                 *fDCAz[2][2];             //!<! *(Data only)* \f$DCA_{z}\f$ distribution for ITS+TPC tracks
  TH2F *fHist2Phi[2]; //! phi vs pt, negative (0) and positive (1): used for monitoring
  TF1 *fTRDboundariesPos[4]; //!<! Function with the phi limits of TRD boundaries as a function of pt
  TF1 *fTRDboundariesNeg[4]; //!<! Function with the phi limits of TRD boundaries as a function of pt
  int  fTRDvintage;          /// TRD configuration (year)
  bool fTRDin;               /// if true only tracks within TRD area are considered
  int  fNanoPIDindexTPC;
  int  fNanoPIDindexTOF;
  int  fRefMult;

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskNucleiYield, 1);
  /// \endcond
};

template<class track_t> void AliAnalysisTaskNucleiYield::TrackLoop(track_t* track, bool nano) {

  if (track->GetID() < 0) return;
  Float_t dca[2] = {0.f,0.f};
  if (!track->TestFilterBit(fFilterBit) && fFilterBit) return;
  bool acceptedTrack = AcceptTrack(track,dca);
  float beta = HasTOF(track,fPID);
  if (beta > 1. - EPS) beta = -1;
  const float m2 = track->P() * track->P() * (1.f / (beta * beta) - 1.f);

  if (!acceptedTrack) return;

  if (fSaveTrees && track->Pt() < 10.) {
    //double mcPt = 0;
    bool good2save{true};
    if (fIsMC) {
      int mcId = std::abs(track->GetLabel());
      AliAODMCParticle *part = (AliAODMCParticle*)MCEvent()->GetTrack(mcId);
      if (part) {
        good2save = std::abs(part->PdgCode()) == fPDG;
        SetSLightNucleus(part, fSimNucleus);
      } else
        good2save = false;
    }
    if (good2save) {
      AliTOFPIDResponse& tofPID = fPID->GetTOFResponse();

      fRecNucleus.pt = track->Pt() * track->Charge();
      fRecNucleus.eta = track->Eta();
      fRecNucleus.dcaxy = dca[0];
      fRecNucleus.dcaz = dca[1];
      fRecNucleus.tpcNsigma = GetTPCsigmas(track);
      fRecNucleus.tofNsigma = GetTOFsigmas(track);
      fRecNucleus.centrality = fCentrality;
      fRecNucleus.tpcPIDcls = track->GetTPCsignalN();
      fRecNucleus.flag = 0;
      fRecNucleus.flag |= !tofPID.GetT0binMask(tofPID.GetMomBin(track->GetTPCmomentum())) ? RLightNucleus::kT0fill : 0;
      if (fIsMC) {
        fRecNucleus.flag |= (fSimNucleus.flag == SLightNucleus::kPrimary) ? RLightNucleus::kPrimary : 0;
        fRecNucleus.flag |= (fSimNucleus.flag == SLightNucleus::kSecondaryWeakDecay) ? RLightNucleus::kSecondaryWeakDecay : 0;
        fRecNucleus.flag |= (fSimNucleus.flag == SLightNucleus::kSecondaryMaterial) ? RLightNucleus::kSecondaryMaterial : 0;
      }
      fRecNucleus.trackingPID = track->GetPIDForTracking();
      fRecNucleus.flag |= (beta > EPS) ? RLightNucleus::kHasTOF : 0;
      fRecNucleus.flag |= (IsSelectedTPCGeoCut(track)) ? RLightNucleus::kActiveLengthStatus : 0;
      fRecNucleus.flag |= (fTriggerMask & AliVEvent::kCentral) ? RLightNucleus::kCentral : 0;
      fRecNucleus.flag |= (fTriggerMask & AliVEvent::kSemiCentral) ? RLightNucleus::kSemiCentral : 0;
      if (std::abs(fRecNucleus.tpcNsigma) < 6.4 && (track->Pt() < fTOFminPtTrees || std::abs(fRecNucleus.tofNsigma) < 6.4))
        fRTree->Fill();
    }
  }

  bool positive = track->Charge() > 0;
  if (fHist2Phi[positive]) fHist2Phi[positive]->Fill(track->Phi() , track->Pt() );

  const int iTof = beta > EPS ? 1 : 0;
  float pT = track->Pt() * fCharge;
  float p_TPC = track->GetTPCmomentum();
  float p = track->P(); 
  int pid_mask = PassesPIDSelection(track);
  bool pid_check = (pid_mask & 7) == 7;
  if (fEnablePtCorrection) PtCorrection(pT,track->Charge() > 0);

  int mcId = TMath::Abs(track->GetLabel());
  if (fPtShape) {
    if (std::find(fRejectedParticles.begin(), fRejectedParticles.end(), mcId) != fRejectedParticles.end()) {
      return;
    }
  }
  if (fIsMC) {
    AliAODMCParticle *part = (AliAODMCParticle*)MCEvent()->GetTrack(mcId);
    /// Workaround: if the AOD are filtered with an AliRoot tag before v5-08-18, hyper-nuclei prongs
    /// are marked as SecondaryFromMaterial.
    const int mother_id = part->GetMother();
    AliAODMCParticle* mother = (mother_id >= 0) ? (AliAODMCParticle*)MCEvent()->GetTrack(mother_id) : 0x0;
    const int mother_pdg = mother ? TMath::Abs(mother->PdgCode()) : 0;
    const bool isFromHyperNucleus = (mother_pdg > 1000000000 && (mother_pdg / 10000000) % 10 != 0);
    if (!part) return;
    const int iC = part->Charge() > 0 ? 1 : 0;
    if (std::abs(part->PdgCode()) == fPDG) {
      for (int iR = iTof; iR >= 0; iR--) {
        if (part->IsPhysicalPrimary()) {
          if (TMath::Abs(dca[0]) <= fRequireMaxDCAxy &&
              (iR || fRequireMaxMomentum < 0 || track->GetTPCmomentum() < fRequireMaxMomentum) &&
              (!iR || pid_check) && (iR || pid_mask & 8))
            fReconstructed[iR][iC]->Fill(fCentrality,pT);
          fDCAPrimary[iR][iC]->Fill(fCentrality,pT,dca[0]);
          if (!iR) {
            fPtCorrection[iC]->Fill(pT,part->Pt()-pT);
            fPcorrectionTPC[iC]->Fill(p_TPC,part->P()-p_TPC);
            } // Fill it only once.
        } else if (part->IsSecondaryFromMaterial() && !isFromHyperNucleus)
          fDCASecondary[iR][iC]->Fill(fCentrality,pT,dca[0]);
        else
          fDCASecondaryWeak[iR][iC]->Fill(fCentrality,pT,dca[0]);
      }
    }
  } else {
    const int iC = (track->Charge() > 0) ? 1 : 0;

    float tpc_n_sigma = GetTPCsigmas(track);
    float tof_n_sigma = iTof ? GetTOFsigmas(track) : -999.f;

    if (iC == 0) {
      if (std::abs(tpc_n_sigma) < 5.)
        fMultDistributionTPC->Fill(fRefMult, pT, tpc_n_sigma);
      if (std::abs(tpc_n_sigma) < 3.5 && std::abs(tof_n_sigma) < 5.)
        fMultDistributionTOF->Fill(fRefMult, pT, tof_n_sigma);
    }

    for (int iR = iTof; iR >= 0; iR--) {

      /// TPC asymmetric cut to avoid contamination from protons in the DCA distributions. TOF sigma cut is set to 4
      /// to compensate for the shift in the sigma (to be rechecked in case of update of TOF PID response)
      if (tpc_n_sigma > -2. && tpc_n_sigma < 3. && (fabs(tof_n_sigma) < 4. || !iTof)) {
        fDCAxy[iR][iC]->Fill(fCentrality, pT, dca[0]);
        fDCAz[iR][iC]->Fill(fCentrality, pT, dca[1]);
      }
    }
    if (TMath::Abs(dca[0]) > fRequireMaxDCAxy) return;
    bool tofCleanUp = fOptionalTOFcleanup < 0 ? true : std::abs(tof_n_sigma) < fOptionalTOFcleanup || beta < 0;
    if ((fRequireMaxMomentum < 0 || track->GetTPCmomentum() < fRequireMaxMomentum) && tofCleanUp &&
        (pid_mask & 8))
      fTPCcounts[iC]->Fill(fCentrality, pT, tpc_n_sigma);

    if (iTof == 0) return;
    if (std::abs(tof_n_sigma) < 4.)
      fTPCsignalTpl[iC]->Fill(fCentrality, pT, tpc_n_sigma);
    else
      fTPCbackgroundTpl[iC]->Fill(fCentrality, pT, tpc_n_sigma);

    if (!pid_check) return;
    /// \f$ m = \frac{p}{\beta\gamma} \f$
    fTOFsignal[iC]->Fill(fCentrality, pT, m2 - fPDGMassOverZ * fPDGMassOverZ);
    fTOFnSigma[iC]->Fill(fCentrality, pT, tof_n_sigma);

    if (fPID->GetTOFResponse().GetStartTimeMask(p) == 0) {
      fTOFT0FillSignal[iC]->Fill(fCentrality, pT, m2 - fPDGMassOverZ * fPDGMassOverZ);
      fTOFT0FillNsigma[iC]->Fill(fCentrality, pT, tof_n_sigma);
    }
    else { 
      fTOFNoT0FillSignal[iC]->Fill(fCentrality, pT, m2 - fPDGMassOverZ * fPDGMassOverZ);
      fTOFNoT0FillNsigma[iC]->Fill(fCentrality, pT, tof_n_sigma);
    }

  }
}

/// This function checks whether a track passes the cuts required in this task
///
/// \param track Track that is going to be checked
/// \param dca[2] Projections on the transverse plane and on z of the distance of closest approach
///               of the track to the primary vertex
/// \return Boolean value: true means that the track has passed all the cuts.
///
template<class track_t>
bool AliAnalysisTaskNucleiYield::AcceptTrack(track_t *track, Float_t dca[2]) {
  fCutVec.SetPtEtaPhiM(track->Pt() * fCharge, track->Eta(), track->Phi(), fPDGMass);
  if (track->Eta() < fRequireEtaMin || track->Eta() > fRequireEtaMax) return false;
  if (fCutVec.Rapidity() < fRequireYmin + fBeamRapidity || fCutVec.Rapidity() > fRequireYmax + fBeamRapidity) return false;
  if (track->Chi2perNDF() > fRequireMaxChi2) return false;
  if (track->GetTPCNcls() < fRequireTPCrecPoints) return false;
  if (track->GetTPCFoundFraction() < fRequireTPCfoundFraction) return false;
  if (track->GetTPCsignalN() < fRequireTPCsignal) return false;
  if (track->GetTPCsignal() < fRequireMinEnergyLoss) return false;
  if (fTRDvintage != 0 && fTRDin != IsInTRD(track->Pt(), track->Phi(), track->Charge())) return false;
  if (fApplyTPCLengthCut && fRequireCutGeoNcrNclGeom1Pt > 0 && fRequireCutGeoNcrNclLength > 0 && fRequireDeadZoneWidth > 0) {
    if(!IsSelectedTPCGeoCut(track)) return false;
  }

  /// ITS related cuts
  dca[0] = 0.;
  dca[1] = 0.;
  if (track->Pt() < fDisableITSatHighPt) {
    int nSPD = 0u, nSDD = 0u, nSSD = 0u;
    int nITS = GetNumberOfITSclustersPerLayer(track, nSPD, nSDD, nSSD);
    if (nITS < fRequireITSrecPoints) return false;
    if (nSPD < fRequireSPDrecPoints) return false;
    if (nSDD < fRequireSDDrecPoints) return false;
    if (fRequireVetoSPD && nSPD > 0) return false;
    track->GetImpactParameters(dca[0], dca[1]);
    if (TMath::Abs(dca[1]) > fRequireMaxDCAz) return false;
    //if (TMath::Abs(dca[0]) > fRequireMaxDCAxy) return false;
  }

  return true;
}


#endif /* defined(__AliAnalysisTaskNucleiYield__) */
