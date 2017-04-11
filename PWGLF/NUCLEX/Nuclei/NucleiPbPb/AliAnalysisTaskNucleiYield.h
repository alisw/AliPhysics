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
#include "AliNuclexEventCuts.h"

class TF1;
class TH2F;
class TH3F;
class AliFlowTrackCuts;
class AliAODTrack;
class AliVVertex;
class AliPIDResponse;
class TList;

class AliAnalysisTaskNucleiYield : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskNucleiYield(TString taskname = "NucleiYieldTask");
  virtual ~AliAnalysisTaskNucleiYield();

  void SetParticleType(AliPID::EParticleType part);
  void SetPDG (int pdg) { fPDG = (pdg >= 0) ? pdg : -pdg; }
  void SetIsMC (bool isMc) { fIsMC = isMc; }
  void SetFillOnlyEventHistos (bool onlyEventHistos) { fFillOnlyEventHistos = onlyEventHistos; }

  void SetRequireITSrefit (bool refit = true) { fRequireITSrefit = refit; }
  void SetRequireTPCrefit (bool refit = true) { fRequireTPCrefit = refit; }
  void SetRequireNoKinks (bool nokinks = true) { fRequireNoKinks = nokinks; }
  void SetRequireITSrecPoints (int rec = 4) { fRequireITSrecPoints = rec; }
  void SetRequireITSsignal (int sig = 3) { fRequireITSsignal = sig; }
  void SetRequireTPCsignal (int sig = 70) { fRequireITSsignal = sig; }
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
  void SetRequireVetoSPD (bool veto) { fRequireVetoSPD = veto; }
  void SetRequireMaxMomentum (float p) { fRequireMaxMomentum = p; }
  void SetEnablePtCorrection (bool cut) { fEnablePtCorrection = cut; }
  void SetDisableITSatHighPt (float pt) { fDisableITSatHighPt = pt; }
  void SetDisableTPCpidAtHighPt (float pt) { fDisableTPCpidAtHighPt = pt; }
  void SetFixForLHC14a6 (bool fix) { fFixForLHC14a6 = fix; }
  void SetForceMassAndZ(float mass, float z = 1) { fPDGMass = mass; fPDGMassOverZ = mass / z; }

  void SetCentBins (Int_t nbins, Float_t *bins);
  void SetDCABins (Int_t nbins, Float_t min, Float_t max);
  void SetDCABins (Int_t nbins, Float_t* bins);
  void SetPtBins (Int_t nbins, Float_t *bins);
  void SetCustomTPCpid (Float_t *par, Float_t sigma);
  void SetTOFBins (Int_t nbins, Float_t min, Float_t max);
  void SetDCAzBins (Int_t nbins, Float_t limit);
  void SetFlatteningProbabilities (Int_t n, Float_t *probs) { fFlatteningProbs.Set(n,probs); }
  void SetUseFlattening (bool useIt) { fEnableFlattening = useIt; }

  static int    GetNumberOfITSclustersPerLayer(AliVTrack *track, unsigned int &nSPD, unsigned int &nSDD, unsigned int &nSSD);
  static float  HasTOF(AliAODTrack *t, AliPIDResponse* pid);

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  AliNuclexEventCuts  fEventCut;
  TArrayD             fTOFfunctionPars;

  UInt_t              fFilterBit;       /// AOD filter bit for the tracks used in this analysis
private:
  AliAnalysisTaskNucleiYield (const AliAnalysisTaskNucleiYield &source);
  AliAnalysisTaskNucleiYield &operator=(const AliAnalysisTaskNucleiYield &source);

  bool   AcceptTrack(AliAODTrack *t, Double_t dca[2]);
  bool   PassesPIDSelection(AliAODTrack *t);
  float  GetTPCsigmas(AliVTrack *t);

  Bool_t Flatten(float cent);
  void PtCorrection(float &pt, bool positiveCharge);

  TF1                  *fTOFfunction;           //!<! TOF signal function

  TList                *fList;                  ///<  Output list
  TLorentzVector        fCutVec;                ///<  Vector used to perform some cuts
  Int_t                 fPDG;                   ///<  PDG code of the particle of interest
  Float_t               fPDGMass;               ///<  PDG mass
  Float_t               fPDGMassOverZ;          ///<  PDG mass over z
  float                 fCharge;                ///<  Charge of the particle of interest (absolute value)
  Bool_t                fIsMC;                  ///<  Switch between MC and data
  Bool_t                fFillOnlyEventHistos;   ///<  Set treu to fill only event related histograms

  AliPIDResponse       *fPID;                   //!<! PID response class
  Float_t               fMagField;              ///<  Magnetic field value for the current event

  Float_t               fDCAzLimit;             ///<  Limits of the \f$DCA_{z}\f$ histograms
  Int_t                 fDCAzNbins;             ///<  Number of bins used for \f$DCA_{z}\f$ distributions

  TArrayF               fPtCorrectionA;         ///<  Array containing the parametrisation of the \f$p_{T}\$ correction for anti-matter
  TArrayF               fPtCorrectionM;         ///<  Array containing the parametrisation of the \f$p_{T}\$ correction for matter

  Float_t               fTOFlowBoundary;        ///<  Lower limit for the TOF mass spectra histograms
  Float_t               fTOFhighBoundary;       ///<  Upper limit for the TOF mass spectra histograms
  Int_t                 fTOFnBins;              ///<  Number of bins used for the TOF mass spectra
  Float_t               fDisableITSatHighPt;    ///<  \f$p_{T}\f$ threshold for ITS cuts
  Float_t               fDisableTPCpidAtHighPt; ///<  \f$p_{T}\f$ threshold for TPC pid cut
  Bool_t                fEnablePtCorrection;    ///<  If true enables the MC based \f$p_{T}\f$ correction
  Bool_t                fRequireITSrefit;       ///<  Cut on tracks: set true to require ITS refit
  Bool_t                fRequireTPCrefit;       ///<  Cut on tracks: set true to require TPC refit
  Bool_t                fRequireNoKinks;        ///<  Cut on tracks: set true to exclude tracks from kink vertices
  UShort_t              fRequireITSrecPoints;   ///<  Cut on tracks: minimum number of required ITS recpoints
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
  Bool_t                fRequireVetoSPD;        ///<  Cut away all the tracks with at least 1 SPD cluster
  Float_t               fRequireMaxMomentum;    ///<  Cut in momentum for TPC only spectrum
  Bool_t                fFixForLHC14a6;         ///<  Switch on/off the fix for the MC centrality distribution

  Bool_t                fEnableFlattening;      ///<  Switch on/off the flattening

  AliPID::EParticleType fParticle;              ///<  Particle specie
  TArrayF               fCentBins;              ///<  Centrality bins
  TArrayF               fDCABins;               ///<  DCA bins
  TArrayF               fPtBins;                ///<  Transverse momentum bins
  TArrayF               fCustomTPCpid;          ///<  Custom parametrisation of the Bethe-Bloch
  TArrayF               fFlatteningProbs;       ///<  Flattening probabilities

  // Event related histograms
  TH1F                 *fFlattenedCentrality;   //!<! Events centrality distribution after the flattening
  TH1F                 *fCentralityClasses;     //!<! Events statistics per centrality classes

  // MC only histograms
  TH1F                 *fProduction;             //!<! *(MC only)* Total number of produced particles
  TH2F                 *fReconstructed[2][2];    //!<! *(MC only)* Positive and negative tracks reconstructed in the acceptance (ITS-TPC,ITS-TPC-TOF)
  TH2F                 *fTotal[2];               //!<! *(MC only)* Positively and negatively charged particles in acceptance
  TH2F                 *fPtCorrection[2];        //!<! *(MC only)* \f$p_{T}^{rec}-p_{T}^{MC}\f$ as a function of \f$p_{T}^{rec}\f$ for positive and negative tracks
  TH3F                 *fDCAPrimary[2][2];       //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of primaries
  TH3F                 *fDCASecondary[2][2];     //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of secondaries from material
  TH3F                 *fDCASecondaryWeak[2][2]; //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of secondaries from Weak Decay

  // Data histograms
  TH3F                 *fTOFsignal[2];           //!<! *(Data only)* TOF signal for (anti-)matter
  TH3F                 *fTPCcounts[2];           //!<! *(Data only)* TPC counts for (anti-)matter
  TH3F                 *fDCAxy[2][2];            //!<! *(Data only)* \f$DCA_{xy}\f$ distribution for ITS+TPC tracks
  TH3F                 *fDCAz[2][2];             //!<! *(Data only)* \f$DCA_{z}\f$ distribution for ITS+TPC tracks
  TH3F                 *fTOFtemplates[5];        //!<! *(Data only)* TOF signal templates for pi/k/p/d/t

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskNucleiYield, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskNucleiYield__) */
