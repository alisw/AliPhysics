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
#include <TArrayF.h>
#include <AliPIDResponse.h>
#include <AliPID.h>

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
  void SetRequireTPCrecPoints (int rec = 70) { fRequireITSrecPoints = rec; }
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
  void SetRequireMagneticField (int cut) { fRequireMagneticField = cut; }
  void SetRequireVetoSPD (bool veto) { fRequireVetoSPD = veto; }
  void SetRequireMaxMomentum (float p) { fRequireMaxMomentum = p; }
  void SetEnablePerformancePlot (bool cut) { fEnablePerformance = cut; }
  void SetEnablePtCorrection (bool cut) { fEnablePtCorrection = cut; }
  void SetDisableITSatHighPt (float pt) { fDisableITSatHighPt = pt; }
  void SetDisableTPCpidAtHighPt (float pt) { fDisableTPCpidAtHighPt = pt; }
  void SetFixForLHC14a6 (bool fix) { fFixForLHC14a6 = fix; }
  void SetRequireTrackLength(float len) { fRequireTrackLength = len; }
  void SetForceMassAndZ(float mass, float z = 1) { fPDGMass = mass; fPDGMassOverZ = mass / z; }

  void SetCentralityLimits(float min, float max) { fRequireMinCentrality = min; fRequireMaxCentrality = max; }

  void SetCentBins (Int_t nbins, Float_t *bins);
  void SetDCABins (Int_t nbins, Float_t min, Float_t max);
  void SetDCABins (Int_t nbins, Float_t* bins);
  void SetPtBins (Int_t nbins, Float_t *bins);
  void SetCustomTPCpid (Float_t *par, Float_t sigma);
  void SetTOFBins (Int_t nbins, Float_t min, Float_t max);
  void SetDCAzBins (Int_t nbins, Float_t limit);
  void SetFlatteningProbabilities (Int_t n, Float_t *probs) { fFlatteningProbs.Set(n,probs); }
  void SetPhiRegions (bool pos,int n, float *regions) { fPhiRegions[int(pos)].Set(n,regions); }
  void SetUseNewCentralityFramework (bool useIt) { fNewCentralityFramework = useIt; }
  void SetUseFlattening (bool useIt) { fEnableFlattening = useIt; }
  void SetTriggerMask (ULong_t mask) { fTriggerMask = mask; }
  
  void SetEnableLogarithmicBinning (bool useit) { fEnableLogAxisInPerformancePlots = useit; }

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

private:
  AliAnalysisTaskNucleiYield (const AliAnalysisTaskNucleiYield &source);
  AliAnalysisTaskNucleiYield &operator=(const AliAnalysisTaskNucleiYield &source);

  Bool_t AcceptTrack(AliAODTrack *t, Double_t dca[2]);
  Bool_t PassesPIDSelection(AliAODTrack *t);
  Float_t HasTOF(AliAODTrack *t);

  Bool_t Flatten(float cent);
  void PtCorrection(float &pt, bool positiveCharge);

  TList                *fList;                  ///<  Output list
  Int_t                 fPDG;                   ///<  PDG code of the particle of interest
  Float_t               fPDGMass;               ///<  PDG mass
  Float_t               fPDGMassOverZ;          ///<  PDG mass over z
  Bool_t                fIsMC;                  ///<  Switch between MC and data
  Bool_t                fFillOnlyEventHistos;   ///<  Set treu to fill only event related histograms

  AliPIDResponse       *fPID;                   //!<! PID response class
  Float_t               fMagField;              ///<  Magnetic field value for the current event
  AliVVertex           *fPrimaryVertex;         //!<! Primary vertex of the current event

  Float_t               fDCAzLimit;             ///<  Limits of the \f$DCA_{z}\f$ histograms
  Int_t                 fDCAzNbins;             ///<  Number of bins used for \f$DCA_{z}\f$ distributions

  TArrayF               fPtCorrectionA;         ///<  Array containing the parametrisation of the \f$p_{T}\$ correction for anti-matter
  TArrayF               fPtCorrectionM;         ///<  Array containing the parametrisation of the \f$p_{T}\$ correction for matter

  Float_t               fTOFlowBoundary;        ///<  Lower limit for the TOF mass spectra histograms
  Float_t               fTOFhighBoundary;       ///<  Upper limit for the TOF mass spectra histograms
  Int_t                 fTOFnBins;              ///<  Number of bins used for the TOF mass spectra
  Float_t               fDisableITSatHighPt;    ///<  \f$p_{T}\f$ threshold for ITS cuts
  Float_t               fDisableTPCpidAtHighPt; ///<  \f$p_{T}\f$ threshold for TPC pid cut
  Bool_t                fEnablePerformance;     ///<  If true enables the task saves in the output the performance plots
  Bool_t                fEnablePtCorrection;    ///<  If true enables the MC based \f$p_{T}\f$ correction
  Bool_t                fRequireITSrefit;       ///<  Cut on tracks: set true to require ITS refit
  Bool_t                fRequireTPCrefit;       ///<  Cut on tracks: set true to require TPC refit
  Bool_t                fRequireNoKinks;        ///<  Cut on tracks: set true to exclude tracks from kink vertices
  UShort_t              fRequireITSrecPoints;   ///<  Cut on tracks: minimum number of required ITS recpoints
  UShort_t              fRequireITSsignal;      ///<  Cut on tracks: minimum number of required ITS PID recpoints
  UShort_t              fRequireSDDrecPoints;   ///<  Cut on tracks: minimum number of required SDD recpoints
  UShort_t              fRequireSPDrecPoints;   ///<  Cut on tracks: minimum number of required SPD recpoints
  UShort_t              fRequireTPCrecPoints;   ///<  Cut on tracks: minimum number of required TPC recpoints
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
  Int_t                 fRequireMagneticField;  ///<  {0 : any magnetic field is fine, -1 : only negative magnetic field, 1 : only positive}
  Bool_t                fRequireVetoSPD;        ///<  Cut away all the tracks with at least 1 SPD cluster
  Float_t               fRequireMaxMomentum;    ///<  Cut in momentum for TPC only spectrum
  Float_t               fRequireTrackLength;    ///<  Cut on the track length
  Bool_t                fFixForLHC14a6;         ///<  Switch on/off the fix for the MC centrality distribution
  Bool_t                fNewCentralityFramework;///<  Use the new centrality framework

  Float_t               fRequireMinCentrality;  ///<  Max centrality
  Float_t               fRequireMaxCentrality;  ///<  Min centrality
  Bool_t                fEnableFlattening;      ///<  Switch on/off the flattening
  ULong_t               fTriggerMask;           ///<  Mask of the accepted triggers
  Bool_t               fEnableLogAxisInPerformancePlots; ///< Switch on/off logarithmic bins

  AliPID::EParticleType fParticle;              ///<  Particle specie
  TArrayF               fCentBins;              ///<  Centrality bins
  TArrayF               fDCABins;               ///<  DCA bins
  TArrayF               fPtBins;                ///<  Transverse momentum bins
  TArrayF               fCustomTPCpid;          ///<  Custom parametrisation of the Bethe-Bloch
  TArrayF               fFlatteningProbs;       ///<  Flattening probabilities
  TArrayF               fPhiRegions[2];            ///<  Azimuthal angle region where the analysis is performed

  // Event related histograms
  TH1F                 *fCentrality;            //!<! Events centrality distribution
  TH1F                 *fFlattenedCentrality;   //!<! Events centrality distribution after the flattening
  TH1F                 *fCentralityClasses;     //!<! Events statistics per centrality classes

  // MC only histograms
  TH1F                 *fProduction;            //!<! *(MC only)* Total number of produced particles
  TH2F                 *fAITS_TPC;              //!<! *(MC only)* Tracks reconstructed in ITS-TPC acceptance (anti-matter)
  TH2F                 *fAITS_TPC_TOF;          //!<! *(MC only)* Tracks reconstructed in ITS-TPC-TOF acceptance (anti-matter)
  TH2F                 *fATotal;                //!<! *(MC only)* Particles in acceptance (anti-matter)
  TH2F                 *fAPtCorrection;         //!<! *(MC only)* \f$p_{T}^{rec}-p_{T}^{MC}\f$ as a function of \f$p_{T}^{rec}\f$ for anti-matter
  TH2F                 *fMITS_TPC;              //!<! *(MC only)* Tracks reconstructed in ITS-TPC acceptance (matter)
  TH2F                 *fMITS_TPC_TOF;          //!<! *(MC only)* Tracks reconstructed in ITS-TPC-TOF acceptance (matter)
  TH2F                 *fMTotal;                //!<! *(MC only)* Particles in acceptance (matter)
  TH2F                 *fMPtCorrection;         //!<! *(MC only)* \f$p_{T}^{rec}-p_{T}^{MC}\f$ as a function of \f$p_{T}^{rec}\f$ for matter
  TH3F                 *fMDCAPrimaryTPC;        //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of primaries for ITS+TPC tracks
  TH3F                 *fMDCASecondaryTPC;      //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of secondaries for ITS+TPC tracks
  TH3F                 *fMDCAPrimaryTOF;        //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of primaries for ITS+TPC+TOF tracks
  TH3F                 *fMDCASecondaryTOF;      //!<! *(MC only)* \f$DCA_{xy}\f$ distribution of secondaries for ITS+TPC+TOF tracks

  // Data histograms
  TH3F                 *fATOFsignal;            //!<! *(Data only)* TOF signal for anti-matter
  TH3F                 *fATPCcounts;            //!<! *(Data only)* TPC counts for anti-matter
  TH3F                 *fATOFphiSignal;         //!<! *(Data only)* TOF signal for anti-matter as a function of \f$\phi\f$
  TH2F                 *fATPCphiCounts;         //!<! *(Data only)* TPC counts for anti-matter as a function of \f$\phi\f$
  TH2F                 *fATPCeLoss;             //!<! *(Data only)* TPC dE/dx for anti-matter
  TH3F                 *fMDCAxyTPC;             //!<! *(Data only)* \f$DCA_{xy}\f$ distribution for ITS+TPC tracks
  TH3F                 *fMDCAzTPC;              //!<! *(Data only)* \f$DCA_{z}\f$ distribution for ITS+TPC tracks
  TH3F                 *fMDCAxyTOF;             //!<! *(Data only)* \f$DCA_{xy}\f$ distribution for ITS+TPC+TOF tracks
  TH3F                 *fMDCAzTOF;              //!<! *(Data only)* \f$DCA_{z}\f$ distribution for ITS+TPC+TOF tracks
  TH3F                 *fMTOFsignal;            //!<! *(Data only)* TOF signal for matter
  TH3F                 *fMTPCcounts;            //!<! *(Data only)* TPC counts for matter
  TH3F                 *fMTOFphiSignal;         //!<! *(Data only)* TOF signal for matter as a function of \f$\phi\f$
  TH2F                 *fMTPCphiCounts;         //!<! *(Data only)* TPC counts for matter as a function of \f$\phi\f$
  TH2F                 *fMTPCeLoss;             //!<! *(Data only)* TPC dE/dx for matter

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskNucleiYield, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskNucleiYield__) */
