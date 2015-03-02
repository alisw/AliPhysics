//TODO: complete the documentation
/// \class AliAnalysisTaskNucleiYield
/// \brief This task fills histograms required to perform the analysis on the light nuclei yield.
///
/// The histograms filled in this tasks are used by the analysis macro
/// ## Monte Carlo
/// There are mainly three items studied here:
/// * The acceptance x efficiency: ..
///
/// \author Maximiliano Puccio <maximiliano.puccio@cern.ch>, University and INFN Torino
/// \date Feb 18, 2015

#ifndef __AliAnalysisTaskNucleiYield__
#define __AliAnalysisTaskNucleiYield__

#include "AliAnalysisTaskSE.h"
#include <Rtypes.h>
#include <TString.h>
#include <TArrayF.h>
#include <TList.h>
#include <AliPIDResponse.h>
class TH2F;
class TH3F;
class AliFlowTrackCuts;
class AliAODTrack;
class AliVVertex;
class AliPIDResponse;

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
  void SetRequireSPDrecPoints (int rec = 1) { fRequireSPDrecPoints = rec; }
  void SetEtaRange (float emin, float emax) { fRequireEtaMin = emin; fRequireEtaMax = emax; }
  void SetYRange (float ymin, float ymax) { fRequireYmin = ymin; fRequireYmax = ymax; }
  void SetRequireMaxChi2 (float maxChi2 = 4.f) { fRequireMaxChi2 = maxChi2; }
  void SetRequireMaxDCAxy (float maxDCA) { fRequireMaxDCAxy = maxDCA; }
  void SetRequireMaxDCAz (float maxDCA) { fRequireMaxDCAz = maxDCA; }
  void SetRequireTPCpidSigmas (float sig) { fRequireTPCpidSigmas = (sig > 0) ? sig : -sig; }
  void SetRequireITSpidSigmas (float sig) { fRequireITSpidSigmas = sig; }
  
  void SetCentBins(Int_t nbins, Float_t *bins);
  void SetDCABins(Int_t nbins, Float_t min, Float_t max);
  void SetPtBins(Int_t nbins, Float_t *bins);
  void SetCustomTPCpid(Float_t *par, Float_t sigma);
  void SetTOFBins(Int_t nbins, Float_t min, Float_t max);
  void SetDCAzBins(Int_t nbins, Float_t limit);
  
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
  
  TList                *fList;                  ///< Output list
  Int_t                 fPDG;                   ///< PDG code of the particle of interest
  Float_t               fPDGMass;               ///< PDG mass
  Float_t               fPDGMassOverZ;          ///< PDG mass over z
  Bool_t                fIsMC;                  ///< Switch between MC and data
  Bool_t                fFillOnlyEventHistos;   ///< Set treu to fill only event related histograms
  
  AliPIDResponse       *fPID;                   //!< PID response class
  Float_t               fMagField;              ///<  Magnetic field value for the current event
  AliVVertex           *fPrimaryVertex;         //!< Primary vertex of the current event
  
  TList                 fMmc;                   ///<
  TList                 fAmc;                   ///<
  
  Float_t               fDCAzLimit;             ///<
  Int_t                 fDCAzNbins;             ///<
  
  TArrayF               fPtCorrectionA;         ///<
  TArrayF               fPtCorrectionM;         ///<
  
  Float_t               fTOFlowBoundary;        ///<
  Float_t               fTOFhighBoundary;       ///<
  Int_t                 fTOFnBins;              ///<
  
  Bool_t                fRequireITSrefit;       ///< Cut on tracks: set true to require ITS refit
  Bool_t                fRequireTPCrefit;       ///< Cut on tracks: set true to require TPC refit
  Bool_t                fRequireNoKinks;        ///< Cut on tracks: set true to exclude tracks from kink vertices
  UShort_t              fRequireITSrecPoints;   ///< Cut on tracks: minimum number of required ITS recpoints
  UShort_t              fRequireITSsignal;      ///< Cut on tracks: minimum number of required ITS PID recpoints
  UShort_t              fRequireSPDrecPoints;   ///< Cut on tracks: minimum number of required SPD recpoints
  UShort_t              fRequireTPCrecPoints;   ///< Cut on tracks: minimum number of required TPC recpoints
  UShort_t              fRequireTPCsignal;      ///< Cut on tracks: minimum number of required TPC PID recpoints
  Float_t               fRequireEtaMin;         ///< Cut on tracks: minimum eta for the track
  Float_t               fRequireEtaMax;         ///< Cut on tracks: maximum eta for the track
  Float_t               fRequireYmin;           ///< Cut on tracks: mimimum y for the track (using PDG mass)
  Float_t               fRequireYmax;           ///< Cut on tracks: maximum y for the track (using PDG mass)
  Float_t               fRequireMaxChi2;        ///< Cut on tracks: maximum TPC \f$\chi^{2}/NDF\f$
  Float_t               fRequireMaxDCAxy;       ///< Cut on tracks: maximum \f$DCA_{xy}\f$ for the track
  Float_t               fRequireMaxDCAz;        ///< Cut on tracks: maximum \f$DCA_{z}\f$ for the track
  Float_t               fRequireTPCpidSigmas;   ///< Cut on TPC PID number of sigmas
  Float_t               fRequireITSpidSigmas;
  
  AliPID::EParticleType fParticle;              ///< Particle specie
  TArrayF               fCentBins;              ///< Centrality bins
  TArrayF               fDCABins;               ///< DCA bins
  TArrayF               fPtBins;                ///< Transverse momentum bins
  TArrayF               fCustomTPCpid;          ///< Custom parametrisation of the Bethe-Bloch
  
  // Event related histograms
  TH1F                 *fCentrality;            //!< Events centrality distribution
  TH1F                 *fFlattenedCentrality;   //!< Events centrality distribution after the flattening
  TH1F                 *fCentralityClasses;     //!< Events statistics per centrality classes
  
  // MC only histograms
  TH1F                 *fProduction;            //!< *MC only* Total number of produced particles
  TH2F                 *fAITS_TPC;              //!< *MC only* Tracks reconstructed in ITS-TPC acceptance
  TH2F                 *fAITS_TPC_TOF;          //!< *MC only* Tracks reconstructed in ITS-TPC-TOF acceptance
  TH2F                 *fATotal;                //!< *MC only*
  TH2F                 *fAPtCorrection;         //!< *MC only* \f$p_{T}^{rec}-p_{T}^{MC}\f$ as a function of \f$p_{T}^{rec}\f$ for anti-matter
  TH2F                 *fMITS_TPC;              //!< *MC only*
  TH2F                 *fMITS_TPC_TOF;          //!< *MC only*
  TH2F                 *fMTotal;                //!< *MC only*
  TH2F                 *fMPtCorrection;         //!< *MC only* \f$p_{T}^{rec}-p_{T}^{MC}\f$ as a function of \f$p_{T}^{rec}\f$ for matter
  TH3F                 *fMDCAPrimary;           //!< *MC only*
  TH3F                 *fMDCASecondary;         //!< *MC only*
  
  // Data histograms
  TH3F                 *fATOFsignal;            //!<
  TH3F                 *fATPCcounts;            //!<
  TH3F                 *fMDCAxy;                //!<
  TH3F                 *fMDCAz;                 //!<
  TH3F                 *fMTOFsignal;            //!<
  TH3F                 *fMTPCcounts;            //!<
  
  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskNucleiYield, 1);
  /// \endcond CLASSDEF
};


#endif /* defined(__AliAnalysisTaskNucleiYield__) */
