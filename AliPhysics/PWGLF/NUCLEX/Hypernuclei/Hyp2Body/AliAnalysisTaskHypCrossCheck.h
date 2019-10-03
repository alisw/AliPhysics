#ifndef ALIANALYSISTASKHYPCROSSCHECK_H
#define ALIANALYSISTASKHYPCROSSCHECK_H


/**************************************************************************
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/



 ///////////////////////////////////////////////////////////////////////////
 // AliAnalysisTaskHypCrossCheck class
 // analysis task for the make some base checks for hypertriton analysis
 //
 // This task is optimized for ESDs.root
 //
 // Author:
 // S. Trogolo, trogolo@to.infn.it
 ///////////////////////////////////////////////////////////////////////////

#include <TROOT.h>

#include "AliAnalysisTaskSE.h"
#include <TString.h>
#include <TArrayF.h>
#include <AliPIDResponse.h>
#include <AliPID.h>


class TChain;
class TH1F;
class TH2F;
class TList;
class TObjArray;
class TTree;

class AliAODVertex;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDVertex;
class AliExternalTrackParam;
class AliPIDResponse;
class AliStack;
class AliVertexerTracks;

class AliAnalysisTaskHypCrossCheck : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskHypCrossCheck(TString taskname = "taskHypertriton");
  virtual ~AliAnalysisTaskHypCrossCheck();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t*);
  virtual void   Terminate(Option_t*);

  void SetReadMC(Bool_t flag = kTRUE) {fMC = flag;}
  void SetFillTree(Bool_t outTree = kFALSE, Bool_t genTree = kFALSE) {fFillTree = outTree; fFillTGen = genTree;}
  void SetRunPeriodSelection(Bool_t run1 = kFALSE, Bool_t run2 = kFALSE) {fRun1PbPb = run1; fRun2PbPb = run2;}
  void SetTriggerConfig(UShort_t trigConf) {fTriggerConfig = trigConf;}
  void SetEvtSpecie(UInt_t evspc = 4){fEvtSpecie = evspc;}
  void SetEmbedEvtSelection(Bool_t evtSel = kFALSE){fEvtEmbedSelection = evtSel;}
  void SetRequestRefit(bool itsR = kFALSE, bool itsRpion = kFALSE) {fRequestITSrefit = itsR; fRequestITSrefitPion = itsRpion;}
  void SetTOFpid(bool reqTOFpid = kFALSE) {fRequestTOFPid = reqTOFpid;}
  void SetCustomTPCpid (Float_t *par, Float_t sigma);
  void SetRequestTPCSigmas(float tpcSgm) {fRequestTPCSigmas = tpcSgm;}
  void SetRequestTOFSigmas(float tofSgm) {fRequestTOFSigmas = tofSgm;}
  void SetChargeTriplet(bool sign_c = kTRUE, bool ls_c = kTRUE) {fMinvSignal = sign_c; fMinvLikeSign = ls_c;}
  void SetMotherType(bool matter = kTRUE, bool antimatter = kTRUE){fChooseMatter = matter; fChooseAntiMatter = antimatter;}
  void SetSideBand(Bool_t sband = kFALSE) {fSideBand = sband;}

  void SetDCAPionPrimaryVtx(double dcapionpv) {fDCAPiPVmin = dcapionpv;}
  void SetDCAzHe3PrimaryVtx(double dcahelium3pv) {fDCAzHe3PVmax = dcahelium3pv;}

  void SetCosinePointingAngle(double mincp) {fCosPointingAngle = mincp;}
  void SetMaxDecayLength(double maxdl) {fMaxDecayLength = maxdl;}
  void SetMinDecayLength(double mindl) {fMinDecayLength = mindl;}
  void SetMinNormalizedDecayLength(double min_norm_dl) {fMinNormalizedDecL = min_norm_dl;}
  void SetMaxLifeTime(double max_ctau) {fMaxLifeTime = max_ctau;}
  void SetMinLifeTime(double min_ctau) {fMinLifeTime = min_ctau;}

  void SetRapidity(double rapid) {fRapidity = rapid;}

  void SetMaxPtMother(double maxpt) {fMaxPtMother = maxpt;}
  void SetMinPtMother(double minpt) {fMinPtMother = minpt;}

  void SetDCAPioDecayVtxXY(double maxpixy) {fDCAPiSVxymax = maxpixy;}
  void SetDCAPioDecayVtxZ(double maxpiz) {fDCAPiSVzmax = maxpiz;}
  void SetDCA3HeDecayVtx(double maxhe3) {fDCAHe3SVmax = maxhe3;}

  void SetDCA3HeliumPion(double max_he3pi) {fDCAhe3pi = max_he3pi;}

  void SetAngle3HeliumPion(double ang_he3pi) {fAnglehe3pi = ang_he3pi;}

  void SetCentrPercentileLimits(double lowc, double highc) {fLowCentrality = lowc; fHighCentrality = highc;}

  void SetConvertedAODVertices(AliESDVertex *ESDvtxp, AliESDVertex *ESDvtxs) const;  //!<! method to set the value for the converted AOD vertices


 private:
  Bool_t PassTriggerSelection(UInt_t PhysSelMask);
  Bool_t PassCentralitySelection();
  Bool_t HasTOF(AliESDtrack *trk, float &beta_tof);
  Float_t GetBetheAlephTPCnsigmas(AliESDtrack *t, AliPID::EParticleType part);
  Bool_t PassPIDSelection(AliESDtrack *trk, int specie, Bool_t isTOFin); // specie according to AliPID enum: 2-pion, 7-helium-
  Bool_t CheckPrimaryDistribution(AliESDtrack *trk, Int_t &pdgParticle);
  Bool_t DaughtersFromSameMother(int lpion, int lhelium3, Double_t &pThyp);
  void CheckGenerated();
  void CombineTwoTracks(Bool_t isMatter, TArrayI arrHe3, TArrayI arrPi, Int_t charge);
  //void SetParametersAtPrimaryVertex(AliESDtrack *trk, AliESDtrack *extpar) const;

  AliESDEvent        *fESDevent;                   ///< ESD event
  AliESDtrackCuts    *fESDtrackCuts;               ///< First set of ESD track cuts
  AliESDtrackCuts    *fESDtrackCutsV0;             ///< Track cuts applied only to \f$\pi^{-}\f$ and \f$\pi^{+}\f$ candidate
  AliESDVertex       *fPrimaryVertex;              //!<! Primary vertex of the current event
  AliPIDResponse     *fPIDResponse;                //!<! PID response class
  AliStack           *fStack;                      //<  Stack used for MC studies and cross-check
  AliVertexerTracks  *fVertexer;                   //!<! Secondary vertex reconstructed with two candidate tracks

  AliAODVertex       *fVtx1;                       //!<! Primary vertex converted from ESD to AOD
  AliAODVertex       *fVtx2;                       //!<! Secondary vertex converted from ESD to AOD

  TObjArray          *fTrkArray;                   //!<! Array containing the two tracks candidated to the secondary vertex reconstruction
  TArrayF             fCustomTPCpid;               ///<  Custom parametrisation of the Bethe-Bloch


  //Variables
  Bool_t             fMC;                          ///< variables for MC selection
  Bool_t             fFillTree;                    ///< variables to fill the Tree
  Bool_t             fFillTGen;                    ///< variables to fill the Generated TTree
  Bool_t             fRun1PbPb;
  Bool_t             fRun2PbPb;
  Bool_t             fEvtEmbedSelection;           ///< kTRUE: embed event selection in centrality estimation; kFALSE(default): event selection done by hand
  UInt_t             fEvtSpecie;
  Float_t            fCentrality;                  ///< Centrality class
  Float_t            fCentralityPercentile;        ///< Centrality percentile
  UShort_t           fTriggerConfig;               ///< select different trigger configuration
  Bool_t             fRequestITSrefit;             ///< flag for switch the ITSrefit request in the track cuts
  Bool_t             fRequestITSrefitPion;         ///< flag for switch the ITSrefit request only for candidate pion track cuts
  Float_t            fRequestTPCSigmas;            ///< number of sigmas for TPC pid
  Bool_t             fRequestTOFPid;               ///< switch on/off TOF pid
  Float_t            fRequestTOFSigmas;            ///< number of sigmas for TOF pid
  Bool_t             fChooseMatter;                ///< flag to switch on/off the study of the hypertriton
  Bool_t             fChooseAntiMatter;            ///< flag to switch on/off the study of the anti-hypertriton
  Bool_t             fMinvSignal;                  ///< flag for correct charge triplet - signal
  Bool_t             fMinvLikeSign;                ///< flag for like-sign charge triplet
  Bool_t             fSideBand;                    ///< select distributions in the side band region where only background

  //Cut variables
  Double_t           fDCAPiPVmin;                  ///< Cut on Min DCA of \f$\pi\f$ from primary vertex
  Double_t           fDCAzHe3PVmax;                ///< Cut on Max DCAz of He3 from primary vertex
  Double_t           fCosPointingAngle;            ///< Cut on Cosine of the pointing angle
  Double_t           fMaxDecayLength;              ///< Cut on maximum Decay length
  Double_t           fMinDecayLength;              ///< Cut on minimum Decay length
  Double_t           fMinNormalizedDecL;           ///< Cut on minimum normalized decay length
  Double_t           fMaxLifeTime;                 ///< Cut on maximum c*\f$\tau\f$
  Double_t           fMinLifeTime;                 ///< Cut on minimum c*\f$\tau\f$
  Double_t           fRapidity;                    ///< Cut on absolute value of mother rapidity y
  Double_t           fMaxPtMother;                 ///< Cut on max mother reconstructed \f$p_{T}\f$
  Double_t           fMinPtMother;                 ///< Cut on min mother reconstructed \f$p_{T}\f$
  Double_t           fDCAPiSVxymax;                ///< Cut on \f$\pi DCA_{xy}\f$ from reconstructed secondary vertex
  Double_t           fDCAPiSVzmax;                 ///< Cut on \f$\pi DCA_{z}\f$ from reconstructed secondary vertex
  Double_t           fDCAHe3SVmax;                 ///< Cut on {3}^He DCA from reconstructed secondary vertex
  Double_t           fDCAhe3pi;                     ///< Cut DCA {3}^He-pion
  Double_t           fAnglehe3pi;                   ///< Cut on the angle between {3}^He-pion
  Double_t           fLowCentrality;               ///< Cut on lower value of centrality class
  Double_t           fHighCentrality;              ///< Cut on high value of centrality class

  //Output list
  TList              *fOutput;                     ///< Output list

  //Histograms
  TH1F               *fHistCount;                  //!<! Total number of events
  TH1F               *fHistCentralityClass;        //!<! Event statistic per centrality class
  TH1F               *fHistCentralityPercentile;   //!<! Event statistic per centrality percentile
  TH1F               *fHistTrigger;                //!<! Event trigger statistic
  TH1F               *fHistMultiplicity;           //!<! Event multiplicity
  TH1F               *fHistZPrimaryVtx;            //!<! Primary vertex Z coordinate
  TH1F               *fHistXPrimaryVtx;            //!<! Primary vertex X coordinate
  TH1F               *fHistYPrimaryVtx;            //!<! Primary vertex Y coordinate
  TH1F               *fHistChi2perTPCcluster;      //!<! TPC \f$\chi^{2}/NDF\f$ tracks distribution
  TH1F               *fHistTrackFlagReco;          //!<! Check of track flags: kTPCrefit, kITSout, kITSrefit


  //PID
  //--> TPC
  TH2F               *fHistTPCpid;                 //!<! TPC dE/dx vs \f$p_{TPC}\f$
  TH2F               *fHistTPCHe3signal;           //!<! TPC PID: helium-3 candidates
  TH2F               *fHistTPCpionsignal;          //!<! TPC PID: \f$\pi^{-}\f$ candidates
  TH2F               *fHistNsigmaHe3;
  TH2F               *fHistNsigmaPion;

  //--> TOF
  TH2F               *fHistTOFsignal;              //!<! TOF \f$\beta\f$ vs \f$p_{TPC}\f$
  TH2F               *fHistTOFHe3signal;           //!<! TOF PID: helium-3 candidates

  //Candidate combination
  // Data and MC histograms
  TH1F               *fHistpionTPCcls;             //!<! TPC clusters distribution of candidate \f$\pi\f$
  //TH2F               *fHistCorrDCAHe3primary;        //!<! Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ helium-3-primary vertex
  //TH2F               *fHistCorrDCApiprimary;       //!<! Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ pion-primary vertex
  TH1F               *fHistDCApiprimary;                    //!<! DCA pion-primary vertex distribution
  TH1F               *fHistDCAXYpiprimary;                  //!<! DCA_xy pion-primary vertex distribution
  TH1F               *fHistDCAZpiprimary;                   //!<! DCA_z pion-primary vertex distribution
  TH1F               *fHistDCAHe3primary;                     //!<! DCA helium-3-primary vertex distribution
  TH1F               *fHistDCAXYHe3primary;                   //!<! DCA_xy helium-3-primary vertex distribution
  TH1F               *fHistDCAZHe3primary;                    //!<! DCA_z helium-3-primary vertex distribution
  TH1F               *fHistDCAhe3pion;                      //!<! DCA helium-3-pion distribution
  TH2F               *fHistDeltaPt_Pion;                    //!<!
  TH2F               *fHistDeltaPt_He3;                     //!<!
  TH1F               *fHistZDecayVtx;                       //!<! Reco secondary vertex Z coordinate
  TH1F               *fHistXDecayVtx;                       //!<! Reco secondary vertex X coordinate
  TH1F               *fHistYDecayVtx;                       //!<! Reco secondary vertex Y coordinate
  TH1F               *fHistDecayLengthH3L;                  //!<! Decay length distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistNormalizedDecayL;                //!<! Normalized decay length distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistDCAXYhe3vtx;                     //!<! \f$DCA_{xy}\f$ candidate helium-3-secondary vertex
  TH1F               *fHistDCAZhe3vtx;                      //!<! \f$DCA_{z}\f$ candidate helium-3-secondary vertex
  TH1F               *fHistDCAhe3vtx;
  TH1F               *fHistDCAXYpionvtx;                    //!<! \f$DCA_{xy}\f$ candidate pion-secondary vertex
  TH1F               *fHistDCAZpionvtx;                     //!<! \f$DCA_{z}\f$ candidate pion-secondary vertex
  TH1F               *fHistDCApionvtx;
  TH2F               *fHistDeltaPt_Hyper;                   //!<!
  TH1F               *fHistLifetime;                        //!<! c*tau distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistAngle_He3_pion;                   //!<! Angle between helium-3 and pion vectors
  TH1F               *fHistPtPion;                          //!<! pion transverse momentum distribution
  TH1F               *fHistPtHelium3;                        //!<! helium-3 transverse momentum distribution
  TH1F               *fHistPtHypertriton;                   //!<! hypertriton transverse momentum distribution
  TH1F               *fHistHyperRapidity;                   //!<! Rapidity distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistCosPointingAngle;                //!<! Cosine of pointing angle distribution of candidate mother particle
  TH1F               *fHistMassHypertriton;                 //!<! Invariant mass distribution of candidate reconstructed \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistMassAntiHypertriton;             //!<! Invariant mass distribution of candidate reconstructed anti-\f$H^{3}_{\Lambda}\f$
  TH2F               *fHistMassVsPt;                     //!<! Invariant mass distribution of candidate reconstructed \f$H^{3}_{\Lambda}\f$ and anti-\f$H^{3}_{\Lambda}\f$ vs pT
  TH2F               *fHistMassHyperVsPt;                //!<! Invariant mass distribution of candidate reconstructed \f$H^{3}_{\Lambda}\f$ vs pT
  TH2F               *fHistMassAntiHyperVsPt;            //!<! Invariant mass distribution of candidate reconstructed anti-\f$H^{3}_{\Lambda}\f$ vs pT

  // MC only histograms
  // Primary distribution
  TH2F               *fHistTPCdeusignal_pdg;       //!<! *(MC only)* deuteron (pdgcode) TPC specific energy
  TH2F               *fHistTPCtrisignal_pdg;       //!<! *(MC only)* triton (pdgcode) TPC specific energy
  TH2F               *fHistTPCHe4signal_pdg;       //!<! *(MC only)* helium-3 (pdgcode) TPC specific energy
  TH2F               *fHistTPCHe3signal_pdg;       //!<! *(MC only)* alpha (pdgcode) TPC specific energy
  TH2F               *fHistTPCHe3signal_3LH;
  TH2F               *fHistTPCpionsignal_3LH;
  TH2F               *fHistNsigmaHe3_3LH;
  TH2F               *fHistNsigmaPion_3LH;
  TH1F               *fHistParticle;               //!<! *(MC only)* Reconstructed particles distribution per species through PDGCode cross-check
  TH1F               *fHistParticle_Mass;
  TH1F               *fHistpionTPCclsMCt;          //!<! *(MC only)* TPC clusters distribution of candidate \f$\pi\f$ through PDGCode cross-check
  TH1F               *fHisthelium3TPCclsMCt;       //!<! *(MC only)* TPC clusters distribution of candidate p through PDGCode cross-check
  TH1F               *fHistpTpionMCt;              //!<! *(MC only)* \f$\p^{T}\f$ distribution of \f$\pi\f$ identified with PDGCode
  TH1F               *fHistpThe3MCt;               //!<! *(MC only)* \f$\p^{T}\f$ distribution of helium-3 identified with PDGCode
  TH1F               *fHistMompionMCt;             //!<! *(MC only)* \f$\p\f$ distribution of \f$\pi\f$ identified with PDGCode
  TH1F               *fHistMomHe3MCt;              //!<! *(MC only)* \f$\p\f$ distribution of helium-3 identified with PDGCode
  TH2F               *fHistCorrDCAHe3primaryMCt;     //!<! *(MC only)* Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ helium-3(PDGCode identified)-primary vertex
  TH2F               *fHistCorrDCApiprimaryMCt;    //!<! *(MC only)* Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ pion(PDGCode identified)-primary vertex
  TH1F               *fHistDCApiprimaryMCt;        //!<! *(MC only)* DCA pion(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCAXYpiprimaryMCt;      //!<! *(MC only)* DCA_xy pion(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCAZpiprimaryMCt;       //!<! *(MC only)* DCA_z pion(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCAHe3primaryMCt;         //!<! *(MC only)* DCA helium-3(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCAXYHe3primaryMCt;       //!<! *(MC only)* DCA_xy helium-3(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCAZHe3primaryMCt;        //!<! *(MC only)* DCA_z helium-3(PDGCode identified)-primary vertex distribution
  // Combining Tracks distributions
  TH1F               *fHistDCAhe3pionMCt;          //!<! *(MC only)* DCA helium-3-pion distribution identified with PdgCode
  TH2F               *fHistDeltaPt_PionMCt;        //!<! *(MC only)*
  TH2F               *fHistDeltaPt_He3MCt;         //!<! *(MC only)*
  TH1F               *fHistZDecayVtxMCt;           //!<! *(MC only)* Reco secondary vertex Z coordinate - candidates daugthers particles identified with PDGCode
  TH1F               *fHistXDecayVtxMCt;           //!<! *(MC only)* Reco secondary vertex X coordinate - candidates daugthers particles identified with PDGCode
  TH1F               *fHistYDecayVtxMCt;           //!<! *(MC only)* Reco secondary vertex Y coordinate - candidates daugthers particles identified with PDGCode
  TH1F               *fHistDecayLengthH3L_MCt;     //!<! *(MC only)* decay length of true hypertriton in 2 body decay
  TH1F               *fHistNormalizedDecayL_MCt;   //!<! *(MC only)* normalized decay length of true hypertriton in 2 body decay
  TH1F               *fHistDCAXYhe3vtxMCt;         //!<! *(MC only)* \f$DCA_{xy}\f$ candidate helium-3(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAZhe3vtxMCt;          //!<! *(MC only)* \f$DCA_{z}\f$ candidate helium-3(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAhe3vtxMCt;
  TH1F               *fHistDCAXYpionvtxMCt;        //!<! *(MC only)* \f$DCA_{xy}\f$ candidate pion(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAZpionvtxMCt;         //!<! *(MC only)* \f$DCA_{z}\f$ candidate pion(PDGCode identified)-secondary vertex
  TH1F               *fHistDCApionvtxMCt;
  TH2F               *fHistDeltaPt_HyperMCt;       //!<!
  TH1F               *fHistLifetime_MCt;           //!<! *(MC only)* c*tau distribution of true \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistAngle_He3_pion_MCt;     //!<! *(MC only)* Angle between helium-3 and pion vectors
  TH1F               *fHistHypertritonMomMCt;      //!<! *(MC only)* hypertriton momentum in the lab rest frame
  TH1F               *fHistPtPionMCt;              //!<! *(MC only)* pion transverse momentum distribution
  TH1F               *fHistPtHelium3MCt;            //!<! *(MC only)* helium-3 transverse momentum distribution
  TH1F               *fHistPtHypertritonMCt;
  TH1F               *fHistHyperRapidityMCt;       //!<! *(MC only)* Rapidity distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistCosPointingAngleMCt;    //!<! *(MC only)* Cosine of pointing angle distribution of candidate mother particle
  TH1F               *fHistMassHypertritonMCt;     //!<! *(MC only)* Invariant mass distribution of reconstructed \f$H^{3}_{\Lambda}\f$ - daughters particles identified with PDGCode
  TH1F               *fHistMassAntiHypertritonMCt; //!<! *(MC only)* Invariant mass distribution of reconstructed \f$H^{3}_{\Lambda}\f$ - daughters particles identified with PDGCode
  TH2F               *fHistMassVsPtMCt;            //!<! *(MC only)* Invariant mass distribution of reconstructed \f$H^{3}_{\Lambda}\f$ and anti-\f$H^{3}_{\Lambda}\f$ vs pT - daughters particles identified with PDGCode
  TH2F               *fHistMassHyperVsPtMCt;       //!<! *(MC only)* Invariant mass distribution of reconstructed \f$H^{3}_{\Lambda}\f$ vs pT - daughters particles identified with PDGCode
  TH2F               *fHistMassAntiHyperVsPtMCt;   //!<! *(MC only)* Invariant mass distribution of reconstructed anti-\f$H^{3}_{\Lambda}\f$ vs pT- daughters particles identified with PDGCode


  TH1F               *fHistHypertritonMomGen;
  TH1F               *fHistHypertritonMomGen_3Body;
  TH1F               *fHistAntiHypertritonMomGen_3Body;
  TH1F               *fHistHypertritonMomGen_isPrimary_3Body;
  TH1F               *fHistHypertritonMomGen_2Body;
  TH1F               *fHistAntiHypertritonMomGen_2Body;
  TH1F               *fHistHypertritonMomGen_isPrimary_2Body;

  TH1F               *fHistHypertritonYGen;
  TH1F               *fHistHypertritonYGen_3Body;
  TH1F               *fHistAntiHypertritonYGen_3Body;
  TH1F               *fHistHypertritonYGen_isPrimary_3Body;
  TH1F               *fHistHypertritonYGen_2Body;
  TH1F               *fHistAntiHypertritonYGen_2Body;
  TH1F               *fHistHypertritonYGen_isPrimary_2Body;

  //TTree
  TTree              *fTTree;                      //!<! Tree used for local tests and cross-check
  Float_t            fTCentralityPerc;
  Bool_t             fTMCtruth;
  // Helium-3
  Float_t            fTchi2NDFhe3;
  UShort_t           fTMhypoTrkhe3;
  UShort_t           fTPCclsPIDhe3;
  Bool_t             fITSrefithe3;
  Float_t            fTpTPChe3;
  Float_t            fTpXhe3;
  Float_t            fTpYhe3;
  Float_t            fTpZhe3;
  Float_t            fTTPCnsigmahe3;
  Float_t            fTDCAXYhe3prvtx;
  Float_t            fTDCAZhe3prvtx;
  // Pion
  Float_t            fTchi2NDFpion;
  UShort_t           fTMhypoTrkpion;
  UShort_t           fTPCclsPIDpion;
  Bool_t             fITSrefitpion;
  Float_t            fTpTPCpion;
  Float_t            fTpXpion;
  Float_t            fTpYpion;
  Float_t            fTpZpion;
  Float_t            fTTPCnsigmapion;
  Float_t            fTDCAXYpioprvtx;
  Float_t            fTDCAZpioprvtx;
  // Triplet
  Float_t            fTDCAhe3pi;

  Float_t            fTDCAXYhevtx;
  Float_t            fTDCAZhevtx;
  Float_t            fTDCAXYpivtx;
  Float_t            fTDCAZpivtx;

  Float_t            fTAngle_he3pi;

  Float_t            fTphe3_gen_X;
  Float_t            fTphe3_gen_Y;
  Float_t            fTphe3_gen_Z;
  Float_t            fTppio_gen_X;
  Float_t            fTppio_gen_Y;
  Float_t            fTppio_gen_Z;

  Int_t              fTpdgHe3;
  Int_t              fTpdgPion;
  Int_t              fTmomidHe3;
  Int_t              fTmomidPi;
  Float_t            fTpdgmomHe3;
  Float_t            fTpdgmomPi;
  Int_t              fTuniqID_he3;
  Int_t              fTuniqID_pion;

  Float_t            fTRapidity;
  Float_t            fTDecayLength;
  Float_t            fTDecayLengthError;
  Float_t            fTCosPA;
  Float_t            fTInvariantMass;

  //TTree
  TTree              *fTGen;                      //!<! Tree used for local test on generated particles
  Float_t            fTCentrality_gen;            //!<! Centrality of the generated event
  Float_t            fTHyper_px;                  //!<! x component of hypertriton momentum
  Float_t            fTHyper_py;                  //!<! y component of hypertriton momentum
  Float_t            fTHyper_pz;                  //!<! z component of hypertriton momentum
  Float_t            fTHyper_eta;                 //!<! eta of hypertriton
  Float_t            fTHyper_rapidity;            //!<! rapidity of hypertriton
  Float_t            fTHyper_pv_x;                //!<! x of production vertex of hypertriton
  Float_t            fTHyper_pv_y;                //!<! y of production vertex of hypertriton
  Float_t            fTHyper_pv_z;                //!<! z of production vertex of hypertriton
  Float_t            fTHyper_charge;              //!<! charge sign of hypertriton
  Float_t            fTHe3_px;                    //!<! x component of helium-3 momentum
  Float_t            fTHe3_py;                    //!<! y component of helium-3 momentum
  Float_t            fTHe3_pz;                    //!<! z component of helium-3 momentum
  Float_t            fTHe3_eta;                   //!<! eta of helium-3
  Float_t            fTHe3_rapidity;              //!<! rapidity of helium-3
  Float_t            fTHe3_pv_x;                  //!<! x of production vertex of helium-3
  Float_t            fTHe3_pv_y;                  //!<! y of production vertex of helium-3
  Float_t            fTHe3_pv_z;                  //!<! z of production vertex of helium-3
  Float_t            fTHe3_charge;                //!<! charge sign of helium-3
  Float_t            fTPion_px;                   //!<! x component of pion momentum
  Float_t            fTPion_py;                   //!<! y component of pion momentum
  Float_t            fTPion_pz;                   //!<! z component of pion momentum
  Float_t            fTPion_eta;                  //!<! eta of pion
  Float_t            fTPion_rapidity;             //!<! rapidity of pion
  Float_t            fTPion_charge;               //!<! charge sign of pion




  AliAnalysisTaskHypCrossCheck(const AliAnalysisTaskHypCrossCheck&); // not implemented
  AliAnalysisTaskHypCrossCheck& operator=(const AliAnalysisTaskHypCrossCheck&); // not implemented

  ClassDef(AliAnalysisTaskHypCrossCheck, 2); // analysisclass

};

#endif
