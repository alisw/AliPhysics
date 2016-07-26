#ifndef ALIANALYSISTASKHYPERTRITON3_H
#define ALIANALYSISTASKHYPERTRITON3_H


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
// AliAnalysisTaskHypertriton3 class
// analysis task for the study of the production of hypertriton
// which decays in 3 prongs: d+p+pi^-
// This task is optimized for ESDs.root
//
// Author:
// S. Trogolo, trogolo@to.infn.it
///////////////////////////////////////////////////////////////////////////

#include <TROOT.h>

#include "AliAnalysisTaskSE.h"
#include <TString.h>

class TChain;
class TH1F;
class TH2F;
class TList;
class TObjArray;
class TTree;

class AliAODVertex;
class AliESDEvent;
class AliESDtrackCuts;
class AliESDtrack;
class AliESDVertex;
class AliPID;
class AliPIDResponse;
class AliVertexerTracks;

class AliAnalysisTaskHypertriton3 : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskHypertriton3(TString taskname = "taskHypertriton");
  virtual ~AliAnalysisTaskHypertriton3();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t*);
  virtual void   Terminate(Option_t*);

  void SetReadMC(Bool_t flag = kTRUE) {fMC = flag;}
  void SetFillTree(Bool_t outTree = kFALSE) {fFillTree = outTree;}
  void SetTriggerConfig(UShort_t trigConf) {fTriggerConfig = trigConf;}
  void SetRequestRefit(bool itsR = kFALSE, bool itsRpion = kFALSE) {fRequestITSrefit = itsR; fRequestITSrefitPion = itsRpion;}
  void SetTOFpid(bool reqTOFpid = kFALSE) {fRequestTOFPid = reqTOFpid;}
  void SetRequestTPCSigmas(float tpcSgm) {fRequestTPCSigmas = tpcSgm;}
  void SetRequestTOFSigmas(float tofSgm) {fRequestTOFSigmas = tofSgm;}
  void SetChargeTriplet(bool sign_c = kTRUE, bool ls_c = kTRUE) {fMinvSignal = sign_c; fMinvLikeSign = ls_c;}
  void SetMotherType(bool matter = kTRUE, bool antimatter = kTRUE){fChooseMatter = matter; fChooseAntiMatter = antimatter;}
  void SetSideBand(Bool_t sband = kFALSE) {fSideBand = sband;}

  void SetDCAPionPrimaryVtx(double dcapionpv) {fDCAPiPVmin = dcapionpv;}
  void SetDCAzProtonPrimaryVtx(double dcaprotonpv) {fDCAzPPVmax = dcaprotonpv;}
  void SetDCAzDeuteronPrimaryVtx(double dcadeuteronpv) {fDCAzDPVmax = dcadeuteronpv;}

  void SetCosinePointingAngle(double mincp) {fCosPointingAngle = mincp;}
  void SetMaxDecayLength(double maxdl) {fMaxDecayLength = maxdl;}
  void SetMinDecayLength(double mindl) {fMinDecayLength = mindl;}
  void SetMinNormalizedDecayLength(double min_norm_dl) {fMinNormalizedDecL = min_norm_dl;}
  void SetMaxLifeTime(double max_ctau) {fMaxLifeTime = max_ctau;}
  void SetMinLifeTime(double min_ctau) {fMinLifeTime = min_ctau;}

  void SetRapidity(double rapid) {fRapidity = rapid;}

  void SetMaxPtMother(double maxpt) {fMaxPtMother = maxpt;}
  void SetMinPtMother(double minpt) {fMinPtMother = minpt;}
  void SetMaxPMotherCM(double maxp_cm){fMaxPMotherCM = maxp_cm;}

  void SetDCAPioDecayVtxXY(double maxpixy) {fDCAPiSVxymax = maxpixy;}
  void SetDCAPioDecayVtxZ(double maxpiz) {fDCAPiSVzmax = maxpiz;}
  void SetDCAProDecayVtx(double maxpro) {fDCAProSVmax = maxpro;}
  void SetDCADeuDecayVtx(double maxdeu) {fDCADeuSVmax = maxdeu;}

  void SetDCADeuteronProton(double maxdp) {fDCAdp = maxdp;}
  void SetDCAPionProton(double maxpip) {fDCApip = maxpip;}
  void SetDCADeuteronPion(double maxdpi) {fDCAdpi = maxdpi;}

  void SetAngleDeuteronProton(double ang_dp) {fAngledp = ang_dp;}
  void SetAngleDeuteronPion(double ang_dpi) {fAngledpi = ang_dpi;}

  void SetCentrPercentileLimits(double lowc, double highc) {fLowCentrality = lowc; fHighCentrality = highc;}

  Double_t GetDCAcut(Int_t part, Double_t dca) const;

  void SetConvertedAODVertices(AliESDVertex *ESDvtxp, AliESDVertex *ESDvtxs) const; //!<! method to set the value for the converted AOD vertices

 private:
  Bool_t HasTOF(AliESDtrack *trk, float &beta_tof);
  Bool_t PassPIDSelection(AliESDtrack *trk, int specie, Bool_t isTOFin); // specie according to AliPID enum: 2-pion, 4-proton, 5-deuteron
  void CombineThreeTracks(Bool_t isMatter, TArrayI arrD, TArrayI arrP, TArrayI arrPi, Bool_t cent0, Bool_t cent1);

  AliESDEvent        *fESDevent;                   ///< ESD event
  AliESDtrackCuts    *fESDtrackCuts;               ///< First set of ESD track cuts
  AliESDtrackCuts    *fESDtrackCutsV0;             ///< Track cuts applied only to \f$\pi^{-}\f$ and \f$\pi^{+}\f$ candidate
  AliESDVertex       *fPrimaryVertex;              //!<! Primary vertex of the current event
  AliPIDResponse     *fPIDResponse;                //!<! PID response class
  AliVertexerTracks  *fVertexer;                   //!<! Secondary vertex reconstructed with three candidate tracks

  AliAODVertex       *fVtx1;                       //!<! Primary vertex converted from ESD to AOD
  AliAODVertex       *fVtx2;                       //!<! Secondary vertex converted from ESD to AOD

  TObjArray          *fTrkArray;                   //!<! Array containing the three tracks candidated to the secondary vertex reconstruction

  //Variables
  Bool_t             fMC;                          ///< variables for MC selection
  Bool_t             fFillTree;                    ///< variables to fill the Tree
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
  Double_t           fDCAzPPVmax;                  ///< Cut on Max DCAz of p from primary vertex
  Double_t           fDCAzDPVmax;                  ///< Cut on Max DCAz of d from primary vertex
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
  Double_t           fDCAProSVmax;                 ///< Cut on proton DCA from reconstructed secondary vertex
  Double_t           fDCADeuSVmax;                 ///< Cut on deuteron DCA from reconstructed secondary vertex
  Double_t           fDCAdp;                       ///< Cut DCA deuteron-proton
  Double_t           fDCApip;                      ///< Cut DCA pion-proton
  Double_t           fDCAdpi;                      ///< Cut DCA deuteron-pion
  Double_t           fAngledp;                     ///< Cut on the angle between deuteron - proton
  Double_t           fAngledpi;                    ///< Cut on the angle between deuteron - pion
  Double_t           fMaxPMotherCM;                ///< Cut on max mother momentum in the CM
  Double_t           fLowCentrality;               ///< Cut on lower value of centrality class
  Double_t           fHighCentrality;              ///< Cut on high value of centrality class

  //Output list
  TList              *fOutput;                     ///< Output list

  //Histograms
  TH1F               *fHistCount;                           //!<! Total number of events
  TH1F               *fHistCentralityClass;                 //!<! Event statistic per centrality class
  TH1F               *fHistCentralityPercentile;            //!<! Event statistic per centrality percentile
  TH1F               *fHistTrigger;                         //!<! Event trigger statistic
  TH1F               *fHistMultiplicity;                    //!<! Event multiplicity
  TH1F               *fHistZPrimaryVtx;                     //!<! Primary vertex Z coordinate
  TH1F               *fHistXPrimaryVtx;                     //!<! Primary vertex X coordinate
  TH1F               *fHistYPrimaryVtx;                     //!<! Primary vertex Y coordinate
  TH1F               *fHistChi2perTPCcluster;               //!<! TPC \f$\chi^{2}/NDF\f$ tracks distribution
  TH1F               *fHistTrackFlagReco;                   //!<! Check of track flags: kTPCrefit, kITSout, kITSrefit


  //PID
  //--> TPC
  TH2F               *fHistTPCpid;                          //!<! TPC dE/dx vs \f$p_{TPC}\f$
  TH2F               *fHistTPCdeusignal;                    //!<! TPC PID: deuteron candidates
  TH2F               *fHistTPCprosignal;                    //!<! TPC PID: proton candidates
  TH2F               *fHistTPCpionsignal;                   //!<! TPC PID: \f$\pi^{-}\f$ candidates

  //--> TOF
  TH2F               *fHistTOFsignal;                       //!<! TOF \f$\beta\f$ vs \f$p_{TPC}\f$
  TH2F               *fHistTOFdeusignal;           //!<! TOF PID: deuteron candidates
  TH2F               *fHistTOFprosignal;           //!<! TOF PID: proton candidates
  //TH1F               *fHistTOFdeumass;                      //!<! TOF mass of deuteron identified with TPC
  //TH1F               *fHistTOFpromass;                      //!<! TOF mass of proton identified with TPC

  //Candidate combination
  // Data and MC histograms
  TH1F               *fHistpionTPCcls;                      //!<! TPC clusters distribution of candidate \f$\pi\f$
  //TH2F               *fHistCorrDCAdprimary;                 //!<! Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ deuteron-primary vertex
  //TH2F               *fHistCorrDCApprimary;                 //!<! Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ proton-primary vertex
  //TH2F               *fHistCorrDCApiprimary;                //!<! Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ pion-primary vertex
  TH1F               *fHistDCApiprimary;                    //!<! DCA pion-primary vertex distribution
  TH1F               *fHistDCAXYpiprimary;                  //!<! DCA_xy pion-primary vertex distribution
  TH1F               *fHistDCAZpiprimary;                   //!<! DCA_z pion-primary vertex distribution
  TH1F               *fHistDCApprimary;                     //!<! DCA proton-primary vertex distribution
  TH1F               *fHistDCAXYpprimary;                   //!<! DCA_xy proton-primary vertex distribution
  TH1F               *fHistDCAZpprimary;                    //!<! DCA_z proton-primary vertex distribution
  TH1F               *fHistDCAdprimary;                     //!<! DCA deuteron-primary vertex distribution
  TH1F               *fHistDCAXYdprimary;                   //!<! DCA_xy deuteron-primary vertex distribution
  TH1F               *fHistDCAZdprimary;                    //!<! DCA_z deuteron-primary vertex distribution
  TH1F               *fHistDCAdeupro;                       //!<! DCA deuteron-proton distribution
  TH1F               *fHistDCApiondeu;	                    //!<! DCA pion-deuteron distribution
  TH1F               *fHistDCApionpro;                      //!<! DCA pion-proton distribution
  TH2F               *fHistDCAdpdpi;
  TH2F               *fHistDCApdppi;
  TH2F               *fHistDCApidpip;
  TH1F               *fHistZDecayVtx;                       //!<! Reco secondary vertex Z coordinate
  TH1F               *fHistXDecayVtx;                       //!<! Reco secondary vertex X coordinate
  TH1F               *fHistYDecayVtx;                       //!<! Reco secondary vertex Y coordinate
  TH1F               *fHistDCApionvtx;                      //!<! \f$DCA\f$ candidate pion-secondary vertex
  TH1F               *fHistDCAXYpionvtx;                    //!<! \f$DCA_{xy}\f$ candidate pion-secondary vertex
  TH1F               *fHistDCAZpionvtx;                     //!<! \f$DCA_{z}\f$ candidate pion-secondary vertex
  TH1F               *fHistDCAprovtx;                       //!<! \f$DCA\f$ candidate proton-secondary vertex
  TH1F               *fHistDCAXYprovtx;                     //!<! \f$DCA_{xy}\f$ candidate proton-secondary vertex
  TH1F               *fHistDCAZprovtx;                      //!<! \f$DCA_{z}\f$ candidate proton-secondary vertex
  TH1F               *fHistDCAdeuvtx;                       //!<! \f$DCA\f$ candidate deuteron-secondary vertex
  TH1F               *fHistDCAXYdeuvtx;                     //!<! \f$DCA_{xy}\f$ candidate deuteron-secondary vertex
  TH1F               *fHistDCAZdeuvtx;                      //!<! \f$DCA_{z}\f$ candidate deuteron-secondary vertex
  TH1F               *fHistDecayLengthH3L;                  //!<! Decay length distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistNormalizedDecayL;                //!<! Normalized decay length distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistLifetime;                        //!<! c*tau distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistAngle_deu_pro;                   //!<! Angle between deuteron and proton vectors
  TH1F               *fHistAngle_deu_pion;                  //!<! Angle between deuteron and pion vectors
  TH1F               *fHistAngle_pro_pion;                  //!<! Angle between proton and pion vectors
  TH2F               *fHistAngleCorr_dp_dpi;                //!<! Correlation between angle_dp vs angle_dpi
  TH2F               *fHistAngleCorr_dp_ppi;                //!<! Correlation between angle_dp vs angle_ppi
  TH2F               *fHistAngleCorr_ppi_dpi;               //!<! Correlation between angle_ppi vs angle_dpi
  TH1F               *fHistHyperRapidity;                   //!<! Rapidity distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistCosPointingAngle;                //!<! Cosine of pointing angle distribution of candidate mother particle
  TH1F               *fHistPtPion;                          //!<! \f$\p^{T}\f$ distribution of candidate \f$\pi\f$
  TH1F               *fHistPtDeuteron;                      //!<! \f$\p^{T}\f$ distribution of candidate d
  TH1F               *fHistPtProton;                        //!<! \f$\p^{T}\f$ distribution of candidate p
  TH1F               *fHistPtHypertriton;                   //!<! hypertriton transverse momentum distribution
  TH1F               *fHistMassHypertriton;                 //!<! Invariant mass distribution of candidate reconstructed \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistMassAntiHypertriton;             //!<! Invariant mass distribution of candidate reconstructed anti-\f$H^{3}_{\Lambda}\f$
  TH1F               *fHistMassHypertriton_Cent;            //!<! Invariant mass distribution of candidate reconstructed \f$H^{3}_{\Lambda}\f$  centrality 0-10%
  TH1F               *fHistMassAntiHypertriton_Cent;        //!<! Invariant mass distribution of candidate reconstructed anti-\f$H^{3}_{\Lambda}\f$  centrality 0-10%
  TH1F               *fHistMassHypertriton_SemiCent;        //!<! Invariant mass distribution of candidate reconstructed \f$H^{3}_{\Lambda}\f$  centrality 10-50%
  TH1F               *fHistMassAntiHypertriton_SemiCent;    //!<! Invariant mass distribution of candidate reconstructed anti-\f$H^{3}_{\Lambda}\f$  centrality 10-50%

  TH1F               *fHistMassHypertriton_LS_Cent;         //!<! Invariant mass distribution of candidate reconstructed \f$H^{3}_{\Lambda}\f$  centrality 0-10%
  TH1F               *fHistMassAntiHypertriton_LS_Cent;     //!<! Invariant mass distribution of candidate reconstructed anti-\f$H^{3}_{\Lambda}\f$  centrality 0-10%
  TH1F               *fHistMassHypertriton_LS_SemiCent;     //!<! Invariant mass distribution of candidate reconstructed \f$H^{3}_{\Lambda}\f$  centrality 10-50%
  TH1F               *fHistMassAntiHypertriton_LS_SemiCent; //!<! Invariant mass distribution of candidate reconstructed anti-\f$H^{3}_{\Lambda}\f$  centrality 10-50%

  // MC only histograms
  TH1F               *fHistParticle;               //!<! *(MC only)* Reconstructed particles distribution per species through PDGCode cross-check
  TH1F               *fHistpionTPCclsMCt;          //!<! *(MC only)* TPC clusters distribution of candidate \f$\pi\f$ through PDGCode cross-check
  TH1F               *fHistpTpionMCt;              //!<! *(MC only)* \f$\p^{T}\f$ distribution of \f$\pi\f$ identified with PDGCode
  TH1F               *fHistpTproMCt;               //!<! *(MC only)* \f$\p^{T}\f$ distribution of proton identified with PDGCode
  TH1F               *fHistpTdeuMCt;               //!<! *(MC only)* \f$\p^{T}\f$ distribution of deuteron identified with PDGCode
  TH1F               *fHistMompionMCt;             //!<! *(MC only)* \f$\p\f$ distribution of \f$\pi\f$ identified with PDGCode
  TH1F               *fHistMomproMCt;              //!<! *(MC only)* \f$\p\f$ distribution of proton identified with PDGCode
  TH1F               *fHistMomdeuMCt;              //!<! *(MC only)* \f$\p\f$ distribution of deuteron identified with PDGCode
  TH2F               *fHistCorrDCAdprimaryMCt;     //!<! *(MC only)* Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ deuteron(PDGCode identified)-primary vertex
  TH2F               *fHistCorrDCApprimaryMCt;     //!<! *(MC only)* Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ proton(PDGCode identified)-primary vertex
  TH2F               *fHistCorrDCApiprimaryMCt;    //!<! *(MC only)* Correlation \f$DCA_{z}\f$ vs \f$DCA_{xy}\f$ pion(PDGCode identified)-primary vertex
  TH1F               *fHistDCApiprimaryMCt;        //!<! *(MC only)* DCA pion(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCApprimaryMCt;         //!<! *(MC only)* DCA proton(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCAdprimaryMCt;         //!<! *(MC only)* DCA deuteron(PDGCode identified)-primary vertex distribution
  TH1F               *fHistDCAdeuproMCt;           //!<! *(MC only)* DCA deuteron-proton distribution identified with PdgCode
  TH1F               *fHistDCApiondeuMCt;          //!<! *(MC only)* DCA pion-deuteron distribution identified with PdgCode
  TH1F               *fHistDCApionproMCt;          //!<! *(MC only)* DCA pion-proton distribution identified with PdgCode
  TH1F               *fHistZDecayVtxMCt;           //!<! *(MC only)* Reco secondary vertex Z coordinate - candidates daugthers particles identified with PDGCode
  TH1F               *fHistXDecayVtxMCt;           //!<! *(MC only)* Reco secondary vertex X coordinate - candidates daugthers particles identified with PDGCode
  TH1F               *fHistYDecayVtxMCt;           //!<! *(MC only)* Reco secondary vertex Y coordinate - candidates daugthers particles identified with PDGCode
  TH1F               *fHistDCAXYdeuvtxMCt;         //!<! *(MC only)* \f$DCA_{xy}\f$ candidate deuteron(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAZdeuvtxMCt;          //!<! *(MC only)* \f$DCA_{z}\f$ candidate deuteron(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAXYprovtxMCt;         //!<! *(MC only)* \f$DCA_{xy}\f$ candidate proton(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAZprovtxMCt;          //!<! *(MC only)* \f$DCA_{z}\f$ candidate proton(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAXYpionvtxMCt;        //!<! *(MC only)* \f$DCA_{xy}\f$ candidate pion(PDGCode identified)-secondary vertex
  TH1F               *fHistDCAZpionvtxMCt;         //!<! *(MC only)* \f$DCA_{z}\f$ candidate pion(PDGCode identified)-secondary vertex
  TH1F               *fHistNormalizedDecayL_MCt;   //!<! *(MC only)* normalized decay length of true hypertriton in 3 body decay
  TH1F               *fHistLifetime_MCt;           //!<! *(MC only)* c*tau distribution of true \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistAngle_deu_pro_MCt;      //!<! *(MC only)* Angle between deuteron and proton vectors
  TH1F               *fHistAngle_deu_pion_MCt;     //!<! *(MC only)* Angle between deuteron and pion vectors
  TH1F               *fHistAngle_pro_pion_MCt;     //!<! *(MC only)* Angle between proton and pion vectors

  TH2F               *fHistAngleCorr_dp_dpi_MCt;   //!<! *(MC only)* Correlation between angle_dp vs angle_dpi
  TH2F               *fHistAngleCorr_dp_ppi_MCt;   //!<! *(MC only)* Correlation between angle_dp vs angle_ppi
  TH2F               *fHistAngleCorr_ppi_dpi_MCt;  //!<! *(MC only)* Correlation between angle_ppi vs angle_dpi
  TH1F               *fHistHypertritonMomMCt;      //!<! *(MC only)* hypertriton momentum in the lab rest frame
  TH1F               *fHistHyperRapidityMCt;       //!<! *(MC only)* Rapidity distribution of candidate \f$H^{3}_{\Lambda}\f$
  TH1F               *fHistMassHypertritonMCt;     //!<! *(MC only)* Invariant mass distribution of reconstructed \f$H^{3}_{\Lambda}\f$ - daughters particles identified with PDGCode
  TH1F               *fHistMassAntiHypertritonMCt; //!<! *(MC only)* Invariant mass distribution of reconstructed \f$H^{3}_{\Lambda}\f$ - daughters particles identified with PDGCode

  //TTree
  TTree              *fTTree;                      //!<! Tree used for local tests and cross-check
  Float_t            fTCentralityPerc;
  Bool_t             fTMCtruth;
  // Deuteron
  Float_t            fTchi2NDFdeu;
  UShort_t           fTPCclsdeu;
  UShort_t           fTPCclsPIDdeu;
  Float_t            fTpTPCdeu;
  Float_t            fTpXdeu;
  Float_t            fTpYdeu;
  Float_t            fTpZdeu;
  Float_t            fTTPCnsigmadeu;
  Float_t            fTTOFmassdeu;
  Float_t            fTDCAXYdeuprvtx;
  Float_t            fTDCAZdeuprvtx;
  // Proton
  Float_t            fTchi2NDFpro;
  UShort_t           fTPCclspro;
  UShort_t           fTPCclsPIDpro;
  Float_t            fTpTPCpro;
  Float_t            fTpXpro;
  Float_t            fTpYpro;
  Float_t            fTpZpro;
  Float_t            fTTPCnsigmapro;
  Float_t            fTTOFmasspro;
  Float_t            fTDCAXYproprvtx;
  Float_t            fTDCAZproprvtx;
  // Pion
  Float_t            fTchi2NDFpion;
  UShort_t           fTPCclspion;
  UShort_t           fTPCclsPIDpion;
  Float_t            fTpTPCpion;
  Float_t            fTpXpion;
  Float_t            fTpYpion;
  Float_t            fTpZpion;
  Float_t            fTTPCnsigmapion;
  Float_t            fTDCAXYpioprvtx;
  Float_t            fTDCAZpioprvtx;
  // Triplet
  Float_t            fTDCAdp;
  Float_t            fTDCAdpi;
  Float_t            fTDCAppi;

  Float_t            fTDCAXYdvtx;
  Float_t            fTDCAZdvtx;
  Float_t            fTDCAXYpvtx;
  Float_t            fTDCAZpvtx;
  Float_t            fTDCAXYpivtx;
  Float_t            fTDCAZpivtx;

  Float_t            fTAngle_dp;
  Float_t            fTAngle_dpi;
  Float_t            fTAngle_ppi;

  Float_t            fTpdeu_gen_X;
  Float_t            fTpdeu_gen_Y;
  Float_t            fTpdeu_gen_Z;
  Float_t            fTppro_gen_X;
  Float_t            fTppro_gen_Y;
  Float_t            fTppro_gen_Z;
  Float_t            fTppio_gen_X;
  Float_t            fTppio_gen_Y;
  Float_t            fTppio_gen_Z;

  Int_t              fTpdgDeu;
  Int_t              fTpdgPro;
  Int_t              fTpdgPion;
  Int_t              fTmomidD;
  Int_t              fTmomidP;
  Int_t              fTmomidPi;
  Float_t            fTpdgmomD;
  Float_t            fTpdgmomP;
  Float_t            fTpdgmomPi;
  Int_t              fTuniqID_deu;
  Int_t              fTuniqID_pro;
  Int_t              fTuniqID_pion;


  Float_t            fTRapidity;
  Float_t            fTDecayLength;
  Float_t            fTDecayLengthError;
  Float_t            fTCosPA;
  Float_t            fTInvariantMass;


  AliAnalysisTaskHypertriton3(const AliAnalysisTaskHypertriton3&); // not implemented
  AliAnalysisTaskHypertriton3& operator=(const AliAnalysisTaskHypertriton3&); // not implemented

  ClassDef(AliAnalysisTaskHypertriton3, 3); // analysisclass

};

#endif
