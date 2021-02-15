#ifndef ALIANALYSISTASKGAMMAHADRON_H
#define ALIANALYSISTASKGAMMAHADRON_H

// $Id$

#include "AliAnalysisTaskEmcal.h"
#include "AliEventPoolManager.h"
#include "AliAODMCHeader.h"
#include "AliEventCuts.h"
#include "AliFiducialCut.h"
#include "AliEMCALRecoUtils.h"
#include <THn.h>
#include <THnSparse.h>
#include <TProfile2D.h>

#include "AliQnCorrectionsManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"

class TH1;
class TH2;
class TH3;
class AliVVZERO;
class AliEvtPoolManager;
class AliAODMCHeader;

class AliQnCorrectionsManager;

using std::vector;

class AliAnalysisTaskGammaHadron : public AliAnalysisTaskEmcal {
public:
	AliAnalysisTaskGammaHadron();
	AliAnalysisTaskGammaHadron(Int_t InputGammaOrPi0,Int_t InputSameEventAnalysis,Bool_t InputMCorData);
	virtual ~AliAnalysisTaskGammaHadron();

	static AliAnalysisTaskGammaHadron* AddTaskGammaHadron(
			Int_t       InputGammaOrPi0 = 0,
			Int_t       InputSeMe       = 0,
			Bool_t      InputMCorData   = 0,
			UInt_t      evtTriggerType  = AliVEvent::kEMCEGA,
			UInt_t      evtMixingType   = AliVEvent::kAnyINT,
			Bool_t      isRun2          = 1,
			Double_t    trackptcut      = 0.15,
			Double_t    clusEcut        = 0.30,
			Bool_t      SavePool        = 0,
			const char *trackName       = "usedefault",
			const char *clusName        = "usedefault",
			const char *taskname        = "AliAnalysisTask",
			const char *suffix          = ""
	);

  //..setters for the analysis
  void                        SetDebug(Int_t input)                                 { fDebug           = input  ; }
  void                        SetCorrectEff(Bool_t input)                           { fCorrectEff      = input  ; }
  void                        SetEventWeightChoice(Int_t input)                     { fEventWeightChoice = input; }
  void                        SetSavePool(Bool_t input)                             { fSavePool        = input  ; }
  void                        SetSaveTriggerPool(Bool_t input)                      { fSaveTriggerPool = input  ; }
  void                        SetDownScaleMixTrigger(Float_t input)                 { fDownScaleMT     = input  ; }
  void                        SetPoolTrackDepth(Int_t input)                        { fTrackDepth      = input  ; }
  void                        SetPlotMore(Int_t input)                              { fPlotQA          = input  ; }
  void                        SetEvtTriggerType(UInt_t input)                       { fTriggerType     = input  ; }
  void                        SetPi0MassSelection(Int_t input);
  void                        SetTriggerPtCut(Double_t input)                       { fTriggerPtCut    = input  ; }
  void                        SetSubDetector(Int_t input)                           { fSubDetector     = input  ; }
  void                        SetEvtMixType(UInt_t input)                           { fMixingEventType = input  ; }
  void                        SetVetoTrigger(UInt_t input)                          { fVetoTrigger     = input  ; }
  void                        SetClEnergyMin(Double_t input)                        { fClEnergyMin     = input  ; }
  void                        SetOpeningAngleCut(Double_t input)                    { fOpeningAngleCut = input  ; }
  void                        SetNLM(Int_t input)                                   { fMaxNLM          = input  ; }
  void                        SetM02(Double_t inputMin,Double_t inputMax)           { fClShapeMin = inputMin; fClShapeMax = inputMax;}
  void                        SetClusterEnergyType(Int_t input)                     { fClusEnergyType = input   ; }
  void                        SetRmvMatchedTrack(Bool_t input, Double_t dEta=-1, Double_t dPhi=-1) { fRmvMTrack  = input; fTrackMatchEta=dEta; fTrackMatchPhi=dPhi;} // dEta, dPhi = -1 or 0 will use pt parametrized cut
  void                        SetHadronicCorrection(Int_t input, Double_t dEta=-1, Double_t dPhi=-1) { // dEta, dPhi = -1 or 0 will use pt parametrized cut
  fHadCorr = input; if(fHadCorr == 0) return;
  fTrackMatchEta=dEta; fTrackMatchPhi=dPhi; fClusEnergyType = AliVCluster::kHadCorr;
  AliWarning("This configuration modifies the HadCorr energy of input clusters. Following tasks beware!"); }
  void                        SetHadronicCorrectionEnergy(Double_t input)           { fHadCorrConstant = input;}
  void                        SetEOverPLimits(Double_t inputMin, Double_t inputMax) { fTrackMatchEOverPLow = inputMin; fTrackMatchEOverPHigh = inputMax; }
  void                        SetUseManualEvtCuts(Bool_t input)                     { fUseManualEventCuts= input;}
  void                        SetDoRotBkg(Bool_t input)                             { fDoRotBkg          = input;}
  void                        SetDoClusMixing(Bool_t input)                         { fDoClusMixing      = input;}
  void                        SetDoPosSwapMixing(Int_t input)                       { fDoPosSwapMixing   = input;}
  void                        SetMakePSHists(Bool_t input)                          { bEnablePosSwapHists = input;}
  void                        SetPSCorrectionLogMode(Bool_t input)                  { bLogPSMod          = input;}
  void                        SetEnableClusPairRot(Bool_t input)                     { if (input) {fDoPosSwapMixing = 2; bEnableClusPairRot=input;} else bEnableClusPairRot = 0;}
  void                        SetClusterDepth(Int_t input)                          { fClusterDepth      = input;}
  void                        SetNRotBkgSamples(Int_t input)                        { fNRotBkgSamples    = input;}
  void                        SetUseParamMassSigma(Bool_t input)                    { fUseParamMassSigma = input;}
  void                        SetPi0NSigma(Float_t input)                           { fPi0NSigma         = input;}
  void                        SetPi0AsymCut(Float_t input)                          { fPi0AsymCut        = input;}
  void                        SetApplyPatchCandCut(Bool_t input)                    { fApplyPatchCandCut = input;}
  void                        SetSidebandChoice(Int_t input)                        { fSidebandChoice    = input;}

  void                        SetMCEmbedReweightMode(Int_t input)                   { fMCEmbedReweightMode = input;}
  void                        SetUseMCReactionPlane(Int_t input)                    { fUseMCReactionPlane  = input;}

  //..Functions for mixed event purposes
  void                        SetExternalEventPoolManager(AliEventPoolManager* mgr) {fPoolMgr = mgr;}
  AliEventPoolManager*        GetEventPoolManager()                                 {return fPoolMgr;}
  TF1*                        GetEffFunction(Int_t no,Int_t cent)                           ;
  //..Set which pools will be saved
  void                        AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPt, Double_t maxPt);
   private:
  AliEventCuts                fEventCuts;                   ///< event selection utility
  AliFiducialCut*             fFiducialCuts;                ///< fiducial cuts for the EMCal and DCal in terms of eta and phi
  AliEMCALRecoUtils*          fFiducialCellCut;             ///< fiducial cut for EMCal+DCal in terms of rows and collumns
  AliQnCorrectionsManager*    fFlowQnVectorMgr;             ///< object for accessing QnVectorCorrections corrected event plane info

   protected:

  void                        InitArrays()                                                 ;
  //..overwritten EMCal base class functions
  Bool_t                      Run()                             	                          ;
  void                        ExecOnce()         									      ;
  Bool_t                      IsEventSelected()											  ;
  void                        UserCreateOutputObjects()        		                      ;

  //..Functions for mixed event purposes
  void                        InitEventMixer(Int_t MixMode = 0); // 0: Init for Mixed Tracks.  1: Init for Mixed Triggers
  void                        InitClusMixer()											  ;
  Int_t                       CalculateEventHash(); // Calculate a hash for the event for classifying and avoiding autocorrelations
  TObjArray*                  CloneToCreateTObjArray(AliParticleContainer* tracks)          ;

  //..Function for event plane purposes
  void                        LoadQnCorrectedEventPlane();


  //..Correlate and fill
  Bool_t                      FillHistograms()                                              ;
  void                        FillTrackHistograms(AliParticleContainer* tracks);
  Int_t                       CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, Double_t Weight);
  Int_t                       CorrelatePi0AndTrack(AliParticleContainer* tracks,TObjArray* bgTracks,Bool_t SameMix, Double_t Weight);
  void                        FillPi0CandsHist(AliTLorentzVector CaloClusterVec,AliTLorentzVector CaloClusterVec2,AliTLorentzVector CaloClusterVecPi0,Double_t fMaxClusM02,Double_t Weight,Int_t isMixed, Int_t mcIndex1 = -1, Int_t mcIndex2 = -1, Int_t PosSwapStatus = 0); // Pos swap status 1 is for conserved energy pair, 2 is for conserved positions pair
  void                        FillTriggerHist(AliTLorentzVector ClusterVec, Int_t CorrMCStatus, Double_t Weight);
  void                        FillGhHistograms(Int_t identifier,AliTLorentzVector ClusterVec,AliVParticle* TrackVec, Int_t CorrMCStatus, Double_t Weight);
  void                        FillQAHistograms(Int_t identifier,AliClusterContainer* clusters,AliVCluster* caloCluster,AliVParticle* TrackVec, Double_t Weight=1);
  Bool_t                      AccClusterForAna(AliClusterContainer* clusters, AliVCluster* caloCluster);
  Bool_t                      AccClusPairForAna(AliVCluster* cluster1, AliVCluster * cluster2, TLorentzVector vecPi0);
  Bool_t                      QuickCheckAccClus(TLorentzVector ClusterVec);
  Bool_t                      DetermineGAPatchCand(AliTLorentzVector CaloClusterVec, AliTLorentzVector CaloClusterVec2);
  Int_t                       DetermineMatchedTrack(AliVCluster* caloCluster,Double_t &etadiff,Double_t & phidiff, Bool_t bApplyHadCorr = 0);

  //..Functions for MC purposes
  
  void                        FillMCPi0Hists(Int_t fMultiplicity);
  Int_t                       FindMCPartForClus(AliVCluster * caloCluster);
  Int_t                       FindMCRootPart(Int_t iMCIndex, Int_t * iMCTreeHeight);
  Int_t                       FindMCLowComAnc(Int_t iMCIndex1, Int_t iMCIndex2);
 // Int_t                       CheckAcceptanceStatus(AliTLorentzVector Vec); // 0 If out of ECal/DCal, 1 if in ECal, 2 if in DCal


  //..Delta phi does also exist in AliAnalysisTaskEmcal. It is overwritten here (ask Raymond)
  Double_t                    DeltaPhi(AliTLorentzVector ClusterVec,AliVParticle* TrackVec) ;
  Double_t                    DeltaPhi(AliTLorentzVector ClusterVec,Double_t phi_EVP)       ;
  Double_t                    GetTrackEff(Double_t pT, Double_t eta)                        ;
  void                        GetDistanceToSMBorder(AliVCluster* caloCluster,Int_t &etaCellDist,Int_t &phiCellDist);
  AliVCluster*                GetLeadingCluster(const char* opt, AliClusterContainer* clusters);

  Int_t                       fGammaOrPi0;               ///< This tells me whether the correlation and the filling of histograms is done for gamma or pi0 or pi0 SB
  Int_t                       fSEvMEv;                   ///< This option performs the analysis either for same event or for mixed event analysis
	Bool_t                      fSaveTriggerPool;          ///< Whether to save an event pool of accepted triggers
  Float_t                     fDownScaleMT;              ///< Downscale factor to restrict Mixed Trigger statistics/runtime. 1.0 -> no downscale. 0.5 -> cut stats in half.
  Int_t                       fSidebandChoice;           ///< This determines which sideband option is used
  Bool_t                      fDebug;			        ///< Can be set for debugging
  Bool_t                      fSavePool;                 ///< Defines whether to save output pools in a root file
  Int_t                       fPlotQA;                   ///< plot additional QA histograms
  Bool_t                      fUseManualEventCuts;       ///< Use manual cuts if automatic setup is not available for the period
  Bool_t                      fCorrectEff;               ///< Correct efficiency of associated tracks
  Bool_t                      fEventWeightChoice;        ///< 0 = no event reweighting, 1 = reweight by GA function

  static const Int_t kNMainCentBins = 4;                 ///< Centrality bins that the analysis is done in (and mass windows are defined in)

  //..Input histograms
  TF1 						 *funcpT_low[4];            ///< input functions for the efficiency correction
  TF1 						 *funcpT_high[4];
  TF1 						 *funcpEta_left[4];
  TF1 						 *funcpEta_right[4];

  //..Constants
  Double_t                    fRtoD;                     ///< conversion of rad to degree
  static const Int_t          kNIdentifier=3;            ///< number of different versions of the same histogram type, can later be used for centrality or mixed event eg.
 // static const Int_t          kNvertBins=20;             ///< vertex bins in which the ME are mixed
  static const Int_t          kNvertBins=10;             ///< vertex bins in which the ME are mixed
  static const Int_t          kNcentBins=8;              ///< centrality bins in which the ME are mixed
  static const Int_t          kNClusVertBins=7;             ///< vertex bins in which the clusters are mixed
  static const Int_t          kNEMCalMultBins=4;              ///< EMCal multiplicity bins in which the clusters are mixed
  static const Int_t          kUsedPi0TriggerPtBins = 5; ///< Number of Bins used for Pi0 Triggers in Correlation mode
  static const Int_t          kNoGammaBins=9;            ///< Bins in gamma pT
  static const Int_t          kNoZtBins=7;               ///< Bins in Zt
  static const Int_t          kNoXiBins=8;               ///< Bins in Xi
  static const Int_t          kNoHPtBins=8;               ///< Bins in hadron pT
  Double_t                    fArray_G_Bins[10];         ///< 10=kNoGammaBins+1
  Double_t                    fArray_ZT_Bins[8];         ///< 8=kNoZtBins+1
  Double_t                    fArray_XI_Bins[9];         ///< 9=kNoXiBins+1
  Double_t                    fArray_HPT_Bins[9];        ///< 9=kNoHPtBins+1
  Double_t                    fArrayNVertBins[21];       ///< 21=kNvertBins+1

  static const Bool_t         bEnableTrackPtAxis = 1;    ///< Whether to swap the xi axis with a track pT axis. Currently must be set here
  static const Bool_t         bEnableEventHashMixing = 1;///< Whether to split events up into 2 classes (odd and even) for event mixing to avoid autocorrelation

  //..cuts
	Int_t                       fSubDetector;              ///< Whether to use all clusters, ECal only, or DCal only
  Double_t                    fTriggerPtCut;             ///< Cut of 5 GeV/c on Trigger Pt
  Double_t                    fMaxPi0Pt;                 ///< Cut of maximum pt for pi0 candidates
  Double_t                    fClShapeMin;               ///< Minimum cluster shape
  Double_t                    fClShapeMax;               ///< Maximum cluster shape
  Double_t                    fClEnergyMin;              ///< Minimum cluster energy for pi0 Reconstruction
  Double_t                    fOpeningAngleCut;          ///< Minimum opening angle cut for pi0 candidates
  Int_t                       fMaxNLM;                   ///< Maximum number of local maxima
  Bool_t                      fRmvMTrack;                ///< Switch to enable removing clusters with a matched track
  Int_t                       fClusEnergyType;           ///< Index of the energy from clusters (nonLinCorr=0,hadCorr=1)
  Int_t                       fHadCorr;                  ///< Switch to enable subtraction of constant energy for matched tracks. 0 = none, 1 = enable subtraction
  Double_t                    fHadCorrConstant;          ///< Energy value that is subtracted from clusters with matched tracks
  Double_t                    fTrackMatchEta;            ///< eta range in which a track is called a match to a cluster
  Double_t                    fTrackMatchPhi;            ///< phi range in which a track is called a match to a cluster
  Double_t                    fTrackMatchEOverPLow;      ///< Minimum E_cluster/p_track to accept cluster track match 
  Double_t                    fTrackMatchEOverPHigh;     ///< Maximum E_cluster/p_track to accept cluster track match (-1 for no cut)
  //..Event pool variables
  TAxis                      *fMixBCent;                 ///< Number of centrality bins for the mixed event
  TAxis                      *fMixBZvtx;                 ///< Number of vertex bins for the mixed event
  TAxis                      *fMixBEMCalMult;            ///< TAxis for EMCAL Multiplicity bins for the mixed clusters
  TAxis                      *fMixBClusZvtx;             ///< TAxis for vertex bins for the mixed clusters
  AliEventPoolManager        *fPoolMgr;                  ///< event pool manager
  Int_t                       fTrackDepth;               ///<  #tracks to fill pool
  Double_t                    fTargetFraction;           ///<  fraction of track depth before pool declared ready
  Int_t                       fClusterDepth;             ///<  #clusters to fill cluster mixing pool
  Int_t                       fPoolSize;                 ///<  Maximum number of events
  vector<vector<Double_t> >   fEventPoolOutputList;      //!<! ???vector representing a list of pools (given by value range) that will be saved
  //..Event selection types
  UInt_t                      fTriggerType;              ///<  Event types that are used for the trigger (gamma or pi0)
	Int_t                       fPi0MassSelection;         ///<  Selection of (mass,sigma) set used for pi0 Mass windows.
  UInt_t                      fMixingEventType;          ///<  Event types that are used for the tracks in the mixed event
  UInt_t                      fCurrentEventTrigger;      //!<! Trigger of the current event
  UInt_t                      fVetoTrigger;              //!<! Trigger that is vetoed in Mixed Events to avoid bias.  Default is EMCAL Gamma Trigger

  Bool_t                      fApplyPatchCandCut;        ///< Add GA Trigger patch candidate status to Pi0Cand THnSparse

  //..Event Plane variables
  Double_t                    fQnCorrEventPlaneAngle;    //!<! Event plane angle corrected by the QnVector framework. Filled by LoadQnCorrectedEventPlane
  Double_t                    fQnCorrEventPlane3Angle;    //!<! Event plane(3rd order) angle corrected by the QnVector framework. Filled by LoadQnCorrectedEventPlane
  Double_t                    fQnCorrEventPlane4Angle;    //!<! Event plane angle (4th order) corrected by the QnVector framework. Filled by LoadQnCorrectedEventPlane

  //..MC stuff
  Bool_t                      fParticleLevel;            ///< Set particle level analysis
  Bool_t                      fIsMC;                     ///< Trigger, MC analysis
  Int_t                       fMCEmbedReweightMode;      ///< Whether to reweight embedded MC particles. 0 = none, 1 = reweight function, 2 = veto embed etas
  Int_t                       fUseMCReactionPlane;       ///< Whether to set the 2nd order event plane to the reaction plane in the MC Header.
  UInt_t                      fAODfilterBits[2];         ///< AOD track filter bit map
  AliAODMCHeader             *fMCHeader;                 //!<! AOD MC Header
  AliMCParticleContainer     *fMCParticles;              ///< Container for MC Particles
  vector<AliAODMCParticle *>  fMCPi0List;                ///< List of MC Pi0s within EMCAL/DCal acceptance
  Double_t                    fMCReactionPlaneAngle;     //!<! Reaction plane angle from MC Truth

  //..Other stuff
  static const int            kCLUS_BUF_SIZE = 3000;     ///< Should cover all clusters in event
  signed char                 fClusterAcceptanceStatus[kCLUS_BUF_SIZE];   ///< List of cluster acceptance statuses (post container cuts) (-1 for rejected, 0 for unexamined
                                                         //   1 for passed, and manual hadronic correction applied)
  TList                      *fEventCutList;           //!<! Output list for event cut histograms
  TList                      *fOutputListQA;           //!<! Output list
  Double_t                    fscaleEta[4];             ///<

  // QnVector info and Event Plane Resolution Profiles
  TH1F            *fEPAngleV0M;              //!<! EP Angle from V0M
  TH1F            *fEPAngleTPCA;             //!<! EP Angle from TPC A (eta > 0)
  TH1F            *fEPAngleTPCC;             //!<! EP Angle from TPC C (eta < 0)

  TH1F            *fEP3AngleV0M;              //!<! EP3 Angle from V0M
  TH1F            *fEP3AngleTPCA;             //!<! EP3 Angle from TPC A (eta > 0)
  TH1F            *fEP3AngleTPCC;             //!<! EP3 Angle from TPC C (eta < 0)

  TH1F            *fEP4AngleV0M;              //!<! EP4 Angle from V0M
  TH1F            *fEP4AngleTPCA;             //!<! EP4 Angle from TPC A (eta > 0)
  TH1F            *fEP4AngleTPCC;             //!<! EP4 Angle from TPC C (eta < 0)


  static const int kNumEPROrders = 6;        ///<  How many orders of EPR to measure
  // For 2nd Order Event Plane
  TProfile2D     **fEPR_CosD1;               //!<! Cos(N[V0 - TPCA]). N is the index
  TProfile2D     **fEPR_CosD2;               //!<! Cos(N[V0 - TPCC])
  TProfile2D     **fEPR_CosD3;               //!<! Cos(N[TPCA - TPC])
  // For 3rd Order Event Plane
  TProfile2D     **fEP3R_CosD1;               //!<! Cos(N[V0 - TPCA]). N is the index
  TProfile2D     **fEP3R_CosD2;               //!<! Cos(N[V0 - TPCC])
  TProfile2D     **fEP3R_CosD3;               //!<! Cos(N[TPCA - TPC])
  // For 4th Order Event Plane
  TProfile2D     **fEP4R_CosD1;               //!<! Cos(N[V0 - TPCA]). N is the index
  TProfile2D     **fEP4R_CosD2;               //!<! Cos(N[V0 - TPCC])
  TProfile2D     **fEP4R_CosD3;               //!<! Cos(N[TPCA - TPC])



  //..Histograms -

  // Hists for MC Efficiency
  TH3             *fHistMCPi0_PtEtaMult;     //!<! MC pi0 dists (pt,eta,mult)
  TH3             *fHistMCPi0_PtEtaEP;       //!<! MC pi0 dists (pt,eta,EP bin)
  TH2             *fEtaPhiMCPion;            //!<! MC pi0 dists (eta,phi)

  TH2             *fHistClusPairInvarMasspT; //!<! Tyler's histogram
  TH1             *fHistPi0;                 //!<! Tyler's histogram
  TH2             *fMAngle;                  //!<! Tyler's histogram
  TH2             *fPtAngle;                 //!<! Tyler's histogram
  TH1             *fMassPionRej;             //!<! Histogram of Mass vs Pt for rejected Pi0 Candidates
  // Event Plane information (reconstructed event plane)
  TH2             *fPtEPAnglePionAcc;        //!<! Histogram of delta Psi_{EP} of accepted pi0 (vs pt)
  TH3             *fPtEPAnglePionAccCent;    //!<! Histogram of delta Psi_{EP} of accepted pi0 (vs pt and cent)
  TH2             *fPtEPAngleMCPion;         //!<! Histogram of delta Psi_{EP} of MC truth pi0 (vs pt)
  TH2             *fPtEPAngleTrueRecMCPion;  //!<! Histogram of delta Psi_{EP} (MC true angle) of properly reconstructed pi0s (vs MC true pt)
  // Event Plane information (reconstructed event plane, 3rd order)
  TH2             *fPtEP3AnglePionAcc;        //!<! Histogram of delta Psi_{EP,3} of accepted pi0 (vs pt)
  TH3             *fPtEP3AnglePionAccCent;    //!<! Histogram of delta Psi_{EP,2} of accepted pi0 (vs pt and cent)
  TH2             *fPtEP3AngleMCPion;         //!<! Histogram of delta Psi_{EP,3} of MC truth pi0 (vs pt)
  TH2             *fPtEP3AngleTrueRecMCPion;  //!<! Histogram of delta Psi_{EP,3} (MC true angle) of properly reconstructed pi0s (vs MC true pt)
  // Event Plane information (reconstructed event plane, 4th order)
  TH2             *fPtEP4AnglePionAcc;        //!<! Histogram of delta Psi_{EP,4} of accepted pi0 (vs pt)
  TH3             *fPtEP4AnglePionAccCent;    //!<! Histogram of delta Psi_{EP,4} of accepted pi0 (vs pt and cent)
  TH2             *fPtEP4AngleMCPion;         //!<! Histogram of delta Psi_{EP,4} of MC truth pi0 (vs pt)
  TH2             *fPtEP4AngleTrueRecMCPion;  //!<! Histogram of delta Psi_{EP,4} (MC true angle) of properly reconstructed pi0s (vs MC true pt)

  TH3             *fHistTrackPsiEPPtCent;    //!<! Histogram of delta Psi_{EP} of accepted tracks (vs pt and centrality)
  TH3             *fHistTrackPsiEP3PtCent;    //!<! Histogram of delta Psi_{EP,3} of accepted tracks (vs pt and centrality)
  TH3             *fHistTrackPsiEP4PtCent;    //!<! Histogram of delta Psi_{EP,4} of accepted tracks (vs pt and centrality)

  // Event Plane information (MC information);
  TH1             *fMCReactionPlane;         //!<! Histogram of distribution of reaction plane angle (MC information)
  TH2             *fPtRPAnglePionAcc;        //!<! Histogram of delta Psi_{RP} of accepted pi0 (vs pt)
  TH2             *fPtRPAngleMCPion;         //!<! Histogram of delta Psi_{RP} of MC truth pi0 (vs pt)
  TH2             *fPtRPAngleTrueRecMCPion;  //!<! Histogram of delta Psi_{RP} (MC true angle) of properly reconstructed pi0s (vs MC true pt)
  TH3             *fHistTrackPsiRPPtCent;    //!<! Histogram of delta Psi_{RP} of accepted tracks (vs pt and centrality)


  TH2             *fEtaPhiPionAcc;           //!<! Histogram of eta,phi location of accepted pions
  TH2             *fMassPtPionAcc;           //!<! Histogram of Mass vs Pt for accepted Pi0 Candidates
  TH2             *fMassPtPionRej;           //!<! Histogram of Mass vs Pt for rejected Pi0 Candidates
  TH3             *fMassPtCentPionAcc;       //!<! Histogram of Mass vs Pt vs Cent for accepted Pi0 Candidates
  TH3             *fMassPtCentPionRej;       //!<! Histogram of Mass vs Pt vs Cent for rejected Pi0 Candidates
  // Track Cluster Matching
  TH2             *fMatchDeltaEtaTrackPt;     //!<! Histogram of Delta eta vs track pt for cluster-track matching
  TH2             *fMatchDeltaPhiTrackPt;     //!<! Histogram of Delta phi vs track pt for cluster-track matching
  TH2             *fMatchCondDeltaEtaTrackPt;     //!<! Histogram of Delta eta vs track pt for cluster-track matching (Requiring delta phi cut)
  TH2             *fMatchCondDeltaPhiTrackPt;     //!<! Histogram of Delta phi vs track pt for cluster-track matching (Requiring delta eta cut)
  TH2             *fClusterEnergyMatchedTracks;   //!<! Histogram counting the number of matched tracks vs the energy of matched clusters
  TH2             *fHistEOverPvE;            //!<! Histogram of E/p vs E_cluster for cluster-track pairs (geometrically matched)
  TH2             *fHistPOverEvE;            //!<! Histogram of p/E vs E_cluster for cluster-track pairs (geometrically matched)
  // Position Swap Corrections
  TH2             *fHistPSDistU;             //!<! Histogram of sqrt((1-cos(theta_A))(1-cos(theta_B))) for Pos Swap Method
  TH2             *fHistPSDistV;             //!<! Histogram of sqrt(E_A/E_B) for Pos Swap Method

  TRandom3        *fRand;                      //!<! Random number generator.  Initialzed by rot background
  TH1             *fClusEnergy;                //!<! Energy of clusters accepted for pi0 analysis

  TH2             *fAccClusEtaPhi;             //!<! Histogram of eta phi distribution of accepted clusters
  TH3             *fAccClusEtaPhiZvtx;         //!<! Histogram of eta phi distribution of accepted clusters binned in Zvtx bins
  Bool_t          bEnableClusPairRot;          //!<! Whether to enable the cluster pair rotation method

  Bool_t          fDoRotBkg;                   ///< Whether or not to calculate the rotational background
  Bool_t          fDoClusMixing;               ///< Whether or not to use event pools to calculate mixed cluster pairs.
  Int_t           fDoPosSwapMixing;            ///< Whether to use Position Swapping ME method.  (0: Not used, 1: Use 1 random cluster, 2: use all avail clusters.)
  Int_t           fNRotBkgSamples;             ///< How many samples to use in the rotational background
  THnSparseF      *fPi0Cands;                  //!<! Michael's THnSparse for pi0 Candidates

  TH1F           *fHistEventHash;            //!<! Histogram tracking the event hash for dividing data

  // Position Swap Correction Histograms
  Bool_t          bEnablePosSwapHists;  ///<  Whether to produce the following histograms for investigating the position swap method (very memory intensive, should have high cluster cut)
  Bool_t          bLogPSMod;            ///<  Whether to store the scaling factors in log form
  // 2D Maps: These store the initial mass,pt and final mass,pt
  THnSparseF      *fPSMassPtMap;                     //!<! Mass,Pt 2D mapping for same E pairs
  THnSparseF      *fESMassPtMap;                     //!<! Mass,Pt 2D mapping for same Pos pairs (Energy Swap)
  // Scaling Versions: These store the scaling of the mass and pt
  THnSparseF      *fUScaleMatrix;                   //!<! Mass and Pt modification scaling matrix for same E Pairs
  THnSparseF      *fVScaleMatrix;                   //!<! Mass and Pt modification scaling matrix for same Pos Pairs


  TH2							*fEMCalMultvZvtx;            //!<! Histogram investigating z-vertex, EMCal Multiplicity for mixed cluster pairs

	// Monte Carlo Histograms
  TH2             *fHistClusMCDE;              //!<! Difference between detector Clus E and MC E (inclusive)
//  TH2             *fHistClusMCDEDRMatch;       //!<! Difference between detector Clus E and MC E (exclusive)
  TH2             *fHistClusMCDPhiDEta;        //!<! 2D Angle Difference between det. Clus and MC Clus. (incl.)
//  TH2             *fHistClusMCDPhiDEtaMatchDE; //!<! 2D Angle Difference between det. Clus and MC Clus. (excl.)

  TH2             *fHistPi0MCDPt;              //!<! Difference between detector pi0 pt and MC pt.
  TH2             *fHistEtaMCDPt;              //!<! Difference between detector Eta pt and MC pt.
  TH2             *fHistPi0MCDPhiDEta;         //!<! 2D Angle Difference between det. pi0 and MC pi0.
  TH2             *fHistEtaMCDPhiDEta;         //!<! 2D Angle Difference between det. Eta and MC Eta.


  Bool_t          fUseParamMassSigma;                           ///< Whether to use parametrized or fixed mass,sigma
  Double_t        fPi0MassFixed[kNMainCentBins][kNoGammaBins];  ///< Fixed Mass Peak per pT bin, per Main Cent bin
  Double_t        fPi0SigmaFixed[kNMainCentBins][kNoGammaBins]; ///< Fixed Sigma per pT bin, per Main Cent bin
  Double_t        fPi0MassFitPars[5];                           ///< Fit Parameters for mass peak
  Double_t        fPi0SigmaFitPars[5];                          ///< Fit Parameters for sigma
  Float_t         fPi0NSigma;                                   ///< number of sigma to cut on peak
  Float_t         fPi0AsymCut;                                  ///< Asymmetry cut for Pi0 candidates

  TH2                      **fEffCorrectionCheck;      //!<! applied efficiency correction in eta-pT to cross check
  TH1 					    *fHistEvsPt;               //!<! E vs pT
  TH1 					   **fHistBinCheckPt;          //!<! plot Pt distribution for ideal binning
  TH1 					   **fHistBinCheckZt;          //!<! plot Zt distribution for ideal binning
  TH1 					   **fHistBinCheckXi;          //!<! plot Xi distribution for ideal binning
  TH1 					   **fHistBinCheckEvtPl;       //!<! plot Event Plane distribution for ideal binning
  TH1 					   **fHistBinCheckEvtPl2;      //!<! plot reduced Event Plane distribution for ideal binning
  TH2                      **fHistDEtaDPhiGammaQA;     //!<! Distribution of gammas in delta phi delta eta
  TH2                      **fHistDEtaDPhiTrackQA;     //!<! Distribution of tracks in delta phi delta eta
  TH2                      **fHistClusterTime;         //!<! Cluster time vs energy
  THnSparseF                *fCorrVsManyThings;        //!<! Thn sparse filled with delta phi, delta eta,Eg,zt,xi,vertex Z,centrality...
  THnSparseF                *fTriggerHist;             //!<! Thn sparse filled with Eg,vertex Z,eventPlane Bin,centrality,
  THnSparseF                *fClusterProp;             //!<! Thn sparse filled with cluster properties
  TH2                 	    *fHPoolReady;             //!<! Check how many Jobs start mixing
  //

   private:
  AliAnalysisTaskGammaHadron(const AliAnalysisTaskGammaHadron&);            // not implemented
  AliAnalysisTaskGammaHadron &operator=(const AliAnalysisTaskGammaHadron&); // not implemented

  ClassDef(AliAnalysisTaskGammaHadron, 12) // Class to analyse gamma hadron correlations
};
#endif
