/*************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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
//
// Select events according to gap conditions, analyze two track events in pp
// collisions at sqrt(s)=13TeV and can be used at lego-train
//
// Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>
//  continued by
//  Taesoo Kim <taesoo.kim@cern.ch>

#ifndef ALIANALYSISTASKCDPWA_H
#define ALIANALYSISTASKCDPWA_H

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#endif

#include "THnSparse.h" // forward declaration is not possible

class AliESDEvent;
class AliVTrack;
class AliPID;
class AliPIDResponse;
class AliPIDCombined;
class AliPhysicsSelection;
class AliESDtrackCuts;
class AliESDAD;
class AliESDVZERO;

class TH1I;
class TH1F;
class TH1D;
class TH2I;
class TH2F;
class TH2D;
class TList;
class TTree;
class TString;
class TParticle;
class THnSparse;
class TObjString;
class TObjArray;

class AliAnalysisTaskCDPWA : public AliAnalysisTaskSE
{
	public:
		AliAnalysisTaskCDPWA(const char* name);
		AliAnalysisTaskCDPWA();
		virtual ~AliAnalysisTaskCDPWA();

		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t *);
		virtual void Terminate(Option_t *);
		void SetIsRun2(Bool_t isRun2) {fIsRun2 = isRun2;}//Run2 mode
		void SetIsMC(Bool_t MC) {fIsMC = MC;}//MC mode
		void SetSaveGapEvents(Bool_t isTree) {fSavemode = isTree;}//kTRUE = Save gap-event, kFALSE = Save 2,4-track
		void SetCombinatoricsMode(Bool_t isComb) {fCombmode = isComb;}//kTRUE = save combinatorics
		void SetSaveGenParticle(Bool_t isSaveGen) {fIsSaveGen = isSaveGen;}//kTRUE = save proton, 2pion info
		void SetIsPythia8(Bool_t ispythia8) {fIsPythia8 = ispythia8;}

	private:
		enum {
			// Gap Information
			kBinDG,
			kBinGC, // single gap c side
			kBinGA, // single gap a side
			kBinNG, // no gap
			kBinAll
		};
		enum {
			// MC event/process types
			kBinND = 0, // used for data as well, "Non-Diffractive"
			kBinCD, // central diffractive
			kBinSD1, // single diffractive
			kBinSD2, // single diffractive
			kBinDD, // double diffractive
			kBinEL, // Elastic event
			kBinMCAll
		};
		enum {
			// Event selection procedure histogram
			kInput = 0, // Total event
			kMCCheck, // For MC only..
			kOnlineTrigger, //Pass MB_OR online trigger
			kOfflineCut, //Offline MB_OR trigger
			kVtxCut, //Vertex cut
			kPileUpCut, //Pile-up cut by secondary SPD vertex
			kClusterCut, //SPD Cluster cut
			kMBOR_V0, //V0A || V0C || SPD
			kMBAND_V0, //V0A && V0C
			kMBOR_AD, //ADA || ADC ||SPD
			kMBAND_AD, //ADA && ADC
			kMBOR_Global, //V0A || V0C || ADA || ADC || SPD
			kMBAND_Global, //(V0A || ADA) && (V0C || ADC)
			kNG_Data,
			kSGA_Data,
			kSGC_Data,
			kDGV0,//!V0A && !V0C
			kDGAD,//!ADA && !ADC
			kDGV0SPD, //!V0A && !V0C && SPD
			kDGADSPD, //!ADA && !ADC && SPD
			kDGV0ADSPD,
			kDGV0ADFMDSPD,
			kDGV0ADFMDZDCSPD,
			kDGV0FMDSPD,
			kDGV0FMDZDCSPD,
			k2Tracks,
			k4Tracks,
			k2Tracks_ITSSA,
			k4Tracks_ITSSA,
			kAll
		};
		enum {
			kTri_CINT5 = 0,
			kTri_CINT7,
			kTri_CINT10,
			kTri_CADAND,
			kTri_C0SMB,
			kTri_CDG6,
			kTri_CDG6_SPD2,
			kTri_CDG7_SPD2,
			kTriAll
		};

		AliAnalysisTaskCDPWA(const AliAnalysisTaskCDPWA  &p);
		AliAnalysisTaskCDPWA& operator=(const AliAnalysisTaskCDPWA  &p);


		// functions called by the UserExec(...), not to be called elsewhere!
		//-------------------------------------------------------------------
		Bool_t CheckInput();// check input file is OK or not
		Bool_t CheckOnlineTrigger(const AliESDEvent *ESDEvent);// check MBOR online trigger is fired or not
		void DoTimeV0AD(const AliESDEvent *ESDEvent, TH1D *v0ahist, TH1D *v0chist, TH1D *adahist, TH1D *adchist);
		void DoSPDCheck(const AliESDEvent *ESDEvent, TH1D *h1, TH1D *h2, TH1D *h3);
		void PostOutputs(); // cares about posting the output before exiting UserExec,
		// WARNING: PostOutputs should only be used directly in the UserExec!!
		static Bool_t SPDLoc2Glo(const Int_t id, const Double_t *loc, Double_t *glo);
		static Int_t CheckChipEta(const Int_t chipKey, const TString scut,
				const Double_t vtxPos[], TH2 *hitMapSPDinner,
				TH2 *hitMapSPDouter);
		Bool_t CutEvent(const AliESDEvent *ESDEvent, TH1* hpriVtxX, TH1* hpriVtxY, TH1* hpriVtxZ);
		static Int_t GetFastORmultiplicity(const AliESDEvent *ESDEvent);
		Bool_t CheckV0Hit(const AliESDEvent *ESDEvent, TH1D* fHitV0A, TH1D* fHitV0C);
		Int_t DetermineGap(const AliESDEvent *ESDEvent, const Bool_t wCent, const Int_t type); // determines the gap of all available detectors
		Int_t DoMCTruth(); // analyses the MCtruth for corrections, returns #Primaries
		void DoCombStudy(const AliESDAD *ad, const AliESDVZERO *vzero); 

		// analyzes track pairs in MC truth
		void DetermineMCprocessType(); // determines the MC process ID

		// analysis task status
		//---------------------
		Double_t fMaxVtxDst; // maximum distance of the track and SPD vertex

		// event information
		//------------------
		AliPIDResponse *fPIDResponse; //! PID Response object
		AliPIDCombined *fPIDCombined; //! PID Combined object
		AliESDEvent *fESDEvent; // esd event object
		AliPhysicsSelection *fPhysicsSelection; // physics selection object
		AliESDtrackCuts *fTrackCuts;
		AliESDtrackCuts *fTrackCuts_ITSSA;
		Double_t fVtxDst; // distance of the primary vertex from tracks and from SPD
		Double_t fVtxZ; // z-position of the primary vertex from tracks
		Int_t fMCprocessType; // MC process type, 0 for data
		Int_t fMCprocess; // detailed MC sub process information

		// information about the trackpair which is currently processed
		Int_t fRun; // number of the run which is about to be processed
		Int_t fPIDmode; // selects set of PID cuts, 0 for 3sigma standard cuts,
		Bool_t fIsRun2;
		Bool_t fSavemode;
		Bool_t fCombmode;
		Bool_t fIsMC;
		Bool_t fIsSaveGen;
		Bool_t fIsPythia8;
		Bool_t fIsPythia;
		Bool_t fIsPhojet;
		Bool_t fIsEPOS;
		// on central activity

		// Output objects-----------------------------------------------------
		TTree *fTree; //! V0 2pion
		TTree *fTree_Comb;//! Combinatorics
		TList *fList; //! List for histogram

		// Tree variables----------------------------------------------------
		Bool_t fSPDFired;//!
		Bool_t fV0Gap;//!
		Bool_t fADGap;//!
		Bool_t fFMDGap;//!
		Bool_t fZDCGap;//!
		Int_t fRunNumber;//!
		Int_t fPeriod;//!
		Bool_t fCheckTwoPion; //! Check that this event has 2 tracks
		Bool_t fCheckFourPion; //! Check that this event has 4 tracks
		Bool_t fCheckTwoPion_ITSSA;//! Check that this event has 2 ITSSA_track events
		Bool_t fCheckFourPion_ITSSA;//! Check that this event has 2 ITSSA_track events
		Int_t fMultiplicity;//!
		Double_t fVertex[3];//! Vertex information

		Double_t fTwoPionTrack[2][10]; //! Two track Momentum, Energy and Sign
		Double_t fTwoPionTPCSigma[2][10]; //! Two track TPC sigma for PID
		Double_t fTwoPionTOFSigma[2][10]; //! Two track TOF sigma for PID
		Double_t fTwoPionITSSigma[2][10]; //! Two track ITS sigma for PID
		UInt_t fTwoPionMask_TPC[2]; //! Two track mask for TPC
		UInt_t fTwoPionMask_TOF[2]; //! Two track mask for TOF
		UInt_t fTwoPionMask_ITS[2]; //! Two track mask for ITS
		UInt_t fTwoPionMask_TRD[2]; //! Two track mask for TRD
		UInt_t fTwoPionMask_tot[2]; //! Two track mask for TRD
		UInt_t fTwoPionDetMask_tot[2]; //! Two track mask for TRD
		Double_t fTwoPionBayesProb_TPC[2][5]; //!Bayesian probabilities for TPC
		Double_t fTwoPionBayesProb_TOF[2][5]; //!Bayesian probabilities for TOF
		Double_t fTwoPionBayesProb_ITS[2][5]; //!Bayesian probabilities for ITS
		Double_t fTwoPionBayesProb_TRD[2][5]; //!Bayesian probabilities for TRD
		Double_t fTwoPionBayesProb_tot[2][5]; //!Bayesian probabilities for tot
		Double_t fTwoPionTrack_ITSSA[2][10]; //! Two track(ITSSA) Momentum, Energy and Sign
		Double_t fTwoPionTPCProb[2][5];//!
		Double_t fTwoPionTOFProb[2][5];//!
		Double_t fTwoPionITSProb[2][5];//!

		Double_t fFourPionTrack[4][10]; //! Four track Momentum, Energy and Sign
		Double_t fFourPionTPCSigma[4][10]; //! Four track Momentum, Energy and Sign
		Double_t fFourPionTOFSigma[4][10]; //! Four track Momentum, Energy and Sign
		Double_t fFourPionITSSigma[4][10]; //! Four track Momentum, Energy and Sign
		UInt_t fFourPionMask_TPC[4]; //! Four track mask for TPC
		UInt_t fFourPionMask_TOF[4]; //! Four track mask for TOF
		UInt_t fFourPionMask_ITS[4]; //! Four track mask for ITS
		UInt_t fFourPionMask_TRD[4]; //! Four track mask for TRD
		UInt_t fFourPionMask_tot[4]; //! Four track mask for TRD
		UInt_t fFourPionDetMask_tot[4]; //! Four track mask for TRD
		Double_t fFourPionBayesProb_TPC[4][5]; //!Bayesian probabilities for TPC
		Double_t fFourPionBayesProb_TOF[4][5]; //!Bayesian probabilities for TOF
		Double_t fFourPionBayesProb_ITS[4][5]; //!Bayesian probabilities for ITS
		Double_t fFourPionBayesProb_TRD[4][5]; //!Bayesian probabilities for TRD
		Double_t fFourPionBayesProb_tot[4][5]; //!Bayesian probabilities for tot
		Double_t fFourPionTrack_ITSSA[4][10];//!
		Double_t fFourPionTPCProb[4][5];//!
		Double_t fFourPionTOFProb[4][5];//!
		Double_t fFourPionITSProb[4][5];//!

		Double_t fMCGenProtonTrack[4][5]; //! Info of generated Proton 1 after scattering
		Double_t fMCGenPionTrack[4][5]; //! Info of generated Pion 1 after scattering

		Double_t fADCharge[16];//!
		Int_t fADANmbBB;//!
		Int_t fADCNmbBB;//!

		// Tree variables for combinatorics----------------------------------
		Bool_t fComb_IsPassMBOR;//! Need for MC
		Bool_t fComb_IsPassVertex;//! Need for MC
		Bool_t fComb_IsPileUp;//!
		Bool_t fComb_IsPassClusterCut;//!
		Int_t fComb_SPDCluster;//!
		Int_t fComb_SPDTracklets;//!
		Bool_t fComb_DetHit[10];//!
		Double_t fComb_V0_Time_Mean[2];//!
		Double_t fComb_V0_Time[64];//!
		Double_t fComb_V0_ADC[64];//!
		Double_t fComb_AD_Time_Mean[2];//!
		Double_t fComb_AD_Time[16];//!
		Double_t fComb_AD_ADC[16];//!
		Int_t fComb_MC_EventProcess;//!
		TLorentzVector fComb_forwardP[2];//!
		TLorentzVector fComb_diffSystem[2];//!


		//Histograms----------------------------------------------------------
		TH1D *fHistEvent; //Histogram for number of event
		TH1D *fHistTrigger;
		TH1D *fHistEventProcesses; //Histogram for number of event(MC)
		TH1D *fHistPrimVtxX; //Histogram for position of primary vertex X
		TH1D *fHistPrimVtxY; //Histogram for position of primary vertex X
		TH1D *fHistPrimVtxZ; //Histogram for position of primary vertex X
		TH1D *fHitV0A; //Histogram for hitting in V0A
		TH1D *fHitV0C; //Histogram for hitting in V0C
		TH1D *fMultRegionsMC; //Eta distribution for detector for MC
		TH1D *fSPDFiredChipsCluster; // Number of SPD fired Chips
		TH1D *fSPDFiredChipsHardware; // Number of SPD fired Chips
		TH1D *fSPDwc; //Which chip is fired
		TH1D *fRunOnline[8];
		TH1D *fRunVsMBOR_V0;
		TH1D *fRunVsMBAND_V0;
		TH1D *fRunVsMBOR_AD;
		TH1D *fRunVsMBAND_AD;
		TH1D *fRunVsMBOR_Global;
		TH1D *fRunVsMBAND_Global;
		TH1D *fRunVsDG_V0;
		TH1D *fRunVsDG_AD;
		TH1D *fRunVsDG_V0SPD;
		TH1D *fRunVsDG_ADSPD;
		TH1D *fRunVsDG_V0ADSPD;
		TH1D *fRunVsDG_V0ADFMDSPD;
		TH1D *fRunVsDG_V0ADFMDZDCSPD;
		TH1D *fRunVsDG_V0FMDSPD;
		TH1D *fRunVsDG_V0FMDZDCSPD;
		TH1D *fRunVs2t;
		TH1D *fRunVs2t_ITSSA;
		TH1D *fRunVs4t;
		TH1D *fRunVs4t_ITSSA;
		TH1D *fMultNG;
		TH1D *fMultNG_MS;
		TH1D *fMultDG;
		TH1D *fMultDG_MS;
		TH1D *fMassNG;
		TH1D *fMassNG_MS;
		TH1D *fMassDG;
		TH1D *fMassDG_MS;
		TH1D *fTrackCutsInfo;
		TH1D *fTrackCutsInfo_ITSSA;
		TH2D *fhClusterVsTracklets_bf;
		TH2D *fhClusterVsTracklets_af;
		TH1D *fMult_ITSSA;
		TH1D *fMult_DG_ITSSA;
		TH1D *fMult_NG_ITSSA;
		TH1D *fMult_ITSSA_MS;
		TH1D *fMult_DG_ITSSA_MS;
		TH1D *fMult_NG_ITSSA_MS;
		TH1D *fADATime_bf;
		TH1D *fADATime_af;
		TH1D *fADCTime_bf;
		TH1D *fADCTime_af;
		TH1D *fV0ATime_bf;
		TH1D *fV0ATime_af;
		TH1D *fV0CTime_bf;
		TH1D *fV0CTime_af;
		TH1D *fADATime_V0SPD;
		TH1D *fADCTime_V0SPD;
		TH1D *fADATime_V0ADSPD;
		TH1D *fADCTime_V0ADSPD;
		TH2D *fTPCSignal;
		TH2D *fTOFSignal;
		TH2D *fITSSignal;
		TH2D *fTRDSignal;
		TH1D *fMC_Eta[7];
		TH1D *fMC_DiffMass[7];
		TH1D *fMC_DiffMass_PDG[7];
		TH1D *fRunFiducial[6];
		TH1D *fMult_Gen;
		TH2D *fMult_Gen_Process;
		TH2D *fMult_Rec_DG_Process;
		TH2D *fMult_Rec_NG_Process;
		// -------------------------------------------------------------------

		ClassDef(AliAnalysisTaskCDPWA, 1);
};


#endif

