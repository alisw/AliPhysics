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
// collisions
//
// Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>
//  continued by
//  Taesoo Kim <taesoo.kim@cern.ch>

#ifndef ALIANALYSISTASKCDTREE_H
#define ALIANALYSISTASKCDTREE_H

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#endif

#include "THnSparse.h" // forward declaration is not possible

class AliESDEvent;
class AliVTrack;
class AliPIDResponse;
class AliPhysicsSelection;
class AliESDtrackCuts;
class AliPIDCombined;

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

class AliAnalysisTaskCDTree : public AliAnalysisTaskSE
{
	public:
		AliAnalysisTaskCDTree(const char* name);
		AliAnalysisTaskCDTree();
		virtual ~AliAnalysisTaskCDTree();

		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t *);
		virtual void Terminate(Option_t *);

	private:
		enum {
			//gapcondition
			kBinDG = 1, // double gap
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
			kBinMCAll
		};
		enum {kInput = 0, kMCCheck, kOfflineCut, kVtxCut, kPileUpCut, kClusterCut,kMBOR,kMBAND,kNG,kGA,kGC,kDG,kNG_wFMD,kGA_wFMD,kGC_wFMD,kDG_wFMD,kCheckDG,kCheckV0Hit,k2Tracks,k4Tracks,k6Tracks,kAll};

		AliAnalysisTaskCDTree(const AliAnalysisTaskCDTree  &p);
		AliAnalysisTaskCDTree& operator=(const AliAnalysisTaskCDTree  &p);


		// functions called by the UserExec(...), not to be called elsewhere!
		//-------------------------------------------------------------------
		Bool_t CheckInput();
		void PostOutputs(); // cares about posting the output before exiting UserExec,
		// WARNING: PostOutputs should only be used directly in the UserExec!!
		static Bool_t SPDLoc2Glo(const Int_t id, const Double_t *loc, Double_t *glo);
		static Int_t CheckChipEta(const Int_t chipKey, const TString scut,
				const Double_t vtxPos[], TH2 *hitMapSPDinner,
				TH2 *hitMapSPDouter);
		static Bool_t CutEvent(const AliESDEvent *ESDEvent, TH1 *hspd, TH1* hfochans,
				TH1* hpriVtxX, TH1* hpriVtxY, TH1* hpriVtxZ, Double_t *hVertex); 
		static Int_t GetFastORmultiplicity(const AliESDEvent *ESDEvent);
		Bool_t CheckV0Hit(const AliESDEvent *ESDEvent, TH1D* fHitV0A, TH1D* fHitV0C);
		Int_t DetermineGap(const AliESDEvent *ESDEvent, const Bool_t wCent, const Bool_t wFMD); // determines the gap of all available detectors
		Int_t DoMCTruth(); // analyses the MCtruth for corrections, returns #Primaries

		// analyzes track pairs in MC truth
		void DetermineMCprocessType(); // determines the MC process ID

		// analysis task status
		//---------------------
		Double_t fMaxVtxDst; // maximum distance of the track and SPD vertex

		// event information
		//------------------
		AliPIDResponse *fPIDResponse; //! PID Response object
		AliPIDCombined *fPIDCombined; //!
		AliESDEvent *fESDEvent; // esd event object
		AliPhysicsSelection *fPhysicsSelection; // physics selection object
		AliESDtrackCuts *fTrackCuts;//For standard
		AliESDtrackCuts *fTrackCuts_ITSSA;//For ITSSA
		Double_t fVtxDst; // distance of the primary vertex from tracks and from SPD
		Double_t fVtxZ; // z-position of the primary vertex from tracks
		Int_t fMCprocessType; // MC process type, 0 for data
		Int_t fMCprocess; // detailed MC sub process information

		// information about the trackpair which is currently processed
		Int_t fRun; // number of the run which is about to be processed
		Int_t fPIDmode; // selects set of PID cuts, 0 for 3sigma standard cuts,
		Int_t fGapInformation; // same as above, without the requirement
		// on central activity

		// Output objects-----------------------------------------------------
		TTree *fTree; //! V0 2pion
		TList *fList; //! List for histogram

		//Two track
		Bool_t fCheckTwoPion; //! Check that this event has 2 tracks
		Bool_t fCheckFourPion; //! Check that this event has 4 tracks
		Bool_t fCheckV0FMD; //! Check that this event verified by V0 and FMD
		Bool_t fCheckTwoPion_ITSSA; //!for 2tracks with ITSSA
		Bool_t fIsMC;//!
		Double_t fVertex[3];//!
		Double_t fTwoPionTrack[2][9];//! First track Momentum, Energy and Sign
		Double_t fTwoPionTPCSigma[2][9];//! TPC PID
		Double_t fTwoPionTOFSigma[2][9];//! TOF PID
		Double_t fTwoPionITSSigma[2][9];//! ITS PID
		UInt_t fTwoPionMask_TPC[2];//! TPC Mask
		UInt_t fTwoPionMask_TOF[2];//! TOF Mask
		UInt_t fTwoPionMask_ITS[2];//! ITS Mask
		UInt_t fTwoPionMask_TRD[2];//! TRD Mask
		UInt_t fTwoPionMask_tot[2];//! Combination Mask
		UInt_t fTwoPionDetMask_tot[2];//!
		Double_t fTwoPionBayesProb_TPC[2][5];//! Bayesian probabilities with TPC
		Double_t fTwoPionBayesProb_TOF[2][5];//! Bayesian probabilities with TOF
		Double_t fTwoPionBayesProb_ITS[2][5];//! Bayesian probabilities with ITS
		Double_t fTwoPionBayesProb_TRD[2][5];//! Bayesian probabilities with TRD
		Double_t fTwoPionBayesProb_tot[2][5];//! Bayesian probabilities with combination
		Double_t fTwoPionTrack_ITSSA[2][9]; //! First track Momentum, Energy and Sign

		//Four track
		Double_t fFourPionTrack[4][9]; //! First track Momentum, Energy and Sign
		Double_t fFourPionTPCSigma[4][9];//! TPC PID
		Double_t fFourPionTOFSigma[4][9];//! TOF PID
		Double_t fFourPionITSSigma[4][9];//! ITS PID
		UInt_t fFourPionMask_TPC[4];//! TPC Mask
		UInt_t fFourPionMask_TOF[4];//! TOF Mask
		UInt_t fFourPionMask_ITS[4];//! ITS Mask
		UInt_t fFourPionMask_TRD[4];//! TRD Mask
		UInt_t fFourPionMask_tot[4];//! Combination Mask
		UInt_t fFourPionDetMask_tot[4];//!
		Double_t fFourPionBayesProb_TPC[4][5];//! Bayesian probabilities with TPC
		Double_t fFourPionBayesProb_TOF[4][5];//! Bayesian probabilities with TOF
		Double_t fFourPionBayesProb_ITS[4][5];//! Bayesian probabilities with ITS
		Double_t fFourPionBayesProb_TRD[4][5];//! Bayesian probabilities with TRD
		Double_t fFourPionBayesProb_tot[4][5];//! Bayesian probabilities with combination
		Double_t fFourPionTrack_ITSSA[4][9]; //! First track Momentum, Energy and Sign

		Double_t fMCGenProtonTrack[2][5]; //! Info of generated Proton 1 after scattering
		Double_t fMCGenPionTrack[2][5]; //! Info of generated Pion 1 after scattering
		Int_t fRunNumber;//!
		Int_t fPeriod;//!

		// Histogram objects---------------------------------------------------
		TH1D *fHistEvent; //Histogram for number of event
		TH1D *fHistEventProcesses; //Histogram for number of event(MC)
		TH1D *fHistPrimVtxX; //Histogram for position of primary vertex X
		TH1D *fHistPrimVtxY; //Histogram for position of primary vertex X
		TH1D *fHistPrimVtxZ; //Histogram for position of primary vertex X
		TH1D *fHitV0A; //Histogram for hitting in V0A
		TH1D *fHitV0C; //Histogram for hitting in V0C
		TH1D *fMultRegionsMC; //Eta distribution for detector for MC
		TH1D *fHistSPDFiredChips; // Number of SPD fired Chips
		TH1D *fRunVsMBOR;
		TH1D *fRunVsMBAND;
		TH1D *fRunVsDG;
		TH1D *fRunVsDG_wFMD;
		TH1D *fRunVs2t;
		TH1D *fRunVs2t_ITSSA;
		TH1D *fRunVs4t;
		TH1D *fRunVs4t_ITSSA;
		TH1D *fMultNG;
		TH1D *fMultNG_MS;
		TH1D *fMultDG;
		TH1D *fMultDG_MS;
		TH1D *fMassNG;
		TH1D *fMassNG_st2t;
		TH1D *fMassNG_MS;
		TH1D *fMassDG;
		TH1D *fMassDG_st2t;
		TH1D *fMassDG_MS;
		TH1D *fTrackCutsInfo;
		TH2D *fSPDTrkvsCls_bf;
		TH2D *fSPDTrkvsCls_af;
		TH1D *fPrimaries;
		// -------------------------------------------------------------------

		ClassDef(AliAnalysisTaskCDTree, 1);
};


#endif

