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
// Description:
//
// Select events according to gap conditions, analyze two-/four-track events in pp
// collisions at sqrt(s)=7 and 13TeV.
//
// This can be used for both lego-train and grid analysis
//
// Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>
//  continued by
//  Taesoo Kim <taesoo.kim@cern.ch>

#ifndef ALIANALYSISTASKCDPWA_H
#define ALIANALYSISTASKCDPWA_H

// ALICE class
class AliESDEvent;
class AliPIDResponse;
class AliPIDCombined;
class AliESDtrackCuts;
class AliESDAD;
class AliESDVZERO;

// ROOT class
class TH1D;
class TH2D;
class TList;
class TTree;
class TObject;
class TArrayI;

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "AliTriggerAnalysis.h"
#include "AliAnalysisUtils.h"
#include "AliMultiplicity.h"
#include "AliPID.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCDPWA : public AliAnalysisTaskSE
{
	public:
		AliAnalysisTaskCDPWA(const char* name);
		AliAnalysisTaskCDPWA();
		virtual ~AliAnalysisTaskCDPWA();

		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t *);
		virtual void Terminate(Option_t *);
		
		//Options for run macro
		void SetIsRun2(Bool_t isRun2) {fIsRun2 = isRun2;}//Run2 mode
		void SetIsMC(Bool_t MC) {fIsMC = MC;}//MC mode
		void SetCombinatoricsMode(Bool_t isComb) {fCombmode = isComb;}//kTRUE = save combinatorics
		void SetIsPythia8(Bool_t ispythia8) {fIsPythia8 = ispythia8;}//kTRUE = Pythia8 mode
		void SetRunSystematic(Bool_t isSys) {fIsSys = isSys;}//kTRUE = Run systematics
		void SetIsPWAMC(Bool_t isPWA) {fIsPWAMC = isPWA;}//kTRUE = for PWA in DRgen/DIME
		void SetMCSystematic(Int_t nSys) {fnSys = nSys;}

		struct EventInfo : public TObject{//Event information for PWA
			EventInfo()//Constructor
				:fPeriod(0)
				 ,fRunNumber(0)
				 ,fTimeStamp(0)
				 ,fSPDFired(0)
				 ,fV0Gap(0)
				 ,fADGap(0)
				 ,fFMDGap(0)
				 ,fZDCGap(0)
				 ,fNtrk_ST(0)
				 ,fNtrk_MS(0) {
					 for (Int_t i = 0; i < 11; i++) {
						 if (i < 3) {
							 fVertexSPD[i] = 0; 
							 fVertexTPC[i] = 0;
							 fVertexTracks[i] = 0;
							 fVertexUsed[i] = 0;
							 fSysVertex[i] = 0;
							 fSysPileUp[i] = 0;
						 }
						 if (i < 5) fSysCluster[i] = 0;
						 fCheckTwoTrack[i] = fCheckFourTrack[i] = 0;
					 }
				 }

			void Fill(const AliESDEvent *);

			UInt_t fPeriod;
			UInt_t fRunNumber;
			UInt_t fTimeStamp;
			Bool_t fSPDFired;
			Bool_t fV0Gap;
			Bool_t fADGap;
			Bool_t fFMDGap;
			Bool_t fZDCGap;
			Int_t fNtrk_ST;//Number of ITS+TPC Standard cut
			Int_t fNtrk_MS;//Number of ITS_TPC Martin's selection
			Bool_t fCheckTwoTrack[11];//Is two-track events? including systematics
			Bool_t fCheckFourTrack[11];//Is four-track events? including systematics
			Double_t fVertexSPD[3];//Vertex info
			Double_t fVertexTPC[3];//
			Double_t fVertexTracks[3];
			Double_t fVertexUsed[3];
			Bool_t fSysVertex[3];
			Bool_t fSysPileUp[3];
			Bool_t fSysCluster[5];

			ClassDef(EventInfo, 1);
		};

		struct CombInfo : public TObject{//Combinatorics
			CombInfo()//Constructor
				:fComb_IsPassMBOR(0)
				 ,fComb_IsPassVertex(0)
				 ,fComb_IsPileUp(0)
				 ,fComb_IsPassClusterCut(0)
				 ,fComb_SPDCluster(0)
				 ,fComb_SPDTracklets(0)
				 ,fComb_MC_EventProcess(0)
				 ,fComb_forwardP("TLorentzVector",2)
				 ,fComb_diffSystem("TLorentzVector",2) {
					 for (Int_t i = 0; i < 10; i++) {
						 fComb_DetHit[i] = 0;
					 }
				 }

			void Fill(const AliESDEvent *);

			Bool_t fComb_IsPassMBOR;//! Need for MC
			Bool_t fComb_IsPassVertex;//! Need for MC
			Bool_t fComb_IsPileUp;//!
			Bool_t fComb_IsPassClusterCut;//!
			UInt_t fComb_SPDCluster;//!
			UInt_t fComb_SPDTracklets;//!
			UInt_t fComb_MC_EventProcess;//!
			Bool_t fComb_DetHit[10];//!
			TClonesArray fComb_forwardP;//!
			TClonesArray fComb_diffSystem;//!

			ClassDef(CombInfo, 1);

		};

		struct TrackInfo : public TObject{//Track information for PWA
			TrackInfo(AliESDtrack *tr=NULL, AliPIDResponse *pidResponse=NULL, AliPIDCombined *pidCombined=NULL)//Constructor
				:TObject()
				 , fSign(0)
				 , fPx(0)
				 , fPy(0)
				 , fPz(0)
				 , fEnergy(0)
				 , fIntegratedLength(0)
				 , fITSSignal(0)
				 , fTPCSignal(0)
				 , fTOFSignal(0)
				 , fTRDSignal(0) {
					 for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
						 fITSSigma[i] = fTPCSigma[i] = fTOFSigma[i] = fITSSignalDelta[i] = fTPCBayesProb[i] = fTOFBayesProb[i] =
							 fITSBayesProb[i] = fTRDBayesProb[i] = fTPCProb[i] = fTOFProb[i] = fITSProb[i] = ftotBayesProb[i] = -999.;
						 fDetMask[i] = 0;
					 }
					 Fill(tr, pidResponse, pidCombined);
				 }

			void Fill(AliESDtrack *, AliPIDResponse *, AliPIDCombined *);

			Double_t fSign;
			Double_t fPx, fPy, fPz, fEnergy, fIntegratedLength;
			Double_t fITSSignal, fTPCSignal, fTOFSignal, fTRDSignal;
			Double_t fITSSigma[AliPID::kSPECIES];
			Double_t fTPCSigma[AliPID::kSPECIES];
			Double_t fTOFSigma[AliPID::kSPECIES];
			Double_t fITSSignalDelta[AliPID::kSPECIES];
			Bool_t fPIDStatus[5];
			Double_t fTPCBayesProb[AliPID::kSPECIES];
			Double_t fTOFBayesProb[AliPID::kSPECIES];
			Double_t fITSBayesProb[AliPID::kSPECIES];
			Double_t fTRDBayesProb[AliPID::kSPECIES];
			Double_t ftotBayesProb[AliPID::kSPECIES];
			Double_t fTPCProb[AliPID::kSPECIES];
			Double_t fTOFProb[AliPID::kSPECIES];
			Double_t fITSProb[AliPID::kSPECIES];
			UInt_t fDetMask[5];

			ClassDef(TrackInfo, 1);

		};

	private:
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
			kDGV0SPD, //!V0A && !V0C && SPD
			kDGADSPD, //!ADA && !ADC && SPD
			kDGV0ADSPD,
			kDGV0ADFMDSPD,
			kDGV0ADFMDZDCSPD,
			kDGV0FMDSPD,
			kDGV0FMDZDCSPD,
			kNG_Data,// No-gap?
			kSGA_Data,//Single-gap A-side
			kSGC_Data,//Single-gap C-side
			k2Tracks,
			k4Tracks,
			k2Tracks_ITSSA,
			k4Tracks_ITSSA,
			kAll
		};
		enum {
			// Trigger selection
			kTri_CINT5 = 0,
			kTri_CINT7,
			kTri_CINT10,
			kTri_CINT11,
			kTri_CADAND,
			kTri_C0SMB,
			kTri_CDG6,
			kTri_CDG6_SPD2,
			kTri_CDG7_SPD2,
			kTriAll
		};

		AliAnalysisTaskCDPWA(const AliAnalysisTaskCDPWA  &p);
		AliAnalysisTaskCDPWA& operator=(const AliAnalysisTaskCDPWA  &p);

		void PostOutputs();
		void SetBranchPWA(TTree *t);
		void SetBranchComb(TTree *t);
		Bool_t CheckInput();
		Bool_t CheckOnlineTrigger(const AliESDEvent *esd);
		void DoTimeMeasurements(const AliESDEvent *esd, const Int_t seq);
		Bool_t DoVertexCut(const AliESDEvent *esd);
		void DoCombinatorics(const AliESDEvent *esd);
		Bool_t DoMCPWA();

		//Member variables
		Bool_t fIsRun2;
		Bool_t fIsMC;
		Bool_t fCombmode;
		Bool_t fIsPythia8;
		Bool_t fIsSys;
		Bool_t fIsPWAMC;
		Bool_t fIsDGTrigger;
		Bool_t fnSys;

		//Determined in the .cxx
		Bool_t fIsPythia;
		Bool_t fIsPhojet;
		Bool_t fIsEPOS;

		// Output objects-----------------------------------------------------
		TTree *fTree; //! V0 2pion
		TTree *fTree_Comb;//! Combinatorics
		TList *fList; //! List for histogram

		TClonesArray fMCTrack;//!
		TClonesArray fTrackInfo;//!
		EventInfo fEventInfo;//!
		CombInfo fCombInfo;//!
		AliTriggerAnalysis fTriggerAnalysis;//!
		AliAnalysisUtils fAnalysisUtils;//!
		AliPIDResponse *fPIDResponse;//!
		AliPIDCombined *fPIDCombined;//!
		AliESDtrackCuts *fTrackCuts;//!
		AliESDEvent *fESDEvent;//!
		AliMultiplicity *fMult;//!

		//Used variable
		UInt_t fRunNumber;

		//Histograms----------------------------------------------------------
		TH1D *fHistEvent; //!Histogram for number of event
		TH1D *fHistTrigger;//!
		TH1D *fHistEventProcesses; //!Histogram for number of event(MC)
		TH1D *fHistPrimVtxX; //!Histogram for position of primary vertex X
		TH1D *fHistPrimVtxY; //!Histogram for position of primary vertex X
		TH1D *fHistPrimVtxZ; //!Histogram for position of primary vertex X
		TH1D *fMultRegionsMC; //!Eta distribution for detector for MC
		TH2D *fSPDFiredChipsClvsHd; //! Number of SPD fired Chips
		TH1D *fSPDwc; //!Which chip is fired
		TH1D *fRunOnline[9];//!
		TH1D *fRunVsMBOR_V0;//!
		TH1D *fRunVsMBAND_V0;//!
		TH1D *fRunVsMBOR_AD;//!
		TH1D *fRunVsMBAND_AD;//!
		TH1D *fRunVsMBOR_Global;//!
		TH1D *fRunVsMBAND_Global;//!
		TH1D *fRunVsDG[10];//!
		TH1D *fRunVs2t;//!
		TH1D *fRunVs2t_ITSSA;//!
		TH1D *fRunVs4t;//!
		TH1D *fRunVs4t_ITSSA;//!
		TH1D *fMultNG_ST;//!
		TH1D *fMultNG_MS;//!
		TH1D *fMultDG_ST;//!
		TH1D *fMultDG_MS;//!
		TH1D *fMassNG_ST_2t;//!
		TH1D *fMassNG_MS_2t;//!
		TH1D *fMassDG_ST_2t;//!
		TH1D *fMassDG_MS_2t;//!
		TH1D *fMassNG_ST_4t;//!
		TH1D *fMassNG_MS_4t;//!
		TH1D *fMassDG_ST_4t;//!
		TH1D *fMassDG_MS_4t;//!
		TH1D *fTrackCutsInfo;//!
		TH1D *fTrackCutsInfo_ITSSA;//!
		TH2D *fhClusterVsTracklets_bf;//!
		TH2D *fhClusterVsTracklets_af;//!
		TH1D *fMult_ITSSA;//!
		TH1D *fMult_DG_ITSSA;//!
		TH1D *fMult_NG_ITSSA;//!
		TH1D *fMult_ITSSA_MS;//!
		TH1D *fMult_DG_ITSSA_MS;//!
		TH1D *fMult_NG_ITSSA_MS;//!
		TH2D *fTPCSignal;//!
		TH2D *fTOFSignal;//!
		TH2D *fITSSignal;//!
		TH2D *fTRDSignal;//!
		TH1D *fMC_Eta[7];//!
		TH1D *fMC_DiffMass[7];//!
		TH1D *fMC_DiffMass_PDG[7];//!
		TH1D *fRunFiducial[6];//!
		TH1D *fMult_Gen;//!
		TH2D *fMult_Gen_Process;//!
		TH2D *fMult_Rec_DG_Process;//!
		TH2D *fMult_Rec_NG_Process;//!
		TH2D *fV0Time[2];//!
		TH2D *fADTime[2];//!
		// -------------------------------------------------------------------

		ClassDef(AliAnalysisTaskCDPWA, 1);
};


#endif

