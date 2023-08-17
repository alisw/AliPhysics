/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskdEdxCalibration_H
#define AliAnalysisTaskdEdxCalibration_H

class AliESDtrackCuts;
class AliESDEvent;
class TList;
class TH1D;
class TH2D;
class TH3D;
class TH3F;
class TH1I;
class TF1;
class TProfile;
class THnSparse;

#include "AliAnalysisTaskSE.h"

#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"

#include <iostream>
#include <string>
#include <map>


class AliAnalysisTaskdEdxCalibration : public AliAnalysisTaskSE {
	public:
		AliAnalysisTaskdEdxCalibration();
		AliAnalysisTaskdEdxCalibration(const char *name);
		virtual ~AliAnalysisTaskdEdxCalibration();

		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t *option);
		virtual void Terminate(Option_t *option);
		double EtaCalibration(const double&);
		double EtaCalibrationEl(const double&);
		double GetFlatenicityV0();
		void dEdxMIP();
		void AnalyzeV0s();
		void SetPtMin(Double_t val) { fPtMin = val; } // Set pT cut for associated particles
		void SetNcl(UShort_t ncl) { fNcl = ncl; } // Set pT cut for associated particles
		void SetUseMC(Bool_t flat_flag = kFALSE) { fUseMC = flat_flag; } // use to analyse MC data
		void IsdEdxCalibrated(Bool_t dEdxCal = kTRUE) { fdEdxCalibrated = dEdxCal; }
		void SetDataPeriod(std::string period = "16l") { fPeriod = period;}
		bool PhiCut(Double_t , Double_t , Double_t , Float_t , TF1* , TF1* );
		bool HasRecVertex();
		bool TOFPID(AliESDtrack*);
		int GetPidCode(Int_t);


	protected:
	private:
		AliESDEvent *fESD; //! input ESD event
		AliEventCuts fEventCuts;
		AliStack *fMCStack; //! MC stack
		AliMCEvent *fMC;    //! MC Event
		bool fUseMC;      // analyze MC events
		AliPIDResponse* fPIDResponse;
		AliAnalysisFilter* fTrackFilter;
		AliAnalysisFilter* fTrackFilterPID;
		TList *fOutputList; //! output list in the root file
		double fEtaCut;
		double fPtMin;
		double fNcl;
		bool fdEdxCalibrated;
		TF1* fEtaCalibrationPos;
		TF1* fEtaCalibrationNeg;
		TF1* fcutLow;
		TF1* fcutHigh;
		TF1* fcutDCAxy;
		std::string fPeriod;
		double fv0mpercentile;
		double fFlat;
		AliMultSelection *fMultSelection;

		TH3F* hPionTPCDCAxyNegData;
		TH3F* hPionTPCDCAxyPosData;
		TH3F* hProtonTPCDCAxyNegData; 
		TH3F* hProtonTPCDCAxyPosData;
		TH3F* hPionTOFDCAxyNegData;
		TH3F* hPionTOFDCAxyPosData;
		TH3F* hProtonTOFDCAxyNegData; 
		TH3F* hProtonTOFDCAxyPosData;

		TH2F* hMIPVsEta;
		TProfile* pMIPVsEta;
		TH2F* hPlateauVsEta;
		TProfile* pPlateauVsEta;
		TH2F* hMIPVsEtaV0s;
		TProfile* pMIPVsEtaV0s;

		//! Histograms for V0s
		TH2F* histPiV0[4];
		TH2F* histPV0[4];
		TH2F* histEV0[4];
		TH2F* histPiTof[4];
		TH2F* histPiTof2[4];
		TH2F* hBetaPi[4];
		TH2F* hBetaP[4];

		AliAnalysisTaskdEdxCalibration(
				const AliAnalysisTaskdEdxCalibration &); // not implemented
		AliAnalysisTaskdEdxCalibration &
			operator=(const AliAnalysisTaskdEdxCalibration &); // not implemented

		ClassDef(AliAnalysisTaskdEdxCalibration, 3);
};

#endif

