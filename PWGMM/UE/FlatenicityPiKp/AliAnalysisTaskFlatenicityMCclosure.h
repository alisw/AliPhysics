/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskFlatenicityMCclosure_H
#define AliAnalysisTaskFlatenicityMCclosure_H

class AliESDtrackCuts;
class AliESDEvent;
class AliESDAD;
class TList;
class TH1D;
class TH2D;
class TH1I;
class TH3F;
class TProfile;
class THnSparse;

#include "AliAnalysisTaskSE.h"

#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"

class AliAnalysisTaskFlatenicityMCclosure : public AliAnalysisTaskSE {
	public:
		AliAnalysisTaskFlatenicityMCclosure();
		AliAnalysisTaskFlatenicityMCclosure(const char *name);
		virtual ~AliAnalysisTaskFlatenicityMCclosure();

		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t *option);
		virtual void Terminate(Option_t *option);
		Double_t GetFlatenicityV0();
		Double_t GetFlatenicityV0EqualALICE();
		Double_t GetFlatenicityTPC();
		Double_t GetFlatenicityMC();
		void ExtractMultiplicities();
		void ExtractMultiplicitiesEqualALICE();
		void ExtractMultiplicitiesMC();
		void MCkinematics();
		void MakeMCanalysis();
		void MakeDataanalysis();
		void GetMCchargedTrueDists();
		int GetPidCode(int pdgCode);
		bool HasRecVertex();

		void SetPtMin(Double_t val) { fPtMin = val; }
		void SetUseMC(Bool_t flat_flag = kFALSE) { fUseMC = flat_flag; }
		void SetDetectorForFlatenicity(TString det = "V0") { fDetFlat = det; }
		void SetRemoveTrivialScaling(Bool_t flat_flag = kFALSE) { fRemoveTrivialScaling = flat_flag; }
		void SetV0Calib(Bool_t calib_flag = kFALSE) { fIsCalib = calib_flag; }
		void SetEqualV0Alice(Bool_t calib_flag = kFALSE) {
			fIsEqualALICE = calib_flag;
		}

	private:
		AliESDEvent *fESD; //! input ESD event
		AliEventCuts fEventCuts;
		AliStack *fMCStack; //! MC stack
		AliMCEvent *fMC;    //! MC Event
		Bool_t fUseMC;      // analyze MC events
		Bool_t fIsCalib;
		Bool_t fIsEqualALICE;
		Float_t fVtxz;
		TF1 *fParVtx;
		Int_t fV0Mindex;
		Float_t fmultTPC;
		int fmultV0A;
		int fmultV0C;
		Float_t fmultADA;
		Float_t fmultADC;
		Float_t fmultTPCmc;
		Float_t fmultV0Amc;
		Float_t fmultV0Cmc;
		Float_t fmultADAmc;
		Float_t fmultADCmc;
		TString fDetFlat;
		Bool_t fRemoveTrivialScaling;
		AliAnalysisFilter *fTrackFilter;
		TList *fOutputList;
		Double_t fEtaCut;
		Double_t fPtMin;
		Double_t fv0mamplitude;
		Double_t fv0mpercentile;
		Float_t fFlat;
		Float_t fFlatMC;
		AliMultSelection *fMultSelection;
		TH2D *hFlatV0vsFlatTPC;
		TH1D *hFlatenicityMC;
		TH2D *hFlatResponse;
		TProfile *hActivityV0DataSectBefore;
		TProfile *hActivityV0DataSect;
		TProfile *hV0vsVtxz;
		TProfile *hActivityV0McSect;
		TH2D *hFlatVsNchMC;
		TH2D *hFlatVsV0M;
		TH2D *hFlatMCVsV0M;
		TH1D *hEtamc;
		TH2D *hMultMCmVsV0M;
		TH2D *hMultMCaVsV0M;
		TH2D *hMultMCcVsV0M;
		TH2D *hMultmVsV0M;
		TH2D *hMultmVsV0Malice;
		TH2D *hMultaVsV0M;
		TH2D *hMultcVsV0M;
		TH1D *hV0MBadruns;
		TH2F* hMultV0AV0CvsFlat_BFTrigSel;
	       	TH2F* hAmpV0AV0CvsFlat_BFTrigSel;
	        TH2F* hPercentileV0MvsFlat_BFTrigSel;
		TH2F* hMultV0AV0CvsFlat_AFTrigSel;
	       	TH2F* hAmpV0AV0CvsFlat_AFTrigSel;
	        TH2F* hPercentileV0MvsFlat_AFTrigSel;
		TH3F *hMultV0AV0CvsFlatvspT_BFTrigSel;
		TH3F *hAmpV0AV0CvsFlatvspT_BFTrigSel;
		TH3F *hPercentileV0MvsFlatvspT_BFTrigSel; 
		TH3F *hMultV0AV0CvsFlatvspT_AFTrigSel;
		TH3F *hAmpV0AV0CvsFlatvspT_AFTrigSel;
		TH3F *hPercentileV0MvsFlatvspT_AFTrigSel; 
		TH2F *hAmpV0AV0CvsFlat;
		TH2F *hPercentileV0MvsFlat; 
		TH3F *hAmpV0AV0CvsFlatvspT;
		TH3F *hPercentileV0MvsFlatvspT; 
		TH1D *hPtAll;
		TH1D *hPtPrimaries;
		TH1D *hPtSecondaries;
		TH2F *hAmpV0AV0CalicevsFlat;
		TH3F *hAmpV0AV0CalicevsFlatvspT; 
		TH2D *hComponentsMult[4];
		TH2D *hCombinedMult[3];
		TH2D *hComponentsMultmc[4];
		TH2D *hCombinedMultmc[3];
		TH2D *hRmCombinedMult[3];
		TH2D *hMultMCmVsFlat[9];
		TH2D *hMultmVsFlat[9];
		TH3F *hAmpV0AV0CvsFlatvspT_pi; 
		TH3F *hAmpV0AV0CvsFlatvspT_k; 
		TH3F *hAmpV0AV0CvsFlatvspT_p; 
		TH3F *hMultV0AV0CvsFlatvspT_pi_BFTrigSel;
		TH3F *hMultV0AV0CvsFlatvspT_k_BFTrigSel;
		TH3F *hMultV0AV0CvsFlatvspT_p_BFTrigSel;

		AliAnalysisTaskFlatenicityMCclosure(
				const AliAnalysisTaskFlatenicityMCclosure &); // not implemented
		AliAnalysisTaskFlatenicityMCclosure &
			operator=(const AliAnalysisTaskFlatenicityMCclosure &); // not implemented

		ClassDef(AliAnalysisTaskFlatenicityMCclosure, 3);
};

#endif
