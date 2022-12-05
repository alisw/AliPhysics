/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskFlatenicityPiKp_H
#define AliAnalysisTaskFlatenicityPiKp_H

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


class AliAnalysisTaskFlatenicityPiKp : public AliAnalysisTaskSE {
	public:
		AliAnalysisTaskFlatenicityPiKp();
		AliAnalysisTaskFlatenicityPiKp(const char *name);
		virtual ~AliAnalysisTaskFlatenicityPiKp();

		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t *option);
		virtual void Terminate(Option_t *option);
		double GetFlatenicityV0();
		double GetMidRapidityMultiplicity();
		double GetFlatenicityTPC();
		double GetFlatenicityMC();
		double EtaCalibration(const double&);
		void MakeMCanalysis();
		void MakeMCanalysisPID();
		void nSigmaContamination();
		void MakeDataanalysis();
		void MakePIDanalysis();
		void SetPtMin(Double_t val) { fPtMin = val; } // Set pT cut for associated particles
		void SetNcl(UShort_t ncl) { fNcl = ncl; } // Set pT cut for associated particles
		void SetUseMC(Bool_t flat_flag = kFALSE) { fUseMC = flat_flag; } // use to analyse MC data
		void SaveDCAxyHistograms(Bool_t saveDCAxy ) { fSaveDCAxyHistograms = saveDCAxy; } 
		void IsV0MCalibrated(Bool_t V0MCal = kTRUE) { fV0MEqualisation = V0MCal; }
		void IsdEdxCalibrated(Bool_t dEdxCal = kTRUE) { fdEdxCalibrated = dEdxCal; }
		void SetDetectorForFlatenicity(TString det = "V0") { fDetFlat = det; }
		void SaveThisMultBin(std::string v0bin = "0_1") { fV0MBin = v0bin; }
		void SetDataPeriod(std::string period = "16l") { fPeriod = period;}
		void SetMCclosureTest(Bool_t flat_flag = kFALSE) { fIsMCclosure = flat_flag; }
		void SetDeltaV0(Bool_t deltav0 = kFALSE) { fDeltaV0 = deltav0; }
		void SetRemoveTrivialScaling(Bool_t flat_flag = kFALSE) { fRemoveTrivialScaling = flat_flag; }
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
		Bool_t fUseMC;      // analyze MC events
		Int_t fV0Mindex;
		Float_t fV0MMultiplicity;
		TString fDetFlat;
		std::string fV0MBin;
		Bool_t fIsMCclosure;
		Bool_t fDeltaV0;
		Bool_t fRemoveTrivialScaling;
		Int_t fnGen;
		AliPIDResponse* fPIDResponse;
		AliAnalysisFilter* fTrackFilter;
		AliAnalysisFilter* fTrackFilterPID;
		TList *fOutputList; //! output list in the root file
		Double_t fEtaCut;
		Double_t fPtMin;
		Double_t fNcl;
		Bool_t fV0MEqualisation;
		Bool_t fdEdxCalibrated;
		Bool_t fSaveDCAxyHistograms;
		TF1* fEtaCalibrationPos;
		TF1* fEtaCalibrationNeg;
		TF1* fV0CCalibration;
		TF1* fV0ACalibration;
		TF1* fcutLow;
		TF1* fcutHigh;
		TF1* fcutDCAxy;
		std::string fPeriod;
		Double_t ftrackmult08;
		Double_t fv0mpercentile;
		Double_t fMidRapidityMult;
		Double_t fFlat;
		Double_t fFlatTPC;
		Float_t fFlatMC;
		AliMultSelection *fMultSelection;
		TH1D *hPtPrimIn;
		TH1D *hPtPrimOut;
		TH1D *hPtSecOut;
		TH1D *hPtOut;
		TH2F *hFlatenicityMC;
		TH2F *hFlatenicityMCRec;
		TH2D *hFlatResponse;
		TH2D *hFlatVsPtMC;
		TH3F *hActivityV0CV0A;
		TProfile *hActivityV0DataSect;
		TProfile *hActivityV0McSect;
		TH2D *hFlatVsNchMC;
		TH2F *hFlatVsV0MVsMult[9];
		TH2F *hFlatVsPtV0M[9];
		TH3F *hNsigmaPiPos[9][4];
		TH3F *hNsigmaKPos[9][4];
		TH3F *hNsigmaPPos[9][4];
		TH2F *hPtTPCEtaPos[9][4];
		TH3F *hNsigmaPiNeg[9][4];
		TH3F *hNsigmaKNeg[9][4];
		TH3F *hNsigmaPNeg[9][4];
		TH2F *hPtTPCEtaNeg[9][4];
		TH3F* hBetaPos[9][4];
		TH2F* hMomentumTOFEtaPos[9][4];
		TH2F* hPtTOFEtaPos[9][4];
		TH3F* hBetaNeg[9][4];
		TH2F* hMomentumTOFEtaNeg[9][4];
		TH2F* hPtTOFEtaNeg[9][4];
		TH3F* hdEdx[9][4];
		TH2F* hPtrTPC[9][4];
		TH2F* hPtVsP[4];

		TH2F* random_cont_in_kaon_h[4];
		TH2F* pion_cont_in_kaon_h[4];
		TH2F* electron_cont_in_kaon_h[4];
		TH2F* nsigma_kaon_h[4];
		TH2F* random_cont_in_proton_h[4];
		TH2F* pion_cont_in_proton_h[4];
		TH2F* electron_cont_in_proton_h[4];
		TH2F* nsigma_proton_h[4];
		TH2F* random_cont_in_pion_h[4];
		TH2F* kaon_cont_in_pion_h[4];
		TH2F* electron_cont_in_pion_h[4];
		TH2F* nsigma_pion_h[4];

		TH2F* hPionTOFDCAxyNeg[3];
		TH2F* hProtonTOFDCAxyNeg[3];
		TH2F* hPionTOFDCAxyPos[3];
		TH2F* hProtonTOFDCAxyPos[3];
		TH2F* hPionTPCDCAxyNeg[3];
		TH2F* hProtonTPCDCAxyNeg[3];	
		TH2F* hPionTPCDCAxyPos[3];	
		TH2F* hProtonTPCDCAxyPos[3];

		TH1F* hMCPtPionPos;
		TH1F* hMCPtKaonPos;
		TH1F* hMCPtProtonPos;
		TH1F* hMCPtPionNeg;
		TH1F* hMCPtKaonNeg;
		TH1F* hMCPtProtonNeg;

		TH1F* hTPCRecTracksPionPos;
		TH1F* hTPCRecTracksKaonPos;
		TH1F* hTPCRecTracksProtonPos;
		TH1F* hTPCRecTracksPionNeg;
		TH1F* hTPCRecTracksKaonNeg;
		TH1F* hTPCRecTracksProtonNeg;

		TH1F* hTOFRecTracksPionPos;
		TH1F* hTOFRecTracksKaonPos;
		TH1F* hTOFRecTracksProtonPos;
		TH1F* hTOFRecTracksPionNeg;
		TH1F* hTOFRecTracksKaonNeg;
		TH1F* hTOFRecTracksProtonNeg;
		TH1F* hrTPCRecTracksPion;
		TH1F* hrTPCRecTracksKaon;
		TH1F* hrTPCRecTracksProton;
		TH2F* hPionTPCDCAxyNegData;
		TH2F* hPionTPCDCAxyPosData;
		TH2F* hProtonTPCDCAxyNegData; 
		TH2F* hProtonTPCDCAxyPosData;
		TH2F* hPionTOFDCAxyNegData;
		TH2F* hPionTOFDCAxyPosData;
		TH2F* hProtonTOFDCAxyNegData; 
		TH2F* hProtonTOFDCAxyPosData;
		TH2F* hMIPVsEta;
		TProfile* pMIPVsEta;
		TH2F* hPlateauVsEta;
		TProfile* pPlateauVsEta;

		AliAnalysisTaskFlatenicityPiKp(
				const AliAnalysisTaskFlatenicityPiKp &); // not implemented
		AliAnalysisTaskFlatenicityPiKp &
			operator=(const AliAnalysisTaskFlatenicityPiKp &); // not implemented

		ClassDef(AliAnalysisTaskFlatenicityPiKp, 3);
};

#endif
