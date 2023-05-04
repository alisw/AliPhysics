/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskMCCorrections_H
#define AliAnalysisTaskMCCorrections_H

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


class AliAnalysisTaskMCCorrections : public AliAnalysisTaskSE {
	public:
		AliAnalysisTaskMCCorrections();
		AliAnalysisTaskMCCorrections(const char *name);
		virtual ~AliAnalysisTaskMCCorrections();

		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t *option);
		virtual void Terminate(Option_t *option);
		double GetFlatenicityV0();
		double GetMidRapidityMultiplicity();
		double GetFlatenicityMC();
		void MakeMCanalysisPID();
		void nSigmaContamination();
		void SetDataPeriod(std::string period="16k") { fPeriod = period; }
		void SetSystVarTrkCuts(const int SystVar=9) { fSystVarTrkCuts = SystVar; }
		void TrueINEL();
		void AccINEL();
		void SetPtMin(Double_t val) { fPtMin = val; } // Set pT cut for associated particles
		void SetNcl(UShort_t ncl) { fNcl = ncl; } // Set pT cut for associated particles
		void SetUseMC(Bool_t flat_flag = kFALSE) { fUseMC = flat_flag; } // use to analyse MC data
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
		Float_t fV0MMultiplicity;
		Bool_t fDeltaV0;
		Bool_t fRemoveTrivialScaling;
		AliPIDResponse* fPIDResponse;
		AliAnalysisFilter* fTrackFilter;
		AliAnalysisFilter* fTrackFilterPID;
		TList *fOutputList; //! output list in the root file
		Double_t fEtaCut;
		Double_t fPtMin;
		Double_t fNcl;
		TF1* fcutLow;
		TF1* fcutHigh;
		TF1* fcutDCAxy;
		int fSystVarTrkCuts;
		Double_t fv0mpercentile;
		Double_t fFlat;
		Double_t fFlatTPC;
		Float_t fFlatMC;
		AliMultSelection *fMultSelection;
		std::string fPeriod;
		TH2F *hFlatenicityMC;
		TH2F *hFlatenicityMCRec;
		TH2F *hFlatResponse;
		TProfile *hActivityV0McSect;
		TH2F *hFlatVsNchMC;
		TH1F* hTrueINEL_vtx;
		TH1F* hAccINEL_vtx;
		TH2F* hTrueINELWithFlat_evts;
		TH2F* hAccINELWithFlat_evts;

		TH2F* nsigma_kaon_h[4];
		TH2F* random_cont_in_kaon_h[4];
		TH2F* nsigma_proton_h[4];
		TH2F* random_cont_in_proton_h[4];
		TH2F* nsigma_pion_h[4];
		TH2F* random_cont_in_pion_h[4];

		TH3F* hPionTOFDCAxyNeg[3];
		TH3F* hProtonTOFDCAxyNeg[3];
		TH3F* hPionTOFDCAxyPos[3];
		TH3F* hProtonTOFDCAxyPos[3];
		TH3F* hPionTPCDCAxyNeg[3];
		TH3F* hProtonTPCDCAxyNeg[3];	
		TH3F* hPionTPCDCAxyPos[3];	
		TH3F* hProtonTPCDCAxyPos[3];

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
		TH3F* hTrueINELWithFlat_pT[4];
		TH3F* hAccINELWithFlat_pT[4];

		AliAnalysisTaskMCCorrections(
				const AliAnalysisTaskMCCorrections &); // not implemented
		AliAnalysisTaskMCCorrections &
			operator=(const AliAnalysisTaskMCCorrections &); // not implemented

		ClassDef(AliAnalysisTaskMCCorrections, 3);
};

#endif

