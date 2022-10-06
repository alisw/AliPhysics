#ifndef AliAnalysisTaskSpectraMC_H
#define AliAnalysisTaskSpectraMC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id$ */


// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TAxis.h>


#include <TProfile.h>
#include <TTreeStream.h>
#include <TRandom.h>
#include <TObject.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliAnalysisFilter.h>
#include <AliVHeader.h>
#include <AliAODMCParticle.h> 
#include <AliESDtrackCuts.h>
#include <AliPIDResponse.h>
#include "AliTPCPIDResponse.h"
#include <AliEventCuts.h>
#include "AliVTrack.h"
#include "TVirtualMCApplication.h"
#include "TVirtualMC.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliGenEventHeader.h"

using namespace std;


class AliAnalysisTaskSpectraMC : public AliAnalysisTaskSE
{
	public:


		AliAnalysisTaskSpectraMC();
		AliAnalysisTaskSpectraMC(const char *name);
		virtual ~AliAnalysisTaskSpectraMC();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);

		Bool_t   GetAnalysisMC() { return fAnalysisMC; }
		Double_t GetEtaCut() { return fEtaCut; }

		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(bool isMC) {fAnalysisMC = isMC;}
		virtual void  SetMCClosure(bool isMCclos) {fIsMCclosure = isMCclos;}
		virtual void  SetTrackID(const int trkID){fTrackID = trkID;}
		virtual void  SetNcl(const Int_t ncl){fNcl = ncl;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		//virtual void  SetTrackCutsType(bool isTPCOnlyTrkCuts) { fSetTPConlyTrkCuts = isTPCOnlyTrkCuts; }
		virtual void  SetHybridTracks(bool isSelectHybrid) { fSelectHybridTracks = isSelectHybrid; }
		AliESDtrack*  SetHybridTrackCuts(AliESDtrack *track, const bool fill1, const bool fill2, const bool fill3);

	private:


		AliESDtrack* GetLeadingTrack();
		void GetLeadingObject(Bool_t isMC);
		void GetMultiplicityDistributions();
		void GetDetectorResponse();
		void GetMCCorrections();
		void GetMultiplicityDistributionsPhi();
		void GetMCCorrectionsPhi();
		virtual Double_t DeltaPhi(Double_t phi, Double_t lphi,
				Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
		short   GetPidCode(Int_t pdgCode) const;

		bool selectVertex2015pp(AliESDEvent* esd, Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
		bool IsGoodSPDvertexRes(const AliESDVertex* spdVertex = NULL);
		bool IsGoodZvertexPos(AliESDEvent *esd);
		bool PhiCut(const double& pt, double phi, const double& q, const float& mag, TF1* phiCutLow, TF1* phiCutHigh);
		float GetMaxDCApTDep( TF1 *fcut, double pt);
	///	virtual void SetTrackCuts(AliAnalysisFilter* fTrackFilter);
		double EtaCalibration(const int &centrality, const double &Eta);
		double EtaCalibrationEl(const int &centrality, const double &Eta);
		bool TOFPID(AliESDtrack* track);

		static const Double_t fgkClight;   // Speed of light (cm/ps)

		AliESDEvent* fESD;                  //! ESD object
		AliEventCuts fEventCuts;
		AliMCEvent*  fMC;                   //! MC object
		AliStack*    fMCStack;              //! MC ESD stack
		TClonesArray* fMCArray;             //! MC array for AOD
		AliPIDResponse* fPIDResponse;       //! Pointer to PIDResponse
		AliESDtrackCuts* fGeometricalCut;    //  Track Filter, set 2010 with golden cuts
		AliESDtrackCuts* fTrackFilter;
		AliESDtrackCuts*   fHybridTrackCuts1;                 //  Track cuts for tracks without SPD hit
	        AliESDtrackCuts*   fHybridTrackCuts2;                 //  Track cuts for tracks witout SPD hit or ITS refit
		AliAnalysisUtils* utils;
		TString     fAnalysisType;        //  "ESD" or "AOD"
		bool        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
		bool        fIsMCclosure;          //  Real(kFALSE) or MC(kTRUE) flag
		TRandom*    fRandom;              //! random number generator

		//
		// Cuts and options
		//

		int fTrackID;
		int        fNcl;
		double     fEtaCut;             // Eta cut used to select particles
		const Double_t fDeDxMIPMin;
		const Double_t fDeDxMIPMax;
		const Double_t fdEdxHigh;
		const Double_t fdEdxLow;
		//bool fSetTPConlyTrkCuts;
		bool fSelectHybridTracks;

		const double fLeadPtCutMin;
		const double fLeadPtCutMax;
		double fGenLeadPhi;
		double fGenLeadPt;
		double fGenLeadEta;
		double fRecLeadEta;
		int    fGenLeadIn;
		double fRecLeadPhi;
		double fRecLeadPt;
		int    fRecLeadIn;
		const double fPtMin;

		//
		// Output objects
		//

		TList*        fListOfObjects;     //! Output list of objects
		TH1F*         fEvents;            //! No of accepted events
		

		// Tracking & Matching efficiencies
		TH1F* hPtGenIn[5];
		TH1F* hPtRecIn[5];
		TH1F* hPtGenPosIn[4];
		TH1F* hPtRecPosIn[4];
		TH1F* hPtGenNegIn[4];
		TH1F* hPtRecNegIn[4];
		TH1F* hPtRecInTOF[4];
		TH1F* hPtRecPosInTOF[4];
		TH1F* hPtRecNegInTOF[4];
		TH1F* hPtrTPCRecIn[4];

		// To Unfold
		TH2F* hNchVsPtDataTPC[4][5];

		// For closure
		TH2F* hNchGenVsPtRec[4][5];
		TH2F* hNchGenVsPtGenPID[4][5];

		TH1F* hPtRec;
		TH1F* hPtRecPion;
		TH1F* hPtRecProton;
		TH1F* hPtPri;
		TH1F* hPtSec;
		TH1F* hInvMassPhi;
		TH2F* fPtLVsNchRec; 
 		TH1F* hDeltaPtLeading;
		TH2F* hPtRecEtaRecLeading;
		TH2F* hEtaGenEtaRecLeading;
		TH2F* hPhiEtaRecLeading;
		TH1F* hPhiLeading;


		TH1F* hPhiGen[4];
		TH1F* hPhiRec[4];

		TH1F* hMultTSGen;
		TH1F* hMultTSRec;
		TH1F* hTrigger;

		// Response matrices
		TH2F* hNchResponse;

		TH2F* hPhiTotal;
		TH1F* hPhiResPhi;
		TH2F* hPhiStandard;
		TH2F* hPhiHybrid1;
		TH2F* hPhiHybrid2;

		TF1* fEtaCalibration;
		TF1* fEtaCalibrationEl;
		TF1* fcutDCAxy;
		TF1* fcutLow;
		TF1* fcutHigh;

		AliAnalysisTaskSpectraMC(const AliAnalysisTaskSpectraMC&);            // not implemented
		AliAnalysisTaskSpectraMC& operator=(const AliAnalysisTaskSpectraMC&); // not implemented

		ClassDef(AliAnalysisTaskSpectraMC, 1);    //Analysis task for high pt analysis
};

#endif



