#ifndef AliAnalysisTaskSpectraRT_H
#define AliAnalysisTaskSpectraRT_H
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


class AliAnalysisTaskSpectraRT : public AliAnalysisTaskSE
{
	public:


		AliAnalysisTaskSpectraRT();
		AliAnalysisTaskSpectraRT(const char *name);
		virtual ~AliAnalysisTaskSpectraRT();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);

		Bool_t   GetAnalysisMC() { return fAnalysisMC; }
		Double_t GetEtaCut() { return fEtaCut; }

////		virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(bool isMC) {fAnalysisMC = isMC;}
		virtual void  SetMCClosure(bool isMCclos) {fIsMCclosure = isMCclos;}
		virtual void  SetNcl(const Int_t ncl){fNcl = ncl;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void  SetPeriod(const char* Period) { fPeriod = Period; }
		virtual void  SetdEdxCalibration(bool isCalibrated) { fdEdxCalibrated = isCalibrated; }
		virtual void  SetTrackCutsType(bool isTPCOnlyTrkCuts) { fSetTPConlyTrkCuts = isTPCOnlyTrkCuts; }
		virtual void  SetHybridTracks(bool isSelectHybrid) { fSelectHybridTracks = isSelectHybrid; }
		AliESDtrack*  SetHybridTrackCuts(AliESDtrack *track, const bool fill1, const bool fill2, const bool fill3);

	private:


		AliESDtrack* GetLeadingTrack();
		void GetLeadingObject(Bool_t isMC);
		void GetMultiplicityDistributions();
		void GetDetectorResponse();
		void GetMCCorrections();
		virtual Double_t DeltaPhi(Double_t phi, Double_t lphi,
				Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
		void ProduceArrayTrksESD();
		short   GetPidCode(Int_t pdgCode) const;

		bool selectVertex2015pp(AliESDEvent* esd, Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
		bool IsGoodSPDvertexRes(const AliESDVertex* spdVertex = NULL);
		bool IsGoodZvertexPos(AliESDEvent *esd);
		bool PhiCut(const double& pt, double phi, const double& q, const float& mag, TF1* phiCutLow, TF1* phiCutHigh);
		float GetMaxDCApTDep( TF1 *fcut, Double_t pt );
		virtual void SetTrackCuts(AliAnalysisFilter* fTrackFilter);
		double EtaCalibration(const int &centrality, const double &Eta);
		double EtaCalibrationEl(const int &centrality, const double &Eta);
		int GetIndex();
		bool TOFPID(AliESDtrack* track);

		static const Double_t fgkClight;   // Speed of light (cm/ps)

		AliESDEvent* fESD;                  //! ESD object
		AliEventCuts fEventCuts;
		AliMCEvent*  fMC;                   //! MC object
		AliStack*    fMCStack;              //! MC ESD stack
		TClonesArray* fMCArray;             //! MC array for AOD
		AliPIDResponse* fPIDResponse;       //! Pointer to PIDResponse
		AliAnalysisFilter* fTrackFilterGolden;    //  Track Filter, set 2010 with golden cuts
		AliAnalysisFilter* fTrackFilter;
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

		int        fNcl;
		double     fEtaCut;             // Eta cut used to select particles
		bool fdEdxCalibrated;
		const Double_t fDeDxMIPMin;
		const Double_t fDeDxMIPMax;
		const Double_t fdEdxHigh;
		const Double_t fdEdxLow;
		TString  fPeriod;
		bool fSetTPConlyTrkCuts;
		bool fSelectHybridTracks;

		const double fLeadPtCutMin;
		const double fLeadPtCutMax;
		double fGenLeadPhi;
		double fGenLeadPt;
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

		TH2D* hNchGenVsPtGenIn[4][4];
		TH2D* hNchGenVsPtRecIn[4][4];
		TH2D* hNchRecVsPtGenIn[4][4];
		TH2D* hNchRecVsPtRecIn[4][4];
		TH2D* hNchGenVsPtGenPosIn[4][4];
		TH2D* hNchGenVsPtGenNegIn[4][4];
		TH2D* hNchGenVsPtRecPosIn[4][4];
		TH2D* hNchGenVsPtRecNegIn[4][4];
		TH2D* hNchGenVsPtRecInTOF[4][4];
		TH2D* hNchGenVsPtRecPosInTOF[4][4];
		TH2D* hNchGenVsPtRecNegInTOF[4][4];
		TH1D* hPtPriGen;
		TH1D* hPtRec;
		TH1D* hSecPtRec;
		TH1D* hSecPtGen;
		TH2D* fPtLVsNchGen; 
		TH2D* fPtLVsNchRec; 


		TH1D* hPhiGen[4];
		TH1D* hPhiRec[4];

		TH1D* hMultTSGen;
		TH1D* hMultTSRec;
		TH1D* hNchTSGen;
		TH1D* hNchTSRec;
		TH1D* hNchTSGenGTZ;
		TH1D* hNchTSCont;
////		TH1D* hNchTSRecPri;
		TH1D* hNchTSRecGTZ;
		TH2D* hNchResponse;
////		TH3D* hNchRMvsPt;
		TH1D* hPtTS;
		TH2D* hPtResponsePID[4];
		TH2D* hNchGenVsPtGenPID[4][4];
		TH2D* hNchGenVsPtRec[4][4];

		TH1D* hNchTSData;
		TH2D* hPtLVsRT;
		TH1D* hPhiData[4];
		TH2D* hPhiTotal;
		TH2D* hPhiStandard;
		TH2D* hPhiHybrid1;

		TH3D* hNchVsPtDataTPC[4][4];
		TH3D* hNchVsPtDataPosPionTPC[4];
		TH3D* hNchVsPtDataNegPionTPC[4];
		TH3D* hNchVsPtDataPosKaonTPC[4];
		TH3D* hNchVsPtDataNegKaonTPC[4];
		TH3D* hNchVsPtDataPosProtonTPC[4];
		TH3D* hNchVsPtDataNegProtonTPC[4];

		TH3D* hNchVsPtDataPosTOF[4][4];
		TH3D* hNchVsPtDataNegTOF[4][4];

		TF1* fEtaCalibration;
		TF1* fEtaCalibrationEl;
		TF1* fcutDCAxy;
		TF1* fcutLow;
		TF1* fcutHigh;

		// Histos rTPC
		
		TH2D* hMIPVsEta; 	
		TProfile* pMIPVsEta;
		TH2D* hPlateauVsEta; 	
		TProfile* pPlateauVsEta;
		TH2D* histPiTof[4];
		TH2D* hMIPVsPhi[4];
		TProfile* pMIPVsPhi[4];
		TH2D* hPlateauVsPhi[4];
		TProfile* pPlateauVsPhi[4];
		TH3D* hDeDxVsP[4][4];
		TH2D* hPtVsP[4];

		AliAnalysisTaskSpectraRT(const AliAnalysisTaskSpectraRT&);            // not implemented
		AliAnalysisTaskSpectraRT& operator=(const AliAnalysisTaskSpectraRT&); // not implemented

		ClassDef(AliAnalysisTaskSpectraRT, 1);    //Analysis task for high pt analysis
};

#endif

