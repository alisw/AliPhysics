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

		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(bool isMC) {fAnalysisMC = isMC;}
		virtual void  SetMCClosure(bool isMCclos) {fIsMCclosure = isMCclos;}
		virtual void  SetNcl(const int ncl){fNcl = ncl;}
		//	virtual void  SetPolarity(const float polarity){fPolarity = polarity;}
		virtual void  SetEtaCut(double etaCut){fEtaCut = etaCut;}
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
		virtual Double_t DeltaEta(Double_t eta, Double_t leta);
		void ProduceArrayTrksESD();
		void ProduceArrayV0ESD();
		short   GetPidCode(Int_t pdgCode) const;

		bool selectVertex2015pp(AliESDEvent* esd, Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
		bool IsGoodSPDvertexRes(const AliESDVertex* spdVertex = NULL);
		bool IsGoodZvertexPos(AliESDEvent *esd);
		bool PhiCut(const double& pt, double phi, const double& q, const float& mag, TF1* phiCutLow, TF1* phiCutHigh);
		float GetMaxDCApTDep( TF1 *fcut, Double_t pt );
		virtual void SetTrackCuts(AliAnalysisFilter* fTrackFilter);
		double EtaCalibration(const double &Eta);
		double EtaCalibrationEl(const double &Eta);
		bool TOFPID(AliESDtrack* track);

		static const Double_t fgkClight;   // Speed of light (cm/ps)

		AliESDEvent* fESD;                  //! ESD object
		AliEventCuts fEventCuts;
		AliMCEvent*  fMC;                   //! MC object
		AliStack*    fMCStack;              //! MC ESD stack
		TClonesArray* fMCArray;             //! MC array for AOD
		AliPIDResponse* fPIDResponse;       //! Pointer to PIDResponse
		AliESDtrackCuts* fGeometricalCut; 
		AliESDtrackCuts* fTrackFilterDaughters;
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

		int fNcl;
		double fEtaCut;
		///		const float fPolarity;
		bool fdEdxCalibrated;
		const double fDeDxMIPMin;
		const double fDeDxMIPMax;
		const double fdEdxHigh;
		const double fdEdxLow;
		TString  fPeriod;
		bool fSetTPConlyTrkCuts;
		bool fSelectHybridTracks;

		const double fLeadPtCutMin;
		const double fLeadPtCutMax;
		double fGenLeadPhi;
		double fGenLeadPt;
		int    fGenLeadIn;
		double fRecLeadEta;
		double fRecLeadPhi;
		double fRecLeadPt;
		int    fRecLeadIn;
		const double fPtMin;

		//
		// Output objects
		//

		TList*        fListOfObjects;     //! Output list of objects
		TH1F*         fEvents;            //! No of accepted events

		TH1F* hNchTSData;
		TH2F* hPhiTotal;
		TH2F* hPhiStandard;
		TH2F* hPhiHybrid1;
		TH2F* hPhiHybrid2;
		TH1F* hPhiLeading;
		TH2F* hDeltaPhiDeltaEta;

		TH2F* hPtVsP[4];
		TH2F* hnSigmaElectrons[4];
		TH1F* hPhiData[3];
		TH2F* hNchVsPtPosTPC[3][4];
		TH2F* hNchVsPtNegTPC[3][4];
		TH2F* hNchVsPtPosTOF[3][4];
		TH2F* hNchVsPtNegTOF[3][4];
		TH2F* hNchVsPPosTOF[3][4];
		TH2F* hNchVsPNegTOF[3][4];

		TH3F* hNchVsPtDataPosPionTPC[3][4];
		TH3F* hNchVsPtDataNegPionTPC[3][4];
		TH3F* hNchVsPtDataPosKaonTPC[3][4];
		TH3F* hNchVsPtDataNegKaonTPC[3][4];
		TH3F* hNchVsPtDataPosProtonTPC[3][4];
		TH3F* hNchVsPtDataNegProtonTPC[3][4];

		TH3F* hNchVsPtDataPosTOF[3][4];
		TH3F* hNchVsPtDataNegTOF[3][4];

		TF1* fEtaCalibrationPos;
		TF1* fEtaCalibrationNeg;
		TF1* fEtaCalibrationPosEl;
		TF1* fEtaCalibrationNegEl;
		TF1* fcutDCAxy;
		TF1* fcutLow;
		TF1* fcutHigh;

		// Histos rTPC

		TH2F* hMIPVsEta; 	
		TProfile* pMIPVsEta;
		TH2F* hPlateauVsEta; 	
		TProfile* pPlateauVsEta;
		TH2F* histEV0[4];
		TH2F* histPV0[4];
		TH2F* histPiV0[4];
		TH2F* hMIPVsEtaV0s;
		TProfile* pMIPVsEtaV0s;
		TH2F* hPhirTPC;
		TH2F* histPiTof[4];
		TH2F* hMIPVsPhi[4];
		TProfile* pMIPVsPhi[4];
		TH2F* hPlateauVsPhi[4];
		TProfile* pPlateauVsPhi[4];
		TH3F* hDeDxVsP[3][4];
		TH2F* hNchVsPrTPC[3][4];
		TH2F* hNchVsPtrTPC[3][4];


		AliAnalysisTaskSpectraRT(const AliAnalysisTaskSpectraRT&);            // not implemented
		AliAnalysisTaskSpectraRT& operator=(const AliAnalysisTaskSpectraRT&); // not implemented

		ClassDef(AliAnalysisTaskSpectraRT, 1);    //Analysis task for high pt analysis
};

#endif


