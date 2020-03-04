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

		virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(bool isMC) {fAnalysisMC = isMC;}
		virtual void  SetMCClosure(bool isMCclos) {fIsMCclosure = isMCclos;}
		virtual void  SetNcl(const Int_t ncl){fNcl = ncl;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void  SetPeriod(const char* Period) { fPeriod = Period; }
		virtual void  SetMeanMultTSdata(const double MeanCh) { fMeanChT = MeanCh; }
		virtual void  SetMeanMultTSMCGen(const double MeanMultTSMCGen) { fMeanMultTSMCGen = MeanMultTSMCGen; }
		virtual void  SetMeanMultTSMCRec(const double MeanMultTSMCRec) { fMeanMultTSMCRec = MeanMultTSMCRec; }
		virtual void  SetTrackCutsType(bool isTPCOnlyTrkCuts) { fSetTPConlyTrkCuts = isTPCOnlyTrkCuts; }

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
		const Double_t fDeDxMIPMin;
		const Double_t fDeDxMIPMax;
		const Double_t fdEdxHigh;
		const Double_t fdEdxLow;
		TString  fPeriod;
		double fMeanChT;
		double fMeanMultTSMCGen;
		double fMeanMultTSMCRec;
		bool fSetTPConlyTrkCuts;

		const double fLeadPtCutMin;
		const double fLeadPtCutMax;
		double fGenLeadPhi;
		double fGenLeadPt;
		int    fGenLeadIn;
		double fRecLeadPhi;
		double fRecLeadPt;
		int    fRecLeadIn;
		double fPtMin;

		//
		// Output objects
		//

		TList*        fListOfObjects;     //! Output list of objects
		TH1F*         fEvents;            //! No of accepted events
		bool       fdEdxCalibrated;
		TH2D* hNchRecVsPtRecOut;
		TH2D* hNchGenVsPtGenIn[4];
		TH2D* hNchRecVsPtGenIn[4];
		TH2D* hNchGenVsPtGenPosIn[4];
		TH2D* hNchGenVsPtGenNegIn[4];
		TH2D* hNchGenVsPtRecIn[4];
		TH2D* hNchGenVsPtRecPosIn[4];
		TH2D* hNchGenVsPtRecNegIn[4];
		TH2D* hNchGenVsPtRecInTOF[4];
		TH2D* hNchGenVsPtRecPosInTOF[4];
		TH2D* hNchGenVsPtRecNegInTOF[4];
		TH2D* hNchRecVsPtGenOut;
		TH1D* hPtPriGen;
		TH1D* hPtRec;
		TH1D* hPtPriRec;
		TH1D* hPtSecRec;
		TH2D* fPtLVsNchGen; 
		TH2D* fPtLVsNchRec; 

		TH2D *hPtVsP[4];

		TH1D* hPhiGen[3];
		TH1D* hPhiRec[3];

		TH1D* hMultTSGen;
		TH1D* hMultTSRec;
		TH1D* hNchTSGen;
		TH1D* hNchTSRec;
		TH1D* hNchTSGen_1;
		TH1D* hNchTSContamination;
		TH1D* hNchTSRecAll;
		TH2D* hNchResponse;
		TH2D* hPtResponsePID[4];
		TH2D* hNchGenVsPtGenPID[3][4];
		TH2D* hNchGenVsPtRec[3][4];

		TH1D* hNchTSGenTest;
		TH1D* hNchTSRecTest;

		TH1D* hNchTSData;
		TH1D* hRTData;
		TH2D* hPtLVsRT;
		TH1D* hPhiData[3];

		TH2D* hNchVsPtDataTPC[3][4];
		TH2D* hNchVsPtDataPosTPC[3][4];
		TH2D* hNchVsPtDataNegTPC[3][4];

		TH2D* hNchVsPtDataTOF[3][4];
		TH2D* hNchVsPtDataPosTOF[3][4];
		TH2D* hNchVsPtDataNegTOF[3][4];

		TF1* fEtaCalibration;
		TF1* fEtaCalibrationEl;
		TF1* fcutDCAxy;
		TF1* fcutLow;
		TF1* fcutHigh;

		AliAnalysisTaskSpectraRT(const AliAnalysisTaskSpectraRT&);            // not implemented
		AliAnalysisTaskSpectraRT& operator=(const AliAnalysisTaskSpectraRT&); // not implemented

		ClassDef(AliAnalysisTaskSpectraRT, 1);    //Analysis task for high pt analysis
};

#endif

