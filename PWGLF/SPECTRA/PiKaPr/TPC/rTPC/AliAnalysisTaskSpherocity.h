#ifndef AliAnalysisTaskSpherocity_H
#define AliAnalysisTaskSpherocity_H
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
#include <AliMCEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliAODMCParticle.h> 
#include <AliESDtrackCuts.h>
#include <AliPIDResponse.h>
#include "AliTPCPIDResponse.h"
//#include <AliSpherocityUtils.h>
#include <AliEventCuts.h>
#include "AliVTrack.h"
#include <vector>
using namespace std;



class AliAnalysisTaskSpherocity : public AliAnalysisTaskSE
{
	public:


		AliAnalysisTaskSpherocity();
		AliAnalysisTaskSpherocity(const char *name);
		virtual ~AliAnalysisTaskSpherocity();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);

		Bool_t   GetAnalysisMC() { return fAnalysisMC; }   
		Double_t GetEtaCut() { return fEtaCut; }     

		virtual void  SetTrackCutsSpherocity(AliAnalysisFilter* fTrackFilter);
		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
		virtual void  SetNcl(const int ncl){fNcl = ncl;}
		virtual void  SetEtaCut(double etaCut){ fEtaCut = etaCut; }
		virtual void  SetMinMult(int multCut){ fMinMult = multCut; }
		virtual void  SetAnalysisTask(bool PostCalib) { fdEdxCalibrated = PostCalib; }
		virtual void  SetTrackCutsSystVars(const int TrackCutVar) { fTrackCuts = TrackCutVar; }
		virtual void  SetPeriod(const char* period) { fPeriod = period; }
		virtual void  SetEstimator(const Bool_t isV0M) { fisV0Mestimator = isV0M; }
		virtual void  SetJettyCutOff(float JettyCutOff) { fJettyCutOff = JettyCutOff; }
		virtual void  SetJettyCutOff_0(float JettyCutOff_0) { fJettyCutOff_0 = JettyCutOff_0; }
		virtual void  SetJettyCutOff_1(float JettyCutOff_1) { fJettyCutOff_1 = JettyCutOff_1; }
		virtual void  SetJettyCutOff_2(float JettyCutOff_2) { fJettyCutOff_2 = JettyCutOff_2; }
		virtual void  SetIsotrCutOff(float IsotrCutOff) { fIsotrCutOff = IsotrCutOff; }
		virtual void  SetIsotrCutOff_0(float IsotrCutOff_0) { fIsotrCutOff_0 = IsotrCutOff_0; }
		virtual void  SetIsotrCutOff_1(float IsotrCutOff_1) { fIsotrCutOff_1 = IsotrCutOff_1; }
		virtual void  SetIsotrCutOff_2(float IsotrCutOff_2) { fIsotrCutOff_2 = IsotrCutOff_2; }

	private:

		virtual void ProduceArrayTrksESD(const float& so);
		///		virtual void ProduceArrayV0ESD(AliESDEvent* event, const Int_t cent, const Int_t sperocity );
		void AnalyseSimDataV0M(const float& ,const float& );
		void AnalyseSimDataCL1(const float& ,const float& );
		Int_t   GetMultiplicityParticles(Double_t etaCut);
		virtual double DeltaPhi(double ,double , double rangeMin = -TMath::Pi()/2, double rangeMax = 3*TMath::Pi()/2 );
		Short_t GetPidCode(Int_t pdgCode) const;
		void ProcessMCTruthV0M(const float& , const float& , const float& );
		void ProcessMCTruthCL1(const float& , const float& , const float& );

//		TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
//		Int_t      FindPrimaryMotherLabel(AliStack* stack, Int_t label);


//		TParticle* FindPrimaryMotherV0(AliStack* stack, Int_t label);
//		int      FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps);
		bool selectVertex2015pp(AliESDEvent* esd, Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
		bool IsGoodSPDvertexRes(const AliESDVertex* spdVertex = NULL);
		bool IsGoodZvertexPos(AliESDEvent *esd);
		bool PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh);
		float GetMaxDCApTDep(TF1* fcut, Double_t pt );
		float GetEventShapeTrue( AliStack* , TH1D*, TH1D* );
		float GetSpherocityMC( TH1D* , TH1D* );
		int ReadMC( vector<Float_t>& ,  vector<Float_t>& , vector<Float_t>& , TH1D* , TH1D* );
		float GetSpherocity( TH1D * hphi, TH1D *heta );
		int ReadESDEvent( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray, TH1D * hphi, TH1D *heta );
		int ReadESDEvent();
		float AnalyseGetSpherocity( const vector<Float_t> &pt, const vector<Float_t> &eta, const vector<Float_t> &phi );
		bool TOFPID(AliESDtrack* track);
		double EtaCalibration(const double& eta);

		void GetNT();
		void GetLeadingObject();
		AliESDtrack* SetHybridTrackCuts(AliESDtrack * , const bool , const bool );

		static const Double_t fgkClight;   // Speed of light (cm/ps)

		AliESDEvent* fESD;                  //! ESD object
		AliAODEvent* fAOD;                  //! AOD object
		AliEventCuts fEventCuts;
		AliMCEvent*  fMC;                   //! MC object
		AliStack*    fMCStack;              //! MC ESD stack
		TClonesArray* fMCArray;             //! MC array for AOD
		AliPIDResponse* fPIDResponse;       //! Pointer to PIDResponse
		AliESDtrackCuts* fTrackFilterGolden;    //  Track Filter, set 2010 with golden cuts
		AliAnalysisFilter* fTrackFilter;
		AliAnalysisUtils* utils;
		TString       fAnalysisType;        //  "ESD" or "AOD"
		Bool_t        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
		Bool_t        fisV0Mestimator; 
		TRandom*      fRandom;              //! random number generator

		//
		// Cuts and options
		//

		int fNcl;                
		double fEtaCut;
		double fPtMinCut;
		double fPtMaxCut;
		const Double_t fDeDxMIPMin;
		const Double_t fDeDxMIPMax;
		int fdEdxHigh;
		int fdEdxLow;
		float fJettyCutOff;
		float fIsotrCutOff;
		float fJettyCutOff_0;
		float fJettyCutOff_1;
		float fJettyCutOff_2;
		float fIsotrCutOff_0;
		float fIsotrCutOff_1;
		float fIsotrCutOff_2;
		int fMinMult;
		int fNrec;
		float fSizeStep;

		//
		// Help variables
		//
		Short_t      fMcProcessType;      // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
		Short_t      fTriggeredEventMB;   // 1 = triggered, 0 = not trigged (MC only)
		Short_t      fVtxStatus;          // -1 = no vtx, 0 = outside cut, 1 = inside cut
		Float_t      fZvtx;               // z vertex
		Float_t      fZvtxMC;             // z vertex MC (truth)

		//
		// Output objects
		//
		TList*        fListOfObjects;     //! Output list of objects
		TH2F*         fEvents;            //! No of accepted events
		bool       fdEdxCalibrated;
		int fTrackCuts;
		TString fPeriod;

		// Histograms for Spherocity

		TH1D *hphiso;
		TH1D *hetaso;
		TH1D *hTruthPhiSo;
		TH1D *hTruthEtaSo;
		TH2F *hSOrvsV0M;
		TH2F *hSOrvsTrks;
		TH3F* hMultPercvsNch;
		TH3F *hSOtvsSOrV0M;
		TH3F *hSOtvsSOrCL1;
		TH3F *hPionSimvsSOV0M;
		TH3F *hPionSimvsSOV0M2;
		TH3F *hPionGenvsSOV0M;
		TH3F *hKaonSimvsSOV0M;
		TH3F *hKaonSimvsSOV0M2;
		TH3F *hKaonGenvsSOV0M;
		TH3F *hProtonSimvsSOV0M;
		TH3F *hProtonSimvsSOV0M2;
		TH3F *hProtonGenvsSOV0M;

		TH3F *hPionSimvsSOCL1;
		TH3F *hPionSimvsSOCL12;
		TH3F *hPionGenvsSOCL1;
		TH3F *hKaonSimvsSOCL1;
		TH3F *hKaonSimvsSOCL12;
		TH3F *hKaonGenvsSOCL1;
		TH3F *hProtonSimvsSOCL1;
		TH3F *hProtonSimvsSOCL12;
		TH3F *hProtonGenvsSOCL1;

		TH3F* hPhiSimV0M;
		TH3F* hPhiGenV0M;
		TH3F* hPhiSimCL1;
		TH3F* hPhiGenCL1;

		// Histograms for PreCalibration

		//TH2D *hMIPVsEta[11][3];
		//TProfile *pMIPVsEta[11][3];
		//TH2D *hMIPVsEtaV0s[11][3];
		//TProfile *pMIPVsEtaV0s[11][3];
		//TH2D *hPlateauVsEta[11][3];
		//TProfile *pPlateauVsEta[11][3];
		//TH2D *hPhi[11];

		//		TH2D *hMIPVsV0M[4];
		//		TProfile *pMIPVsV0M[4];
		//		TH2D *hMIPVsNch[4];
		//		TProfile *pMIPVsNch[4];

		//		TH2D     *hMIPVsPhi[11][4][3];
		//		TProfile *pMIPVsPhi[11][4][3];
		//		TH2D     *hPlateauVsPhi[11][4][3];
		//		TProfile *pPlateauVsPhi[11][4][3];


		// Histograms for PostCalibration


		TH2F *hPtAll[2][4];
		TH2D *hPtVsP[4];

		TH3F *hDeDxVsP[2][4][4];

		TH3F *hnSigmaPiPos[2][4][4];
		TH3F *hnSigmaKPos[2][4][4];
		TH3F *hnSigmaPPos[2][4][4];
		TH3F *hnSigmaPiNeg[2][4][4];
		TH3F *hnSigmaKNeg[2][4][4];
		TH3F *hnSigmaPNeg[2][4][4];

		TH3F *hBetavsPneg[2][4][4];
		TH3F *hBetavsPpos[2][4][4];

		TH2F *hPtneg_TPC_Eta[2][4][4];
		TH2F *hPtpos_TPC_Eta[2][4][4];
		TH2F *hPtneg_TPC[2][4];
		TH2F *hPtpos_TPC[2][4];
		TH2F *hPtneg_TOF_Eta[2][4][4];
		TH2F *hPtpos_TOF_Eta[2][4][4];
		TH2F *hPneg_TOF_Eta[2][4][4];
		TH2F *hPpos_TOF_Eta[2][4][4];
		TH2F *hPtneg_TOF[2][4];
		TH2F *hPtpos_TOF[2][4];

		TH2F* hPtAllSoInt;
		TH2F* hPtpos_TPCSoInt;
		TH2F* hPtneg_TPCSoInt;
		TH2F* hPtpos_TOFSoInt;
		TH2F* hPtneg_TOFSoInt;

		TH3F* hDeDxVsPSoInt[4];
		TH3F* hnSigmaPiPosSoInt[4];
		TH3F* hnSigmaKPosSoInt[4];
		TH3F* hnSigmaPPosSoInt[4];
		TH3F* hnSigmaPiNegSoInt[4];
		TH3F* hnSigmaKNegSoInt[4];
		TH3F* hnSigmaPNegSoInt[4];

		TH3F* hBetavsPposSoInt[4];
		TH3F* hBetavsPnegSoInt[4];

		TH2F* hPtpos_TPC_EtaSoInt[4];
		TH2F* hPtneg_TPC_EtaSoInt[4];
		TH2F* hPtpos_TOF_EtaSoInt[4];
		TH2F* hPtneg_TOF_EtaSoInt[4];
		TH2F* hPpos_TOF_EtaSoInt[4];
		TH2F* hPneg_TOF_EtaSoInt[4];

		TF1* fEtaCalibrationNeg;
		TF1* fEtaCalibrationPos;
		TF1* fcutDCAxy;
		TF1* fcutLow;
		TF1* fcutHigh;

                double fRecLeadPhi;
                double fRecLeadPt;
                int fRecLeadIn;
                int fNT;
                TH2F* hPhiStandard;
                TH2F* hPhiHybrid1;
                TH2F* hPhiTotal;
                AliESDtrackCuts* fGeometricalCut;
                AliESDtrackCuts* fHybridTrackCuts1;
                TH3F* hV0MvsSOvsNT;
                TH3F* hCL1vsSOvsNT;
                TH3F* hV0MvsSOvsNoNT;
                TH3F* hCL1vsSOvsNoNT;
                TH3F* hV0MvsNoSOvsNT;
                TH3F* hCL1vsNoSOvsNT;
                TH3F* hV0MvsNoSOvsNoNT;
                TH3F* hCL1vsNoSOvsNoNT;
		TH1F* hPhiToward;
		TH1F* hPhiAway;
		TH1F* hPhiTransverse;

		AliAnalysisTaskSpherocity(const AliAnalysisTaskSpherocity&);            // not implemented
		AliAnalysisTaskSpherocity& operator=(const AliAnalysisTaskSpherocity&); // not implemented
		ClassDef(AliAnalysisTaskSpherocity, 1);    //Analysis task for high pt analysis 
};

#endif


