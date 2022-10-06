#ifndef AliAnalysisTaskPPvsMult_H
#define AliAnalysisTaskPPvsMult_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id$ */


// ROOT includes
#include <TList.h>
#include <TH1.h>
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
#include <AliEventCuts.h>
#include "AliVTrack.h"



class AliAnalysisTaskPPvsMult : public AliAnalysisTaskSE
{
public:


		AliAnalysisTaskPPvsMult();
		AliAnalysisTaskPPvsMult(const char *name);
		virtual ~AliAnalysisTaskPPvsMult();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);

		Bool_t   GetAnalysisMC() { return fAnalysisMC; }   
		Double_t GetEtaCut() { return fEtaCut; }     

////		virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
		virtual void  SetTrackFilterTPC(AliAnalysisFilter* trackF) {fTrackFilterTPC = trackF;}
////		virtual void  SetTrackFilter2015PbPb(AliAnalysisFilter* trackF) {fTrackFilter2015PbPb = trackF;}
		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
		virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void  SetNcl(const Int_t ncl){fNcl = ncl;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
                virtual void  SetPeriod(const char* Period) { fPeriod = Period; }
//		virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
//		virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
//		virtual void  SetAnalysisPbPb(Bool_t isanaPbPb) { fAnalysisPbPb = isanaPbPb; }
		virtual void  SetAnalysisTask(Bool_t PostCalib) { fdEdxCalibrated = PostCalib; }
		virtual void  SetTrackCutsSystVars(const int TrackCutVar) { fTrackCuts = TrackCutVar; }
/////		virtual void  SetAnalysisPID(Bool_t makePid) { fMakePid = makePid; }
/////		virtual void  SetAddLowPt(Bool_t addlowpt) { fLowPt = addlowpt; }
		virtual void  SetPeriod(Int_t isLHC16l) { fLHC16l = isLHC16l; }

	private:
		virtual Float_t GetVertex(const AliVEvent* event) const;
		virtual void AnalyzeESD(AliESDEvent* esd); 
		virtual void ProduceArrayTrksESD(AliESDEvent* event, const Int_t cent);
		virtual void ProduceArrayV0ESD(AliESDEvent* event, const Int_t cent );
		Short_t   GetPidCode(Int_t pdgCode) const;
		void      ProcessMCTruthESD( const Int_t cent );

		Short_t   GetPythiaEventProcessType(Int_t pythiaType);
		Short_t   GetDPMjetEventProcessType(Int_t dpmJetType);
		ULong64_t GetEventIdAsLong(AliVHeader* header) const;

		TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
		Int_t      FindPrimaryMotherLabel(AliStack* stack, Int_t label);

		TParticle* FindPrimaryMotherV0(AliStack* stack, Int_t label);
		Int_t      FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps);
		Bool_t selectVertex2015pp(AliESDEvent* esd, Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
                Bool_t IsGoodSPDvertexRes(const AliESDVertex* spdVertex = NULL);
                Bool_t IsGoodZvertexPos(AliESDEvent *esd);
		Bool_t PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh);
		Float_t GetMaxDCApTDep( TF1 *fcut, Double_t pt );
		Double_t EtaCalibration(const double& Eta);
		Double_t EtaCalibrationEl(const double& Eta);
 		bool TOFPID(AliESDtrack* track);

		static const Double_t fgkClight;   // Speed of light (cm/ps)

		AliESDEvent* fESD;                  //! ESD object
		AliAODEvent* fAOD;                  //! AOD object
		AliEventCuts fEventCuts;
		AliMCEvent*  fMC;                   //! MC object
		AliStack*    fMCStack;              //! MC ESD stack
		TClonesArray* fMCArray;             //! MC array for AOD
		AliPIDResponse* fPIDResponse;       //! Pointer to PIDResponse
		AliESDtrackCuts* fTrackFilterGolden;    //  Track Filter, set 2010 with golden cuts
		AliAnalysisFilter* fTrackFilterTPC; // track filter for TPC only tracks
		AliAnalysisUtils* utils;
		TString       fAnalysisType;        //  "ESD" or "AOD"
		Bool_t        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
///		Bool_t        fAnalysisPbPb;        //  true you want to analyze PbPb data, false for pp
		TRandom*      fRandom;              //! random number generator

		//
		// Cuts and options
		//

		Double_t     fVtxCut;             // Vtx cut on z position in cm
		Int_t        fNcl;                
		Double_t     fEtaCut;             // Eta cut used to select particles
		Int_t        cent; //minimum centrality
//		Float_t      fMinCent; //minimum centrality
//		Float_t      fMaxCent; //maximum centrality
                const Double_t fDeDxMIPMin;
                const Double_t fDeDxMIPMax;
                const Double_t fdEdxHigh;
                const Double_t fdEdxLow;
		TString fPeriod;

		//
		// Help variables
		//
		Short_t      fMcProcessType;      // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
		Short_t      fTriggeredEventMB;   // 1 = triggered, 0 = not trigged (MC only)
		Short_t      fVtxStatus;          // -1 = no vtx, 0 = outside cut, 1 = inside cut
		Float_t      fZvtx;               // z vertex
		Float_t      fZvtxMC;             // z vertex MC (truth)
		Int_t        fRun;                // run no
		ULong64_t    fEventId;            // unique event id

		//
		// Output objects
		//
		TList*        fListOfObjects;     //! Output list of objects
		TH2F*         fEvents;            //! No of accepted events
		TH1F*         fV0M;            //! No of accepted events
		TH1F*         fVtxMC;             //! Event vertex info for ALL MC events
		Bool_t       fdEdxCalibrated;
		int fTrackCuts;
////		Bool_t       fMakePid;
////		Bool_t       fLowPt;
		Int_t  fLHC16l;
		TH1F* fcent;
		TH1F* fcentAfterPrimaries;
		TH1F* fcentAfterV0s;


		// Histograms for PreCalibration

		TH2D *hMIPVsEta[11];
		TProfile *pMIPVsEta[11];
		TH2D *hMIPVsEtaV0s[11];
		TProfile *pMIPVsEtaV0s[11];
		TH2D *hPlateauVsEta[11];
		TProfile *pPlateauVsEta[11];
		TH2D *hPhi[11];

		TH2D *hMIPVsV0M[4];
		TProfile *pMIPVsV0M[4];
		TH2D *hMIPVsNch[4];
		TProfile *pMIPVsNch[4];

		TH2D     *hMIPVsPhi[11][4];
		TProfile *pMIPVsPhi[11][4];
		TH2D     *hPlateauVsPhi[11][4];
		TProfile *pPlateauVsPhi[11][4];


		// Histograms for PostCalibration


		TH1D *hPt_TPC[11];
		TH1D *hPtpos_TPC[11];
		TH1D *hPtneg_TPC[11];
		TH1D *hPtpos_TPC_Eta[11][4];
		TH1D *hPtneg_TPC_Eta[11][4];
		TH2D *hPtVsP[11][4];

		TH2D *hDeDxVsP[11][4];
		TH2D *hnSigmaPiPos[11][4];
		TH2D *hnSigmaKPos[11][4];
		TH2D *hnSigmaPPos[11][4];
		TH2D *hnSigmaPiNeg[11][4];
		TH2D *hnSigmaKNeg[11][4];
		TH2D *hnSigmaPNeg[11][4];

                TH2D *hBetavsPneg[11][4];
                TH1D *hPtneg_TOF_Eta[11][4];
                TH1D *hPneg_TOF_Eta[11][4];
                TH1D *hPtneg_TOF[11];
                TH2D *hBetavsPpos[11][4];
                TH1D *hPtpos_TOF_Eta[11][4];
                TH1D *hPpos_TOF_Eta[11][4];
                TH1D *hPtpos_TOF[11];            

		TH2D* histPiV0[11][4];
		TH1D* histpPiV0[11][4];
		TH2D* histPV0[11][4];
		TH1D* histpPV0[11][4];
		TH2D* histPiTof[11][4];
		TH1D* histpPiTof[11][4];
		TH2D* histEV0[11][4];

		TH1D* hMcInNeg[3];
		TH1D* hMcInPos[3];
		TH1D* hMcOutNeg[3];
		TH1D* hMcOutPos[3];
		TH1D* hMcOutNegTOF[3];
		TH1D* hMcOutPosTOF[3];

		TH3F* hDCApTPrim[2][2];
		TH3F* hDCApTWDec[2][2];
		TH3F* hDCApTMate[2][2];

		TH3F* hDCApTPrim_TOF[2][2];
		TH3F* hDCApTWDec_TOF[2][2];
		TH3F* hDCApTMate_TOF[2][2];

		TH3F* hDCAxyVsPtPiNeg_TPC;
		TH3F* hDCAxyVsPtPNeg_TPC;
		TH3F* hDCAxyVsPtPiNeg_TOF;
		TH3F* hDCAxyVsPtPNeg_TOF;

		TH3F* hDCAxyVsPtPiPos_TPC;
		TH3F* hDCAxyVsPtPPos_TPC;
		TH3F* hDCAxyVsPtPiPos_TOF;
		TH3F* hDCAxyVsPtPPos_TOF;

		TF1* fEtaCalibrationPos;
		TF1* fEtaCalibrationNeg;
		TF1* felededxfitPos;
		TF1* felededxfitNeg;
		TF1* fcutDCAxy;
		TF1* fcutLow;
		TF1* fcutHigh;


		AliAnalysisTaskPPvsMult(const AliAnalysisTaskPPvsMult&);            // not implemented
		AliAnalysisTaskPPvsMult& operator=(const AliAnalysisTaskPPvsMult&); // not implemented

		//TTree*        fTree;              //! Debug tree 

		ClassDef(AliAnalysisTaskPPvsMult, 1);    //Analysis task for high pt analysis 
};

#endif
