#ifndef ALIANALYSISTASKRTPCPPMULTIPLICITY_H
#define ALIANALYSISTASKRTPCPPMULTIPLICITY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id$ */


// ROOT includes
#include <TList.h>
#include <TH1.h>
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



class AliAnalysisTaskrTPCPPMultiplicty : public AliAnalysisTaskSE {
	public:


		AliAnalysisTaskrTPCPPMultiplicty();
		AliAnalysisTaskrTPCPPMultiplicty(const char *name);
		virtual ~AliAnalysisTaskrTPCPPMultiplicty();






		//AliAnalysisTaskQAHighPtDeDx(const char *name="<default name>");
		//virtual ~AliAnalysisTaskQAHighPtDeDx() { /*if (fOutputList) delete fOutputList;*/}//;

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);

		Bool_t   GetAnalysisMC() { return fAnalysisMC; }   
		Double_t GetVtxCut() { return fVtxCut; }   
		Double_t GetEtaCut() { return fEtaCut; }     
		//Double_t GetMinPt() { return fMinPt; }   
		//Int_t    GetTreeOption() { return fTreeOption; }  

		virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
		virtual void  SetTrackFilterTPC(AliAnalysisFilter* trackF) {fTrackFilterTPC = trackF;}
		virtual void  SetTrackFilter2015PbPb(AliAnalysisFilter* trackF) {fTrackFilter2015PbPb = trackF;}
		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
		virtual void  SetAnalysislowpT(Bool_t islowpT) {flowpT = islowpT;}
		virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void  SetNcl(const Int_t ncl){fNcl = ncl;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
		virtual void  SetMinMult(Float_t minvalc) {fMinMult = minvalc;}
		virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
		virtual void  SetMaxMult(Float_t maxvalc) {fMaxMult = maxvalc;}
		virtual void  SetStoreMcIn(Bool_t value) {fStoreMcIn = value;}
		virtual void  SetAnalysisPbPb(Bool_t isanaPbPb) {fAnalysisPbPb = isanaPbPb;}
		virtual void  SetAnalysisTask(Bool_t PostCalib) {fdEdxCalibrated = PostCalib;}
		virtual void  SetAnalysisPID(Bool_t makePid) {fMakePid = makePid;}

	private:
		virtual Float_t GetVertex(const AliVEvent* event) const;
		virtual void AnalyzeESD(AliESDEvent* esd); 
		virtual void AnalyzeAOD(AliAODEvent* aod); 
		virtual void ProduceArrayTrksESD(AliESDEvent* event, const Int_t cent);
		virtual void ProduceArrayV0ESD(AliESDEvent* event, const Int_t cent );
		virtual void ProduceArrayTrksAOD(AliAODEvent* event);
		virtual void ProduceArrayV0AOD(AliAODEvent* event);
		Short_t   GetPidCode(Int_t pdgCode) const;
		void      ProcessMCTruthESD( const Int_t cent );
		void      ProcessMCTruthAOD(); 

		Short_t   GetPythiaEventProcessType(Int_t pythiaType);
		Short_t   GetDPMjetEventProcessType(Int_t dpmJetType);
		ULong64_t GetEventIdAsLong(AliVHeader* header) const;

		TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
		Int_t      FindPrimaryMotherLabel(AliStack* stack, Int_t label);

		AliAODMCParticle* FindPrimaryMotherAOD(AliAODMCParticle* startParticle);

		TParticle* FindPrimaryMotherV0(AliStack* stack, Int_t label);
		Int_t      FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps);
		Bool_t PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh);
		Float_t GetMaxDCApTDep( TF1 *fcut, Double_t pt );
		Double_t EtaCalibrationNeg(const Int_t centrality, const Double_t Eta);
		Double_t EtaCalibrationPos(const Int_t centrality, const Double_t Eta);
		Double_t EtaCalibrationNegEl(const Int_t centrality, const Double_t Eta);
		Double_t EtaCalibrationPosEl(const Int_t centrality, const Double_t Eta);



		AliAODMCParticle* FindPrimaryMotherAODV0(AliAODMCParticle* startParticle, Int_t& nSteps);



		static const Double_t fgkClight;   // Speed of light (cm/ps)
		//  static const Int_t nCent=10;   // Speed of light (cm/ps)

		AliESDEvent* fESD;                  //! ESD object
		AliAODEvent* fAOD;                  //! AOD object
		AliEventCuts fEventCuts;
		AliMCEvent*  fMC;                   //! MC object
		AliStack*    fMCStack;              //! MC ESD stack
		TClonesArray* fMCArray;             //! MC array for AOD
		AliPIDResponse* fPIDResponse;       //! Pointer to PIDResponse
		AliAnalysisFilter* fTrackFilter2015PbPb;    //  Track Filter, set 2010 with golden cuts
		AliAnalysisFilter* fTrackFilterGolden;    //  Track Filter, set 2010 with golden cuts
		AliAnalysisFilter* fTrackFilterTPC; // track filter for TPC only tracks
		TString       fAnalysisType;        //  "ESD" or "AOD"
		Bool_t        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
		Bool_t        flowpT;               //  true to include analysis of low pT
		Bool_t        fAnalysisPbPb;        //  true you want to analyze PbPb data, false for pp
		TRandom*      fRandom;              //! random number generator



		//
		// Cuts and options
		//

		Double_t     fVtxCut;          // Vtx cut on z position in cm
		Int_t  	     fNcl;             // Ncl cut: Default 70
		Double_t     fEtaCut;          // Eta cut used to select particles
		Int_t        cent; 	       
		Float_t      fMinCent;         //Minimum centrality accepted
		Float_t      fMinMult;         //Minimum multiplicity accepted
		Float_t      fMaxCent;         //Maximum centrality accepted
		Float_t      fMaxMult;         //Maximum multiplicity accepted
		Bool_t       fStoreMcIn;       // Store MC input tracks

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
		//  TH1I*         fEvents;            //! No of accepted events
		//  TH1I*         fVtx;               //! Event vertex info
		TH1F*         fVtxMC;             //! Event vertex info for ALL MC events
		//  TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
		//  TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts
		//  TH1F* fn1;
		//  TH1F* hEvents;
		Bool_t       fdEdxCalibrated;
		Bool_t       fMakePid;
		TH1F* fcent;
		TH1F* fcentAfterPrimaries;
		TH1F* fcentAfterV0s;
		TH1I* fHistEventCounter;


		// Histograms for PreCalibration

		TH2D *hMIPVsEta[11];
		TProfile *pMIPVsEta[11];
		TH2D *hMIPVsEtaV0s[11];
		TProfile *pMIPVsEtaV0s[11];
		TH2D *hPlateauVsEta[11];
		TProfile *pPlateauVsEta[11];
		TH2D *hPhi[11];

		TH2D     *hMIPVsNch[11][4];
		TProfile *pMIPVsNch[11][4];

		TH2D     *hMIPVsPhi[11][4];
		TProfile *pMIPVsPhi[11][4];
		TH2D     *hPlateauVsPhi[11][4];
		TProfile *pPlateauVsPhi[11][4];


		// Histograms for PostCalibration


		TH1D *hPtAll[11];
		TH2D *hPtVsP[11][4];

		TH2D *hDeDxVsP[11][4];

		TH2D* histPiV0[11][4];
		TH1D* histpPiV0[11][4];

		TH2D* histPV0[11][4];
		TH1D* histpPV0[11][4];

		TH2D* histPiTof[11][4];
		TH1D* histpPiTof[11][4];

		TH2D* histEV0[11][4];

		TH1D* hMcIn[11][7];
		TH1D* hMcOut[11][7];
		TH1D* hMcInNeg[11][7];
		TH1D* hMcInPos[11][7];
		TH1D* hMcOutNeg[11][7];
		TH1D* hMcOutPos[11][7];

		//    [Cent][Pid][Charge: 0:neutral 1:Neg 2:Pos]
		TH2D* hDCApTPrim[10][7][3];
		TH2D* hDCApTWDec[10][7][3];
		TH2D* hDCApTMate[10][7][3];

		TH2D* hDCApTPrim2[10][7][3];
		TH2D* hDCApTWDec2[10][7][3];
		TH2D* hDCApTMate2[10][7][3];

		TF1* fEtaCalibrationNeg;
		TF1* fEtaCalibration;
		TF1* felededxfitPos;
		TF1* felededxfitNeg;
		TF1* fcutDCAxy;


		AliAnalysisTaskrTPCPPMultiplicty(const AliAnalysisTaskrTPCPPMultiplicty&);            // not implemented
		AliAnalysisTaskrTPCPPMultiplicty& operator=(const AliAnalysisTaskrTPCPPMultiplicty&); // not implemented

		//TTree*        fTree;              //! Debug tree 

		ClassDef(AliAnalysisTaskrTPCPPMultiplicty, 1);    //Analysis task for high pt analysis 
};

#endif
