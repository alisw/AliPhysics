#ifndef AliAnalysisTaskPPvsMult_H
#define AliAnalysisTaskPPvsMult_H
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



class AliAnalysisTaskPPvsMult : public AliAnalysisTaskSE
{
public:


		AliAnalysisTaskPPvsMult();
		AliAnalysisTaskPPvsMult(const char *name);
		virtual ~AliAnalysisTaskPPvsMult();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);

		Bool_t   GetAnalysisMC() { return fAnalysisMC; }   
//		Double_t GetVtxCut() { return fVtxCut; }   
		Double_t GetEtaCut() { return fEtaCut; }     
		//Double_t GetMinPt() { return fMinPt; }   
		//Int_t    GetTreeOption() { return fTreeOption; }  

//		virtual void  SetTrigger(UInt_t ktriggerInt) {ftrigBit = ktriggerInt;}
		virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
		virtual void  SetTrackFilterTPC(AliAnalysisFilter* trackF) {fTrackFilterTPC = trackF;}
		virtual void  SetTrackFilter2015PbPb(AliAnalysisFilter* trackF) {fTrackFilter2015PbPb = trackF;}
//		virtual void  SetCentralityEstimator(const char * centEst) {fCentEst = centEst;}
		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
		virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void  SetNcl(const Int_t ncl){fNcl = ncl;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
//		virtual void  SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;}   
		virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
		virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
//		virtual void  SetStoreMcIn(Bool_t value) {fStoreMcIn = value;}
		virtual void  SetAnalysisPbPb(Bool_t isanaPbPb) { fAnalysisPbPb = isanaPbPb; }
		virtual void  SetAnalysisTask(Bool_t PostCalib) { fdEdxCalibrated = PostCalib; }
		virtual void  SetAnalysisPID(Bool_t makePid) { fMakePid = makePid; }
		virtual void  SetAddLowPt(Bool_t addlowpt) { fLowPt = addlowpt; }
		virtual void  SetPeriod(Int_t isLHC16l) { fLHC16l = isLHC16l; }

	private:
		virtual Float_t GetVertex(const AliVEvent* event) const;
		virtual void AnalyzeESD(AliESDEvent* esd); 
//		virtual void AnalyzeAOD(AliAODEvent* aod); 
		virtual void ProduceArrayTrksESD(AliESDEvent* event, const Int_t cent);
		virtual void ProduceArrayV0ESD(AliESDEvent* event, const Int_t cent );
//		virtual void ProduceArrayTrksAOD(AliAODEvent* event);
//		virtual void ProduceArrayV0AOD(AliAODEvent* event);
		Short_t   GetPidCode(Int_t pdgCode) const;
		void      ProcessMCTruthESD( const Int_t cent );
//		void      ProcessMCTruthAOD(); 

		Short_t   GetPythiaEventProcessType(Int_t pythiaType);
		Short_t   GetDPMjetEventProcessType(Int_t dpmJetType);
		ULong64_t GetEventIdAsLong(AliVHeader* header) const;

		TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
		Int_t      FindPrimaryMotherLabel(AliStack* stack, Int_t label);

//		AliAODMCParticle* FindPrimaryMotherAOD(AliAODMCParticle* startParticle);

		TParticle* FindPrimaryMotherV0(AliStack* stack, Int_t label);
		Int_t      FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps);
		Bool_t selectVertex2015pp(AliESDEvent* esd, Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
                Bool_t IsGoodSPDvertexRes(const AliESDVertex* spdVertex = NULL);
                Bool_t IsGoodZvertexPos(AliESDEvent *esd);
		Bool_t PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh);
		Float_t GetMaxDCApTDep( TF1 *fcut, Double_t pt );
		Double_t EtaCalibrationNeg(const Int_t centrality, const Double_t Eta);
		Double_t EtaCalibrationPos(const Int_t centrality, const Double_t Eta);
		Double_t EtaCalibrationNegEl(const Int_t centrality, const Double_t Eta);
		Double_t EtaCalibrationPosEl(const Int_t centrality, const Double_t Eta);



//		AliAODMCParticle* FindPrimaryMotherAODV0(AliAODMCParticle* startParticle, Int_t& nSteps);



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
		AliAnalysisUtils* utils;
//		TString       fCentEst;             // V0A , V0M, 
		TString       fAnalysisType;        //  "ESD" or "AOD"
		Bool_t        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
		Bool_t        fAnalysisPbPb;        //  true you want to analyze PbPb data, false for pp
		//  Boolt_t       fAnalysisTask;
//		UInt_t        ftrigBit;
		TRandom*      fRandom;              //! random number generator
//		Bool_t        fPileUpRej;           // kTRUE is pile-up is rejected



		//
		// Cuts and options
		//

		Double_t     fVtxCut;             // Vtx cut on z position in cm
		Int_t        fNcl;                
		Double_t     fEtaCut;             // Eta cut used to select particles
		Int_t        cent; //minimum centrality
		Float_t      fMinCent; //minimum centrality
		Float_t      fMaxCent; //maximum centrality
                const Double_t fDeDxMIPMin;
                const Double_t fDeDxMIPMax;
                const Double_t fdEdxHigh;
                const Double_t fdEdxLow;
//		Bool_t       fStoreMcIn;          // Store MC input tracks

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
		//  TH1I*         fVtx;               //! Event vertex info
		TH1F*         fVtxMC;             //! Event vertex info for ALL MC events
		//  TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
		//  TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts
		//  TH1F* fn1;
		//  TH1F* hEvents;
		Bool_t       fdEdxCalibrated;
		Bool_t       fMakePid;
		Bool_t       fLowPt;
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


		TH1D *hPtAll[11];
		TH1D *hPtAllPos[11];
		TH1D *hPtAllNeg[11];
		TH1D *hPtPos[11][4];
		TH1D *hPtNeg[11][4];
		TH2D *hPtVsP[11][4];

		TH2D *hDeDxVsP[11][4];

		TH2D *hnSigmaPiPos[11][4];
		TH2D *hnSigmaKPos[11][4];
		TH2D *hnSigmaPPos[11][4];

		TH2D *hnSigmaPiNeg[11][4];
		TH2D *hnSigmaKNeg[11][4];
		TH2D *hnSigmaPNeg[11][4];

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

		TH2D* hPiondEdx[11];
                TH2D* hKaondEdx[11];
                TH2D* hProtondEdx[11];


		TH2D* hDCAxyVsPtPiNeg[11];
		TH2D* hDCAxyVsPtPiNegC[11];
		TH2D* hDCAxyVsPtKNeg[11];
		TH2D* hDCAxyVsPtKNegC[11];
		TH2D* hDCAxyVsPtPNeg[11];
		TH2D* hDCAxyVsPtPNegC[11];
		TH2D* hDCAxyVsPtPiPos[11];
		TH2D* hDCAxyVsPtPiPosC[11];
		TH2D* hDCAxyVsPtKPos[11];
		TH2D* hDCAxyVsPtKPosC[11];
		TH2D* hDCAxyVsPtPPos[11];
		TH2D* hDCAxyVsPtPPosC[11];

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
		TF1* fcutLow;
		TF1* fcutHigh;


		AliAnalysisTaskPPvsMult(const AliAnalysisTaskPPvsMult&);            // not implemented
		AliAnalysisTaskPPvsMult& operator=(const AliAnalysisTaskPPvsMult&); // not implemented

		//TTree*        fTree;              //! Debug tree 

		ClassDef(AliAnalysisTaskPPvsMult, 1);    //Analysis task for high pt analysis 
};

#endif
