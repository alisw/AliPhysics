#ifndef AliAnalysisTaskPPvsMultINEL0_H
#define AliAnalysisTaskPPvsMultINEL0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */
//Authors: Sergio Iga ,sergio.iga@correo.nucleares.unam.mx

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TProfile.h>
#include <TTreeStream.h>
#include <TRandom.h>
#include <TObject.h>
#include <TH2.h>
#include <TH3.h>

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

#include <AliMCEvent.h>
#include <AliStack.h>


#include "AliPPVsMultUtils.h"
#include "AliAnalysisUtils.h"

#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliMultSelectionTask.h"

//#include "AliTransverseEventShape.h"

/*
class AliPPVsMultUtils;
class AliTransverseEventShape;
*/

class AliAnalysisTaskPPvsMultINEL0 : public AliAnalysisTaskSE {
	public:


		AliAnalysisTaskPPvsMultINEL0();
		AliAnalysisTaskPPvsMultINEL0(const char *name);
		virtual ~AliAnalysisTaskPPvsMultINEL0();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);

		Double_t GetVtxCut() { return fVtxCut; }   
		Double_t GetEtaCut() { return fEtaCut; }     


		virtual void  SetTrigger(UInt_t ktriggerInt) {ftrigBit = ktriggerInt;}  
		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
		virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void  SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;}


	private:
			
		virtual void AnalyzeESD(AliESDEvent* esd); 
		virtual void AnalyzeAOD(AliAODEvent* aod);
		virtual void AnalyzeMC(AliMCEvent* mc);
		virtual void AnalyzeESDforDCA(AliESDEvent* esdEvent);
		virtual Bool_t selectVertex2015pp(AliESDEvent *esd,Bool_t checkSPDres, Bool_t requireSPDandTrk,Bool_t checkProximity); 
		virtual Bool_t IsGoodSPDvertexRes(const AliESDVertex * spdVertex);
		virtual Bool_t isMCEventTrueINEL0(AliMCEvent* fMCEvent);
		virtual ULong64_t GetEventIdAsLong(AliVHeader* header);

		AliMCEvent *fMCEvent;
		AliStack   *fMCStack;

		AliESDEvent* fESD;                  //! ESD object
		AliAODEvent* fAOD;                  //! AOD object
		Bool_t fAnalysisMC;
		AliAnalysisFilter *fTrackFilterDCA;
		TString       fAnalysisType;        //  "ESD" or "AOD"
		UInt_t       ftrigBit;
		TRandom*      fRandom;              //! random number generator
		Bool_t        fPileUpRej;           // kTRUE is pile-up is rejected
		AliPPVsMultUtils *fPPVsMultUtils;
		AliMultSelection *fMultSelection;

		//
		// Cuts and options
		//

		Double_t     fVtxCut;             // Vtx cut on z position in cm
		Double_t     fEtaCut;             // Eta cut used to select particles
		//
		// Help variables
		//
		Short_t      fTriggeredEventMB;   // 1 = triggered, 0 = not trigged (MC only)
		Short_t      fVtxStatus;          // -1 = no vtx, 0 = outside cut, 1 = inside cut
		Float_t      fZvtx;               // z vertex
		Float_t      fZvtx_SPD;               // z vertex
		Int_t        fRun;                // run no
		ULong64_t    fEventId;            // unique event id
		Float_t fdcaxy;
		Float_t fdcaz;
		Bool_t fisPS;
		Bool_t fisTracklet;
		Bool_t fisMCvtxInZcut;

		//
		// Output objects
		//
		TList*        fListOfObjects;     //! Output list of objects
		TH1I*         fEvents;            //! No of accepted events
		TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
		TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts
		TH1F* 		fn1;

		Double_t fNchargedTrue;      //
		Double_t ftrackmult08;       //       
		Double_t fv0mpercentile;     //       
		Bool_t isINEL0Rec;
		Bool_t isINEL0True;
		//

		TH2F *refMultvsV0mINEL0;

		TH2F *refMultvsV0mTrigINEL0;
		TH2F *refMultvsV0mTrueINEL0;

		TH3F *ptvstrackletsvsdcaData;
		TH3F *ptvstrackletsvsdcaPrim;
		TH3F *ptvstrackletsvsdcaDecs;
		TH3F *ptvstrackletsvsdcaMatl;
		
		TH3F *ptvsv0mvsdcaData;
		TH3F *ptvsv0mvsdcaPrim;
		TH3F *ptvsv0mvsdcaDecs;
		TH3F *ptvsv0mvsdcaMatl;

		TH3F *ptvstrackletsvsdcacentralData;
		TH3F *ptvstrackletsvsdcacentralPrim;
		TH3F *ptvstrackletsvsdcacentralDecs;
		TH3F *ptvstrackletsvsdcacentralMatl;

		TH3F *ptvsv0mvsdcacentralData;
		TH3F *ptvsv0mvsdcacentralPrim;
		TH3F *ptvsv0mvsdcacentralDecs;
		TH3F *ptvsv0mvsdcacentralMatl;

		TH3F *effcomputationGen;
		TH3F *sigLossTrueINEL0;
		TH3F *sigLossTrigINEL0;
		TH2F *nchtruevsrefmult08;

		//for particle composition
		TH3F *effcomputationGenPi;
		TH3F *effcomputationGenK;
		TH3F *effcomputationGenP;
		TH3F *effcomputationGenSm;
		TH3F *effcomputationGenSp;
		TH3F *effcomputationGenO;
		TH3F *effcomputationGenXi;
		TH3F *effcomputationGenL;
		TH3F *effcomputationGenRest;

		TH3F *effcomputationRecPi;
		TH3F *effcomputationRecK;
		TH3F *effcomputationRecP;
		TH3F *effcomputationRecSm;
		TH3F *effcomputationRecSp;
		TH3F *effcomputationRecO;
		TH3F *effcomputationRecXi;
		TH3F *effcomputationRecL;
		TH3F *effcomputationRecRest;

		TH2F *fPS_MC;
		TH2F *fVtxPS_MC;
		TH2F *fPS;
		TH2F *fVtxPS;
		
		TH2F *refMultvsZvtx;
		TH1F *v0mPercentileQA;

		TH3F *secondaries[18];
		TH3F *primariesTrackFilter[18];
		TH2F* ptvsv0m[18];
		TH2F* ptvstracklets[18];
		TH3F *effcomputationRec[18];
		AliAnalysisFilter* fTrackFilter[18];// track filter





		AliAnalysisTaskPPvsMultINEL0(const AliAnalysisTaskPPvsMultINEL0&);            // not implemented
		AliAnalysisTaskPPvsMultINEL0& operator=(const AliAnalysisTaskPPvsMultINEL0&); // not implemented

		//TTree*        fTree;              //! Debug tree 

		ClassDef(AliAnalysisTaskPPvsMultINEL0, 1);    //Analysis task for high pt analysis 
};

#endif