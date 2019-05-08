#ifndef AliAnalysisTaskPPvsUE_H
#define AliAnalysisTaskPPvsUE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */
//Authors: Sergio Iga ,sergio.iga@correo.nucleares.unam.mx
//         Valentina Zaccolo, valentina.zaccolo@cern.ch

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
#include <AliInputEventHandler.h>
#include <AliMCEventHandler.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliVEvent.h>
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


class AliPPVsMultUtils;
class AliTransverseEventShape;


class AliAnalysisTaskPPvsUE : public AliAnalysisTaskSE {
	public:


		AliAnalysisTaskPPvsUE();
		AliAnalysisTaskPPvsUE(const char *name);
		virtual ~AliAnalysisTaskPPvsUE();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);

		Double_t GetVtxCut() { return fVtxCut; }   
		Double_t GetEtaCut() { return fEtaCut; }     


		virtual void  SetTrigger(UInt_t ktriggerInt = AliVEvent::kINT7) {ftrigBit = ktriggerInt;}  
		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
		virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void  SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;}
		virtual void  SetAveMultiInTrans(Double_t value) {fAveMultiInTrans = value;}
		
		TObjArray*    FindLeadingObjects(TObjArray *array );
		void          QSortTracks(TObjArray &a, Int_t first, Int_t last);
                TObjArray*    SortRegions(const AliVParticle* leading, TObjArray *array);
		TObjArray*    GetMinMaxRegion(TList *transv1, TList *transv2);
	
	private:
			
		virtual void AnalyzeESD(AliESDEvent* esd); 
		virtual void AnalyzeMC(AliMCEvent* mc);
		virtual void AnalyzeESDforDCA(AliESDEvent* esdEvent);
		virtual Bool_t isMCEventTrueINEL0(AliMCEvent* fMCEvent);
		virtual ULong64_t GetEventIdAsLong(AliVHeader* header);

		AliMCEvent *fMCEvent;               //! MC event

		AliESDEvent* fESD;                  //!  ESD object
		AliAODEvent* fAOD;                  //!  AOD object
		Bool_t fAnalysisMC;                 //
		AliAnalysisFilter *fTrackFilterDCA; //!
		TString       fAnalysisType;        // "ESD" or "AOD"
		UInt_t       ftrigBit;		    //
		TRandom*      fRandom;              //! random number generator
		Bool_t        fPileUpRej;           // kTRUE is pile-up is rejected
		AliPPVsMultUtils *fPPVsMultUtils;   //!	

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
		Float_t      fZvtx_SPD;           // z vertex
		Int_t        fRun;                // run no
		ULong64_t    fEventId;            // unique event id
		Float_t fdcaxy;			  //
		Float_t fdcaz;		          //
		Bool_t fisPS; 			  //
		Bool_t fisTracklet;		  //
		Bool_t fisMCvtxInZcut;		  //
		Double_t fAveMultiInTrans;      //
		//
		// Output objects
		//
		TList*        fListOfObjects;     //! Output list of objects
		TH1I*         fEvents;            //! No of accepted events
		TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
		TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts
		TH1F* 		fn1;		  //!

		Double_t fNchargedTrue;      	
		Bool_t isINEL0Rec;		
		Bool_t isINEL0True;
		//

		TH1F *INEL0;      		//!

		TH1F *TrigINEL0;		//!
		TH1F *TrueINEL0;		//!

		TH2F *ptvsdcaData;		//!
		TH2F *ptvsdcaPrim;		//!
		TH2F *ptvsdcaDecs;		//!
		TH2F *ptvsdcaMatl;		//!
		

		TH2F *ptvsdcacentralData;	//!
		TH2F *ptvsdcacentralPrim;	//!
		TH2F *ptvsdcacentralDecs;	//!
		TH2F *ptvsdcacentralMatl;	//!

		TH1F *effcomputationGen;	//!
		TH1F *sigLossTrueINEL0;		//!
		TH1F *sigLossTrigINEL0;		//!
		TH1F *nchtrue;			//!

		//for particle composition
		TH1F *effcomputationGenPi;	//!
		TH1F *effcomputationGenK;	//!
		TH1F *effcomputationGenP;	//!
		TH1F *effcomputationGenSm;	//!
		TH1F *effcomputationGenSp;	//!
		TH1F *effcomputationGenO;	//!
		TH1F *effcomputationGenXi;	//!
		TH1F *effcomputationGenL;	//!
		TH1F *effcomputationGenRest;	//!
	
		TH1F *effcomputationRecPi;	//!
		TH1F *effcomputationRecK;	//!
		TH1F *effcomputationRecP;	//!
		TH1F *effcomputationRecSm;	//!
		TH1F *effcomputationRecSp;	//!
		TH1F *effcomputationRecO;	//!
		TH1F *effcomputationRecXi;	//!
		TH1F *effcomputationRecL;	//!
		TH1F *effcomputationRecRest;	//!
	
		TH1F *fPS_MC;			//!
		TH1F *fVtxPS_MC;		//!
		TH1F *fPS;			//!
		TH1F *fVtxPS;			//!

                TH1F *Zvtx;			//!
			
		TH1F *fhRT;                  //!		

		TH1F *secondaries[18];		//!
		TH1F *primariesTrackFilter[18];	//!
		TH1F *effcomputationRec[18];	//!
		TH1F *pti[18];			//!
		TH1F *pti1[18];                  //!
		TH1F *pti2[18];                  //!
		TH1F *pti3[18];                  //!
		TH1F *pti4[18];                  //!
		TH1F *pti5[18];                  //!
		TH1F *ptiMB[18];		 //!
		AliAnalysisFilter* fTrackFilter[18];//! track filter





		AliAnalysisTaskPPvsUE(const AliAnalysisTaskPPvsUE&);            // not implemented
		AliAnalysisTaskPPvsUE& operator=(const AliAnalysisTaskPPvsUE&); // not implemented

		ClassDef(AliAnalysisTaskPPvsUE, 1);    //Analysis task for high pt analysis 
};

#endif
