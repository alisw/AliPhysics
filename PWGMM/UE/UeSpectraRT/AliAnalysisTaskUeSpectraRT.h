#ifndef AliAnalysisTaskUeSpectraRT_H
#define AliAnalysisTaskUeSpectraRT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */
//Authors: Sergio Iga ,sergio.iga@correo.nucleares.unam.mx
//         Valentina Zaccolo, valentina.zaccolo@cern.ch
//         Aditya Nath Mishra, Aditya.Nath.Mishra@cern.ch

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


// #include "AliPPVsMultUtils.h"
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


// class AliPPVsMultUtils;
class AliTransverseEventShape;
//class AliMCSpectraWeights;

class AliAnalysisTaskUeSpectraRT : public AliAnalysisTaskSE {
	public:


		AliAnalysisTaskUeSpectraRT();
		AliAnalysisTaskUeSpectraRT(const char *name);
		virtual ~AliAnalysisTaskUeSpectraRT();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);

		Double_t GetVtxCut() { return fVtxCut; }
		Double_t GetEtaCut() { return fEtaCut; }
                Double_t GetLeadMin() { return fLeadMin; }

		virtual void  SetTrigger(UInt_t ktriggerInt = AliVEvent::kINT7) {ftrigBit = ktriggerInt;}
		virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
		virtual void  SetAnalysisCorr(Bool_t isCorr) {fAnalysisCorr = isCorr;}
		virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
                virtual void  SetPtLeadMin(Double_t leadMin){fLeadMin = leadMin;}
		virtual void  SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;}
		virtual void  SetAveMultiInTrans(Double_t value) {fAveMultiInTrans = value;}
                virtual void  SetAveRecMultiInTrans(Double_t value) {fAveRecMultiInTrans = value;}
		virtual void  SetAveGenMultiInTrans(Double_t value) {fAveGenMultiInTrans = value;}
		//virtual void  SetMCSpectraWeightObject(AliMCSpectraWeights *obj) {fMCSpectraWeights = obj;}

		TObjArray*    FindLeadingObjects(TObjArray *array );
		void          QSortTracks(TObjArray &a, Int_t first, Int_t last);
                TObjArray*    SortRegions(const AliVParticle* leading, TObjArray *array);
		TObjArray*    GetMinMaxRegion(TList *transv1, TList *transv2);
		TObjArray*    GetRegionAwTow(TList *region);
  		void          FillRTResponseMatrix(const AliVParticle* leadingMC, const AliVParticle* leading, TList* listMaxMC, TList* listMax, TList* listMinMC, TList* listMin);
	private:

		virtual void AnalyseDataRT(AliESDEvent* esd);
		virtual void CorrectionsDataRT(AliESDEvent* esd, Bool_t isVtxGood);
		virtual void CorrectionsMCRT(AliMCEvent* mc, AliESDEvent* esd);
		virtual void AnalyzeESD(AliESDEvent* esd);
		virtual void AnalyzeMC(AliMCEvent* mc);
		virtual void AnalyzeESDforDCA(AliESDEvent* esdEvent);
		virtual Bool_t selectVertex2015pp(AliESDEvent *esd,Bool_t checkSPDres, Bool_t requireSPDandTrk,Bool_t checkProximity);
		virtual Bool_t IsGoodSPDvertexRes(const AliESDVertex * spdVertex);
		virtual Bool_t isMCEventTrueINEL0(AliMCEvent* fMCEvent);
		virtual ULong64_t GetEventIdAsLong(AliVHeader* header);

		AliMCEvent *fMCEvent;               //! MC event

		AliESDEvent* fESD;                  //!  ESD object
		AliAODEvent* fAOD;                  //!  AOD object
		Bool_t fAnalysisMC;                 //
		Bool_t fAnalysisCorr;                 //
		AliAnalysisFilter *fTrackFilterDCA; //!
		AliAnalysisFilter *fTrackFilterMatchEff; //!
		TString       fAnalysisType;        // "ESD" or "AOD"
		UInt_t       ftrigBit;		    //
		TRandom*      fRandom;              //! random number generator
		Bool_t        fPileUpRej;           // kTRUE is pile-up is rejected
//		AliPPVsMultUtils *fPPVsMultUtils;   //!
		//AliMCSpectraWeights *fMCSpectraWeights;//!
		AliMultSelection *fMultSelection;  //!

		//
		// Cuts and options
		//

		Double_t     fVtxCut;             // Vtx cut on z position in cm
		Double_t     fEtaCut;             // Eta cut used to select particles
                Double_t     fLeadMin;             // Low limit to pT of the leading particle
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
		Double_t fAveRecMultiInTrans;   //
                Double_t fAveGenMultiInTrans;   //
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
                TH1F *INEL0Gen;                    //!

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
                TH1F *effcomputationGen1;        //!
                TH1F *effcomputationGen2;        //!
                TH1F *effcomputationGen3;        //!
                TH1F *effcomputationGen4;        //!
                TH1F *effcomputationGen5;        //!
		TH1F *sigLossTrueINEL0;		//!
		TH1F *sigLossTrigINEL0;		//!
                TH1F *sigLossTrueINEL01;         //!
                TH1F *sigLossTrigINEL01;         //!
                TH1F *sigLossTrueINEL02;         //!
                TH1F *sigLossTrigINEL02;         //!
                TH1F *sigLossTrueINEL03;         //!
                TH1F *sigLossTrigINEL03;         //!
                TH1F *sigLossTrueINEL04;         //!
                TH1F *sigLossTrigINEL04;         //!
                TH1F *sigLossTrueINEL05;         //!
                TH1F *sigLossTrigINEL05;         //!
		TH1F *sigLossTrueINEL06;         //!
                TH1F *sigLossTrigINEL06;         //!
                TH1F *sigLossTrueINEL07;         //!
                TH1F *sigLossTrigINEL07;         //!
                TH1F *sigLossTrueINEL08;         //!
                TH1F *sigLossTrigINEL08;         //!
                TH1F *sigLossTrueINEL09;         //!
                TH1F *sigLossTrigINEL09;         //!
                TH1F *sigLossTrueINEL010;         //!
                TH1F *sigLossTrigINEL010;         //!
		TH1F *sigLossTrueINEL011;         //!
                TH1F *sigLossTrigINEL011;         //!
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

                TH1F *effcomputationRec;    //!
                TH1F *effcomputationRec1;    //!
                TH1F *effcomputationRec2;    //!
                TH1F *effcomputationRec3;    //!
                TH1F *effcomputationRec4;    //!
                TH1F *effcomputationRec5;    //!
		TH1F *effcomputationRec6;    //!
		TH1F *effcomputationRec7;    //!
		TH1F *effcomputationRec8;    //!
		TH1F *effcomputationRec9;    //!
		TH1F *effcomputationRec10;    //!
		TH1F *effcomputationRec11;    //!

		TH1F *fPS_MC;			//!
		TH1F *fVtxPS_MC;		//!
		TH1F *fPS;			//!
		TH1F *fVtxPS;			//!

                TH1F *Zvtx;			//!
  		TH2F *Phi;                      //!
  		TH2F *Eta;                      //!
  		TH2F *EtaPhi;                   //!

		TH1F *fhRTData;                  //!
                TH1F *fhRTReco;                  //!
                TH1F *fhRTTrue;                  //!

		TH2F *fhRTResponse;              //!

		TH1F *secondaries[18];		//!
		TH1F *primariesTrackFilter[18];	//!
		TH1F *pti[18];
		TH1F *ptli1;                  //!
		TH1F *ptli2;                  //!
		TH1F *ptli3;                  //!
		TH1F *ptli4;                  //!
		TH1F *ptli5;                  //!
		TH1F *ptli6;                  //!
		TH1F *ptli7;                  //!
		TH1F *ptli8;                  //!
		TH1F *ptli9;                  //!
		TH1F *ptli10;                  //!
		TH1F *ptli11;                  //!
		TH1F *ptliMB;		       //!
		TH1F *pti1Trans[18];                  //!
		TH1F *pti2Trans[18];                  //!
		TH1F *pti3Trans[18];                  //!
		TH1F *pti4Trans[18];                  //!
		TH1F *pti5Trans[18];                  //!
		TH1F *pti6Trans[18];                  //!
		TH1F *pti7Trans[18];                  //!
		TH1F *pti8Trans[18];                  //!
		TH1F *pti9Trans[18];                  //!
		TH1F *pti10Trans[18];                  //!
		TH1F *pti11Trans[18];                  //!
		TH1F *ptiMBTrans[18];		 //!
		TH1F *pti1Tow[18];                  //!
                TH1F *pti2Tow[18];                  //!
                TH1F *pti3Tow[18];                  //!
                TH1F *pti4Tow[18];                  //!
                TH1F *pti6Tow[18];                  //!
		TH1F *pti5Tow[18];                  //!
		TH1F *pti7Tow[18];                  //!
		TH1F *pti8Tow[18];                  //!
		TH1F *pti9Tow[18];                  //!
		TH1F *pti10Tow[18];                  //!
		TH1F *pti11Tow[18];                  //!
		TH1F *ptiMBTow[18];                 //!
                TH1F *pti1Aw[18];                  //!
                TH1F *pti2Aw[18];                  //!
                TH1F *pti3Aw[18];                  //!
                TH1F *pti4Aw[18];                  //!
                TH1F *pti5Aw[18];                  //!
		TH1F *pti6Aw[18];                  //!
		TH1F *pti7Aw[18];                  //!
		TH1F *pti8Aw[18];                  //!
		TH1F *pti9Aw[18];                  //!
		TH1F *pti10Aw[18];                  //!
		TH1F *pti11Aw[18];                  //!
		TH1F *ptiMBAw[18];                 //!
		TH1F *ptiNoLead[18];             //!

		TH1F *primariesTrackFilterME;	 //!
		TH1F *ptiME;                      //!

		TH1F *pTGen;			//!
                TH1F *pTGenTrans;                    //!
                TH1F *pTGenTow;                    //!
                TH1F *pTGenAw;                    //!
		TH1F *ptliGen1;                  //!
		TH1F *ptliGen2;                  //!
		TH1F *ptliGen3;                  //!
		TH1F *ptliGen4;                  //!
		TH1F *ptliGen5;                  //!
		TH1F *ptliGen6;                  //!
		TH1F *ptliGen7;                  //!
		TH1F *ptliGen8;                  //!
		TH1F *ptliGen9;                  //!
		TH1F *ptliGen10;                  //!
		TH1F *ptliGen11;                  //!
		TH1F *pTGenTrans1;                    //!
                TH1F *pTGenTrans2;                    //!
                TH1F *pTGenTrans3;                    //!
                TH1F *pTGenTrans4;                    //!
                TH1F *pTGenTrans5;                    //!
		TH1F *pTGenTrans6;                    //!
		TH1F *pTGenTrans7;                    //!
		TH1F *pTGenTrans8;                    //!
		TH1F *pTGenTrans9;                    //!
		TH1F *pTGenTrans10;                    //!
		TH1F *pTGenTrans11;                    //!
                TH1F *pTGenTow1;                    //!
                TH1F *pTGenTow2;                    //!
                TH1F *pTGenTow3;                    //!
                TH1F *pTGenTow4;                    //!
                TH1F *pTGenTow5;                    //!
		TH1F *pTGenTow6;                    //!
		TH1F *pTGenTow7;                    //!
		TH1F *pTGenTow8;                    //!
		TH1F *pTGenTow9;                    //!
		TH1F *pTGenTow10;                    //!
		TH1F *pTGenTow11;                    //!
                TH1F *pTGenAw1;                    //!
                TH1F *pTGenAw2;                    //!
                TH1F *pTGenAw3;                    //!
                TH1F *pTGenAw4;                    //!
                TH1F *pTGenAw5;                    //!
		TH1F *pTGenAw6;                    //!
		TH1F *pTGenAw7;                    //!
		TH1F *pTGenAw8;                    //!
		TH1F *pTGenAw9;                    //!
		TH1F *pTGenAw10;                    //!
		TH1F *pTGenAw11;                    //!

		//Mult
		TH1F * hRefMult08;
		TH1F * hV0Mmult;

		Double_t ftrackmult08;       //
		Double_t fv0mpercentile;     //


		AliAnalysisFilter* fTrackFilter[18];//! track filter





		AliAnalysisTaskUeSpectraRT(const AliAnalysisTaskUeSpectraRT&);            // not implemented
		AliAnalysisTaskUeSpectraRT& operator=(const AliAnalysisTaskUeSpectraRT&); // not implemented

		ClassDef(AliAnalysisTaskUeSpectraRT, 1);    //Analysis task for high pt analysis
};

#endif
