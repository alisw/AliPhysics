#ifndef ALIANALYSISTASKSPECTRADPHI_H
#define ALIANALYSISTASKSPECTRADPHI_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH3D.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <TRandom.h>
#include <TObject.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliAnalysisFilter.h>
#include <AliVHeader.h>
#include <AliESDtrackCuts.h>
#include <AliEventCuts.h>
#include <AliMultSelection.h>
#include <AliCentrality.h>
#include <AliMultSelectionTask.h>


// class AliPPVsMultUtils;

class AliAnalysisTaskSpectraDPhi : public AliAnalysisTaskSE {
  
 public:
  
  
		AliAnalysisTaskSpectraDPhi();
		AliAnalysisTaskSpectraDPhi(const char *name);
		virtual ~AliAnalysisTaskSpectraDPhi();

		Double_t GetVtxCut() { return fVtxCut; }   
		Double_t GetEtaCut() { return fEtaCut; }     
                Double_t GetLeadMin() { return fLeadMin; }
		
		virtual void  SetTrigger(UInt_t ktriggerInt = AliVEvent::kINT7) {ftrigBit = ktriggerInt;}  
		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);

		virtual void     AnalysisSpectra(AliESDEvent* esd);
		//	virtual void     AnalyseDataRT(AliESDEvent* esd);
		virtual void     SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void     SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
		virtual void     SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void     SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void     SetPtLeadMin(Double_t leadMin){fLeadMin = leadMin;}
		virtual void     SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;}
		//	virtual void     SetAveMultiInTrans(Double_t value) {fAveMultiInTrans = value;}
	
		virtual void     SetTrackCuts(AliAnalysisFilter* fTrackFilter);
		virtual Double_t DeltaPhi(Double_t phia, Double_t phibs);
	       	virtual Bool_t   selectVertex2015pp(AliESDEvent *esd, Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
		virtual Bool_t   IsGoodSPDvertexRes(const AliESDVertex * spdVertex);
		virtual ULong64_t GetEventIdAsLong(AliVHeader* header);
			
        protected:

		void		 Exit(const char *msg);
		
	private:
		AliESDEvent* fESD;                  //!  ESD object
		AliAODEvent* fAOD;                  //!  AOD object
	  
		AliAnalysisFilter *fTrackFilter; //!
		AliAnalysisFilter *fTrackFilterDCA; //!
		//	AliAnalysisFilter *fTrackFilterMatchEff; //!

		TString      fAnalysisType;        // "ESD" or "AOD"
		UInt_t       ftrigBit;		   //
		Bool_t       fPileUpRej;           // kTRUE is pile-up is rejected
		Bool_t       fAnalysisMC;
		Bool_t       isINEL0Rec;
			
		AliEventCuts     fEventCuts;
		AliMultSelection *fMultSelection;

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
		Float_t      fZvtx;               // z vertex
		Int_t        fRun;                // run no
		ULong64_t    fEventId;            // unique event id
		Float_t fdcaxy;			  //
		Float_t fdcaz;		          //
		Bool_t fisPS; 			  //
		//	Double_t fAveMultiInTrans;      //
		TRandom*      fRandom;              //! random number generator


		//
		// Output objects
		//
		TList *  fListOfObjects;     //! Output list of objects
		TH1I  *  fEvents;            //! No of accepted events
		TH1D  *  fVtxBeforeCuts;     //! Vertex z dist before cuts
		TH1D  *  fVtxAfterCuts;      //! Vertex z dist after cuts
		TH1D  *  fEventCounter;	     //!
		
		TH1D * hpT;
		TH1D * hPhi;
		TH1D * hEta;
			
		TH1D * hPtL;
		TH1D * hEtaL;
		TH1D * hPhiL;
		
		TH1D * hDphi;
		TH1D * hRefMult08;
		TH1D * hV0Mmult;

		TH2D *hpTvsDphiOA;
		TH2D *hpTvsDphiSA;

		TH1D *hpTDphiBinsOA[100];
		TH1D *hpTDphiBinsSA[100];
		TH1D *hMultDphiBinsOA[100];
		TH1D *hMultDphiBinsSA[100];
				
		Double_t ftrackmult08;       //       
		Double_t fv0mpercentile;     //

		ClassDef(AliAnalysisTaskSpectraDPhi, 11);    //Analysis task for high pt analysis 

};

#endif
