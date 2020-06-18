#ifndef ALIANALYSISTASKUESPECTRADPHI_H
#define ALIANALYSISTASKUESPECTRADPHI_H

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
#include <AliMCEvent.h>
#include <AliMCEventHandler.h>
#include <AliESDVertex.h>
#include <AliAnalysisFilter.h>
#include <AliVHeader.h>
#include <AliESDtrackCuts.h>
#include <AliEventCuts.h>
#include <AliMultSelection.h>
#include <AliCentrality.h>
#include <AliMultSelectionTask.h>

class AliAnalysisTaskUeSpectraDphi : public AliAnalysisTaskSE {

 public:
		AliAnalysisTaskUeSpectraDphi();
		AliAnalysisTaskUeSpectraDphi(const char *name);
		virtual ~AliAnalysisTaskUeSpectraDphi();

		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);

		virtual void     AnalyzeESD(AliESDEvent* esd);
		virtual void     AnalyzeMC();
		//virtual void     SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void     SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
		virtual void     SetTrackCuts(AliAnalysisFilter* fTrackFilter);
		virtual Double_t DeltaPhi(Double_t phia, Double_t phib,
				 Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );

		Double_t GetVtxCut() { return fVtxCut; }
		Double_t GetEtaCut() { return fEtaCut; }
		Double_t GetTPCNclCut() { return fNcl; }
		virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void  SetTPCNclCut(Double_t nclCut){fNcl = nclCut;}

		virtual void  SetTrigger(UInt_t ktriggerInt = AliVEvent::kINT7) {ftrigBit = ktriggerInt;}
		virtual void  SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;}

 protected:
		void		 Exit(const char *msg);

 private:
		virtual Bool_t selectVertex2015pp(AliESDEvent *esd,Bool_t checkSPDres, Bool_t requireSPDandTrk,Bool_t checkProximity);
		virtual Bool_t IsGoodSPDvertexRes(const AliESDVertex * spdVertex);
		virtual Bool_t isMCEventTrueINEL0(AliMCEvent* fMCEvent);

		//
		// Cuts and options
		//
		Int_t        fNcl;
		Double_t     fVtxCut;             // Vtx cut on z position in cm
		Double_t     fEtaCut;             // Eta cut used to select particles

		AliMCEvent *        fMCEvent;               //! MC event
		AliStack   *        fMCStack;
		AliESDEvent*        fESD;
		AliEventCuts        fEventCuts;
		AliAnalysisFilter*  fTrackFilter;
		//TString       fAnalysisType;
		AliMultSelection *fMultSelection;

		UInt_t        ftrigBit;		    //
		Bool_t        fPileUpRej;           // kTRUE is pile-up is rejected
		Bool_t        fAnalysisMC;

		//
		// Help variables
		//
		Short_t      fTriggeredEventMB;   // 1 = triggered, 0 = not trigged (MC only) */
	        Double_t     fZvtx_SPD;           // z vertex */
		Bool_t       fisPS; 			  //
		Bool_t       fisTracklet;		  //
		Bool_t       isINEL0Rec;
		Bool_t       isINEL0True;
		Bool_t       fMCVzCut;

		//
		// Output objects
		//
		TList* fListOfObjects;
		TH1D * fVtxBeforeCuts;
		TH1D * fVtxAfterCuts;
		TH1D * hSelEv;
		TH1D * hINEL0;

		TH1D * hPS;			//!
		TH1D * hVtxPS;			//!
		
		TH1D * hpT;
		TH1D * hEta;
		TH1D * hPhi;
		TH1D * hPtL;
		TH1D * hEtaL;
		TH1D * hPhiL;
		TH1D * hDphi;
		TH1D * hRefMult08;
		TH1D * hV0Mmult;

		TH2D * hpTvsDphiOA;
		TH2D * hpTvsDphiSA;
		TH2D * hMultvsDphiOA;
		TH2D * hMultvsDphiSA;
		TH3D * hMultvspTvsDphi;

		TH1D * hpTDphiBinsSAWLP; //With leading particle

		TH1D * hpTDphiBinsOA[18];
		TH1D * hpTDphiBinsSA[18];
		TH1D * hMultDphiBinsOA[18];
		TH1D * hMultDphiBinsSA[18];

		//MC....
		TH1D * hINEL0MCTrig;
		TH1D * hINEL0MCTrue;
		TH1D * hPS_MC;			//!
		TH1D * hVtxPS_MC;		//!
		TH1D * hpTMCTrue;
		TH1D * hEtaMCTrue;
		TH1D * hPhiMCTrue;
		TH1D * hPtLMCTrue;
		TH1D * hEtaLMCTrue;
		TH1D * hPhiLMCTrue;
		TH1D * hDphiMCTrue;
		TH2D * hpTvsDphiOAMCTrue;
		TH2D * hpTvsDphiSAMCTrue;
		TH2D * hMultvsDphiOAMCTrue;
		TH2D * hMultvsDphiSAMCTrue;
		TH3D * hMultvspTvsDphiMCTrue;

		TH1D *hpTDphiBinsSAWLPMCTrue; //With leading particle
		
		TH1D *hpTDphiBinsOAMCTrue[18];
		TH1D *hpTDphiBinsSAMCTrue[18];
		TH1D *hMultDphiBinsOAMCTrue[18];
		TH1D *hMultDphiBinsSAMCTrue[18];
		
		Double_t ftrackmult08;       //
		Double_t fv0mpercentile;     //

		ClassDef(AliAnalysisTaskUeSpectraDphi, 1);    //Analysis task for high pt analysis

};

#endif
