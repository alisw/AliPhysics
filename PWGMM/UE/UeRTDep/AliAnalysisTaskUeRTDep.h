#ifndef ALIANALYSISTASKUERTDEP_H
#define ALIANALYSISTASKUERTDEP_H

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
#include <AliCentrality.h>

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

class AliAnalysisTaskUeRTDep : public AliAnalysisTaskSE {
  
 public:
  
  
		AliAnalysisTaskUeRTDep();
		AliAnalysisTaskUeRTDep(const char *name);
		virtual ~AliAnalysisTaskUeRTDep();

		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);

		virtual void     DataAnaMult( Double_t etaCut );
		virtual void     SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
	
		virtual void     SetTrackCuts(AliAnalysisFilter* fTrackFilter);
		virtual Double_t DeltaPhi(Double_t phia, Double_t phib,
				 Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );

		//02/08/2019  ************
		Double_t GetVtxCut() { return fVtxCut; }   
		Double_t GetEtaCut() { return fEtaCut; }

//	virtual void  SetTrigger(UInt_t ktriggerInt) {ftrigBit = ktriggerInt;}
                virtual void  SetTrigger(UInt_t ktriggerInt = AliVEvent::kINT7) {ftrigBit = ktriggerInt;}  
		virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void  SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;}
		// ********************
	
        protected:

		void		 Exit(const char *msg);
		
	private:
		virtual Bool_t selectVertex2015pp(AliESDEvent *esd,Bool_t checkSPDres, Bool_t requireSPDandTrk,Bool_t checkProximity);
		virtual Bool_t IsGoodSPDvertexRes(const AliESDVertex * spdVertex);

		AliESDEvent*        fESD;
		AliEventCuts        fEventCuts;
		AliAnalysisFilter*  fTrackFilter;

		//02/08/2019  ************
		Double_t     fVtxCut;             // Vtx cut on z position in cm
		Double_t     fEtaCut;             // Eta cut used to select particles
		UInt_t       ftrigBit;
		Short_t      fTriggeredEventMB;   // 1 = triggered, 0 = not trigged (MC only)
		Bool_t fisPS;
		Bool_t fisTracklet;
		Bool_t isINEL0Rec;
		Bool_t fPileUpRej;           // kTRUE is pile-up is rejected
                AliPPVsMultUtils *fUeRTUtils;
                AliMultSelection *fMultSelection;
		
		// ********************
		
		TString       fAnalysisType;
		Int_t fNcl;
              	
		//
		// Output objects
		//

		TList  * fListOfObjects;
		TH1D   * fHistEventCounter;
		TH1D*         fVtxBeforeCuts;     //! Vertex z dist before cuts
		TH1D*         fVtxAfterCuts;      //! Vertex z dist after cuts
		TH1D * hPhi;
		TH1D * hpT;
		TH1D * hDphi;
		TH1D * hDphiNS;
		TH1D * hDphiAS;
		TH1D * hDphiTS;
		TH1D * hPtL;
		TH1D * hEtaL;
		TH1D * hPhiL;
		TH1D * hRefMult08;
		TH1D * hV0Mmult;
		
		TProfile * ProfpTLvsNch;
		TProfile * ProfpTLvsNchNS;
		TProfile * ProfpTLvsNchAS;
		TProfile * ProfpTLvsNchTS;

		TH2D * hpTvsNch;
		TH2D * hpTLvsNch;
		TH2D * hpTLvsNchNS;
		TH2D * hpTLvsNchAS;
		TH2D * hpTLvsNchTS;
		TH2D * hpTvsRefMult08;
		TH2D * hpTLvsRefMult08;
		TH2D * hpTvsV0Mmult;
		TH2D * hpTLvsV0Mmult;
		TH2D * hRefMultvsV0Mmult;
		
		TH3D * hpTvsV0MmultvsRefMult08;
		TH3D * hpTLvsV0MmultvsRefMult08;
		
		TH3D * hpTLvsRefMult08vsDphi;
		TH3D * hpTLvsRefMult08vsNch;
		TH3D * hpTLvsRefMult08vsNchNS;
		TH3D * hpTLvsRefMult08vsNchAS;
		TH3D * hpTLvsRefMult08vsNchTS;
		
		TH3D * hpTLvsV0MmultvsDphi;
		TH3D * hpTLvsV0MmultvsNch;
		TH3D * hpTLvsV0MmultvsNchNS;
		TH3D * hpTLvsV0MmultvsNchAS;
		TH3D * hpTLvsV0MmultvsNchTS;
		
		TH3D * hpTvspTLvsRefMult08;
		TH3D * hpTvspTLvsRefMult08NS;
		TH3D * hpTvspTLvsRefMult08AS;
		TH3D * hpTvspTLvsRefMult08TS;
		
		TH3D * hpTvspTLvsV0Mmult;
		TH3D * hpTvspTLvsV0MmultNS;
		TH3D * hpTvspTLvsV0MmultAS;
		TH3D * hpTvspTLvsV0MmultTS;

		//SumpT vs Nch plot...
		TH2D * hSumptVsRefMult08;
		      		
		Double_t ftrackmult08;       //       
		Double_t fv0mpercentile;     //

		//Plots for RT....
		TH1D * hTrackRT;
		TH1D * hINEL0;
		TH1D * hpTLRT[15];
		TH1D * hpTRT[15];
		TH1D * hpTNSRT[15];
		TH1D * hpTASRT[15];
		TH1D * hpTTSRT[15];
		TH2D * hDEtaDPhiRT[15];
		TH2D * hDeltaphiVspTLRT[15];

		TH2D * hTrackRTvsRefMult08;
		TH2D * hTrackRTvsmultTS;
		TH3D * hpTLvsRefMult08vsRT;

		//Plots for RT for pTL > 5....
		TH1D * hTrackRTpTL5;
		TH1D * hINEL0pTL5;
		TH1D * hpTLRTpTL5[15];
		TH1D * hpTRTpTL5[15];
		TH1D * hpTNSRTpTL5[15];
		TH1D * hpTASRTpTL5[15];
		TH1D * hpTTSRTpTL5[15];
		TH2D * hDEtaDPhiRTpTL5[15];
		TH2D * hDeltaphiVspTLRTpTL5[15];

		TH2D * hTrackRTvsRefMult08pTL5;
		TH2D * hTrackRTvsmultTSpTL5;
		TH3D * hpTLvsRefMult08vsRTpTL5;
		
		ClassDef(AliAnalysisTaskUeRTDep, 11);    //Analysis task for high pt analysis 

};

#endif
