#ifndef ALIANALYSISTASKUEMULTDEP_H
#define ALIANALYSISTASKUEMULTDEP_H

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

class AliAnalysisTaskUeMultDep : public AliAnalysisTaskSE {
  
 public:
  
  
		AliAnalysisTaskUeMultDep();
		AliAnalysisTaskUeMultDep(const char *name);
		virtual ~AliAnalysisTaskUeMultDep();

		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);

		virtual void     MakeAnalysis( Double_t etaCut );
		virtual void     SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
	
		virtual void     SetTrackCuts(AliAnalysisFilter* fTrackFilter);
		virtual Double_t DeltaPhi(Double_t phia, Double_t phib,
				 Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
	
        protected:

		void		 Exit(const char *msg);
		
	private:

		AliESDEvent*        fESD;
		AliEventCuts        fEventCuts;
		AliAnalysisFilter*  fTrackFilter;

		TString       fAnalysisType;
		Int_t fNcl;
		AliMultSelection *fMultSelection;
		//
		// Output objects
		//

		TList  * fListOfObjects;
		TH1D   * fHistEventCounter;
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
		
		Double_t ftrackmult08;       //       
		Double_t fv0mpercentile;     //

		ClassDef(AliAnalysisTaskUeMultDep, 11);    //Analysis task for high pt analysis 

};

#endif
