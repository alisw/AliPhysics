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
		virtual void     SetPeriod(const TString period) {fdata_set = period;}

		virtual void     SetTrackCuts(AliAnalysisFilter* fTrackFilter);
		virtual Bool_t   PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh);
		virtual Double_t EtaCalibrationNeg(const Int_t centrality, const Double_t Eta);
		virtual Double_t EtaCalibrationPos(const Int_t centrality, const Double_t Eta);
		virtual Int_t    GetIndexEta(Double_t etain);
		virtual Double_t DeltaPhi(Double_t phia, Double_t phib,
				 Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );
		virtual Int_t    GetIndexPeriod(TString lPeriod);



        protected:

		void		 Exit(const char *msg);
		
	private:

		AliESDEvent*        fESD;
		TString             fdata_set;
		AliEventCuts        fEventCuts;
		AliAnalysisFilter*  fTrackFilter;

		TString       fAnalysisType;
		TF1* fcutLow;
		TF1* fcutHigh;
		const Double_t fDeDxMIPMin;
		const Double_t fDeDxMIPMax;
		TF1 *fEtaCalibrationNeg;
		TF1 *fEtaCalibrationPos;
		Int_t fNcl;
		AliMultSelection *fMultSelection;
		//
		// Output objects
		//

		TList  * fListOfObjects;
		TH1D   * fHistEventCounter;
		TH1I   * fEvents;
		TH2D * hPhi;
		TH2D * hMIPVsEta;
		TProfile * pMIPVsEta;
		TH1D * hpT;
		TH1D * hDphi;
		TH1D * hDphiNNS; // R<0.5
		TH1D * hDphiNS;
		TH1D * hDphiAS;
		TH1D * hDphiTS;

		TH1D * hPtL[5];
		TH2D * hDeDxL[5];
		TH1D * hEtaL[5];
		TH1D * hPhiL[5];

		TProfile * pNNS[5]; // R<0.5
		TProfile * pNS[5];
		TProfile * pAS[5];
		TProfile * pTS[5];

		TH2D * hNNS[5]; // R<0.5
		TH2D * hNS[5];
		TH2D * hAS[5];
		TH2D * hTS[5];

		TH3D * h3DNNS[5]; // R<0.5
		TH3D * h3DNS[5];
		TH3D * h3DAS[5];
		TH3D * h3DTS[5];
		TH2D * hptVSp[5];

		TH1D * hRefMult08;
		TH1D * hV0Mmult;
		TH2D * hpTvsNch;
		TH2D * hpTLvsNch;
		TH2D * hpTvsRefMult08;
		TH2D * hpTLvsRefMult08;
		TH2D * hpTvsV0Mmult;
		TH2D * hpTLvsV0Mmult;
		TH2D * hRefMultvsV0Mmult;
		TH3D * hpTvsV0MmultvsRefMult08;
		TH3D * hpTLvsV0MmultvsRefMult08;
		TH3D * hpTLvsRefMult08vsDphi;
		TH3D * hpTLvsV0MmultvsDphi;
		
		TH1D * hDphiPtLBins[10];
		TH2D * hDetaDphiPtLBins[10];

		THnSparseD * hptLvsv0MvsRefMultvsDphivsDeta;       //

		Double_t ftrackmult08;       //       
		Double_t fv0mpercentile;     //

		ClassDef(AliAnalysisTaskUeMultDep, 11);    //Analysis task for high pt analysis 

};

#endif

