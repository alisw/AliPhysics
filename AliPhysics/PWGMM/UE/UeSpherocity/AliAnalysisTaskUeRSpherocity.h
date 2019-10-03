#ifndef ALIANALYSISTASKUERSPHEROCITY_H
#define ALIANALYSISTASKUERSPHEROCITY_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TObject.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliAnalysisFilter.h>
#include <AliVHeader.h>
#include <AliESDtrackCuts.h>
#include <AliSpherocityUtils.h>
#include <AliEventCuts.h>

class AliAnalysisTaskUeRSpherocity : public AliAnalysisTaskSE {

	public:


		AliAnalysisTaskUeRSpherocity();
		AliAnalysisTaskUeRSpherocity(const char *name);
		virtual ~AliAnalysisTaskUeRSpherocity();

		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);
		virtual int      IndexLeadingTrack( Double_t etaCut );
		virtual bool     MakeAnalysis( Int_t index_sample, Int_t index_leading, Double_t etaCut );
		virtual void     SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void     SetPeriod(const TString period) {fdata_set = period;}

		virtual void     SetTrackCuts(AliAnalysisFilter* fTrackFilter);
		virtual Double_t DeltaPhi(Double_t phia, Double_t phib,
				Double_t rangeMin = -TMath::Pi()/2, Double_t rangeMax = 3*TMath::Pi()/2 );


	protected:

		void		 Exit(const char *msg);

	private:

		AliESDEvent*        fESD;
		TString             fdata_set;
		AliEventCuts        fEventCuts;
		AliAnalysisFilter*  fTrackFilter;
		AliSpherocityUtils* fSpheroUtils;


		TString       fAnalysisType;
		//
		// Output objects
		//

		TList  * fListOfObjects;
		TH1D   * fHistEventCounter;
		TH1D * hcounter;
		TH1D   * hphiso;
		TH1D   * hetaso;
		TH1D * hR[5];
		TH1D * hPtLeading[5];
		TH2D * hRDPhi[7][5];
		TH2D * hPtVsR[7][5];

		AliAnalysisTaskUeRSpherocity(const AliAnalysisTaskUeRSpherocity&);            // not implemented
		AliAnalysisTaskUeRSpherocity& operator=(const AliAnalysisTaskUeRSpherocity&); // not implemented
		ClassDef(AliAnalysisTaskUeRSpherocity, 11);    //Analysis task for high pt analysis 

};

#endif

