#ifndef ALIANALYSISTASKESA_H
#define ALIANALYSISTASKESA_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TRandom.h>
#include <TObject.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliMCEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliESDtrackCuts.h>
#include <AliSpherocityUtils.h>
#include <AliEventCuts.h>

class AliAnalysisUeSpherocityTask : public AliAnalysisTaskSE {

	public:


		AliAnalysisUeSpherocityTask();
		AliAnalysisUeSpherocityTask(const char *name);

		virtual ~AliAnalysisUeSpherocityTask();

		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);

		virtual void     MakeESDAnalysis( AliESDEvent *event, Double_t etaCut );
		virtual void     SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void     SetSpherocityPercBinning(double* soB) {fsoB = soB;}
		virtual void     SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}

		virtual Int_t    GetMultiplicityParticles(Double_t etaCut);
		virtual Double_t GetSpheroPercentile( Double_t valES, Int_t valMult );
		virtual Int_t    GetSpheroPercentileBin( Double_t valES );
		virtual void     SetTrackCuts(AliAnalysisFilter* fTrackFilter);

	protected:

		void		 Exit(const char *msg);

	private:

		AliESDEvent*        fESD;                  //! ESD object
		AliEventCuts        fEventCuts;
		AliMCEvent*         fMC;                   //! MC object
		AliStack*           fMCStack;              //! MC ESD stack
		AliSpherocityUtils* fSpheroUtils; 
		AliAnalysisFilter*  fTrackFilter;

		TString       fAnalysisType;        //  "ESD" or "AOD"
		Bool_t        fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
		Int_t         fnMultbins;           // number of multiplicity bins
		Int_t         fnRefGlobal;
		double*       fsoB;
		Int_t         fnsoB;            // number of bins spherocity percentile
		Int_t         fbinSom;          // spherocity bin

		//
		// Helper histograms
		//
		TH1D * hSOMPerc[100];
		TH2D * hSOMAux[100];
		//
		// Output objects
		//

		TList  * fListOfObjects;     //! Output list of objects
		TH1D   * fHistEventCounter; //!
		TH1I   * fEvents;            //! No of accepted events
		TH2D   * hSOGlobal08;
		TH2D   * hSOGlobal[200];    // spherocity distributions for different percentiles
		TH2D   * hptRefGlobal[200];
		TH2D   * hptRefGlobal08;    // pt vs multiplicity, MB 
		TH2D   * hSOTrue08;
		TH2D   * hRMGlobal;
		TH1D   * hphiso;
		TH1D   * hetaso;
		TH1F   * fn1;

		ClassDef(AliAnalysisUeSpherocityTask, 11);    //Analysis task for high pt analysis 

};

#endif

