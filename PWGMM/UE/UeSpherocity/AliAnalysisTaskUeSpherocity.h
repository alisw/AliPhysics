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

class AliAnalysisTaskUeSpherocity : public AliAnalysisTaskSE {

	public:


		AliAnalysisTaskUeSpherocity();
		AliAnalysisTaskUeSpherocity(const char *name);
		virtual ~AliAnalysisTaskUeSpherocity();

		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);

		virtual void     MakeESDAnalysis( AliESDEvent *event, Double_t etaCut );
		virtual void     SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void     SetHisto(TH1D *hBining) {fSoBining = hBining;}
		virtual void     SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}

		virtual Int_t    GetMultiplicityParticles(Double_t etaCut);
		virtual Double_t GetSpheroPercentile( Double_t valES, Int_t valMult );
		virtual Int_t    GetSpheroPercentileBin( Double_t valES );
		virtual void     SetTrackCuts(AliAnalysisFilter* fTrackFilter);

	protected:

		void		 Exit(const char *msg);

	private:

		AliESDEvent*        fESD;
		AliEventCuts        fEventCuts;
		AliMCEvent*         fMC;
		AliStack*           fMCStack;
		AliSpherocityUtils* fSpheroUtils;
		AliAnalysisFilter*  fTrackFilter;

		TString       fAnalysisType;
		Bool_t        fAnalysisMC;
		Int_t         fnMultbins;
		Int_t         fnRefGlobal;
		TH1D        * fSoBining;
		Int_t         fnsoB;
		Int_t         fbinSom;

		//
		// Output objects
		//

		TList  * fListOfObjects;
		TH1D   * fHistEventCounter;
		TH1I   * fEvents;
		TH2D   * hSOGlobal08;
		TH2D   * hSOGlobal[200];
		TH2D   * hptRefGlobal[200];
		TH2D   * hptRefGlobal08; 
		TH2D   * hSOTrue08;
		TH2D   * hRMGlobal;
		TH1D   * hphiso;
		TH1D   * hetaso;
		TH1F   * fn1;

		ClassDef(AliAnalysisTaskUeSpherocity, 11);    //Analysis task for high pt analysis 

};

#endif

