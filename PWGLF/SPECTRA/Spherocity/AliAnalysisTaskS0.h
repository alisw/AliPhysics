#ifndef ALIANALYSISTASKS0_H
#define ALIANALYSISTASKS0_H

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

class AliAnalysisTaskS0 : public AliAnalysisTaskSE {

	public:


		AliAnalysisTaskS0();
		AliAnalysisTaskS0(const char *name);
		virtual ~AliAnalysisTaskS0();

		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);

		virtual void     MakeESDAnalysis( AliESDEvent *event, Double_t etaCut );
		virtual void     MakeESDAnalysisMC( AliMCEvent* fMCEvent, Double_t etaCut );
		virtual void     SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
		virtual void     SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}

		virtual Int_t    GetMultiplicityParticles(Double_t etaCut);
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
		Int_t         fnRefGlobal;
		Double_t      fS0R;
		Double_t      fS0T;
		
		//
		// Output objects
		//

		TList  * fListOfObjects;
		TH1D   * fHistEventCounter;
		TH1D   * hS0;
		TH1D   * hRefMult;
		TH1D   * hpt;
		TH1D   * hphiso;
		TH1D   * hetaso;
		TH2D   * hS0RefMult;
		TH2D   * hS0pt;
		TH2D   * hS0phi;
		TH2D   * hS0eta;

		//True...
		TH1D   * hMultTrue;
		TH1D   * hS0True;
		TH2D   * hS0Truept;
		TH2D   * hS0Truephi;
		TH2D   * hS0Trueeta;
		TH2D   * hS0TrueMult;
		
       		ClassDef(AliAnalysisTaskS0, 11);    //Analysis task for high pt analysis 

};

#endif

