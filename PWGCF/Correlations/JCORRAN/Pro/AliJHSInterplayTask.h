#ifndef AliJHSInterplay_cxx
#define AliJHSInterplay_cxx

#include "AliAnalysisTaskSE.h"
#include "AliJFFlucAnalysis.h"
#include "AliJHistos.h"
#include "AliJCard.h"
#include "AliJJetTask.h"
#include "AliJFJTask.h"
#include "AliJCatalystTask.h"

const int kNESE = 6;
//#include <fstream>

class AliJHSInterplayTask : public AliAnalysisTaskSE {

	public:
		AliJHSInterplayTask();
		AliJHSInterplayTask(const char *name);
		AliJHSInterplayTask(const AliJHSInterplayTask& a); // not implemented
		AliJHSInterplayTask& operator=(const AliJHSInterplayTask& ap); // not implemented

		virtual ~AliJHSInterplayTask(); 

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);
		virtual Bool_t UserNotify();

		void SetDebugMode( int debug) { fDebugMode = debug; };
		AliJCard *GetCard() { return fCard; }
		void SetCard( AliJCard *c ) { fCard = c; }
		void RegisterList(TClonesArray* listToFill, TObjArray* listFromToFill,double lpt, double hpt);
		void SetPtHardMin( double pthardmin ){fPtHardMin = pthardmin; };
		void SetPtHardMax( double pthardmax ){fPtHardMax = pthardmax; };
		void SetDiJetAsymMin( double min ){fDiJetAsymMin = min; };
		void SetJetTaskName(TString name){ fJetTaskName=name; }
		void SetJetSel(int iS){ fJetSel=iS; }
		 // Methods specific for this class
  		void SetJCatalystTaskName(TString name){ fJCatalystTaskName=name; } // Setter for filter task name
 		void SetJFFlucAnalysis(AliJFFlucAnalysis *ana){ fFFlucAna=ana; } // Setter for analysis
 		void SetJFJTaskName(TString name){ fJFJTaskName=name; } 
 		void Setjettask_tagging(int ijettask) { jettask_tagging=ijettask; }
 		void ESETagging(int itask, int iESE, double lpPT);


	private:
		TDirectory           *fOutput;     // Output
		AliJCatalystTask *fJCatalystTask;  // 
 		TString           fJCatalystTaskName; // Name for JCatalyst task
		AliJFFlucAnalysis *fFFlucAna;
		AliJJetTask           * fJetTask;
		TString fJetTaskName;
		int fJetSel;
		AliJFJTask			  *fJFJTask;
		TString fJFJTaskName;
		AliJCard * fCard;
		AliJHistos *fHistos;
		TClonesArray * fInputListSpectra;
		TClonesArray * fInputListFlow;
		int fVnMethod; // 0; RunFlow 1 : JFluc
		int jettask_tagging;
		int fESMethod; // 0; all pt 1 : leading pt
		int cBin;
		int zBin;
		double zVert;
		Int_t fevt; // event number
		int fDebugMode;
		double fPtHardMin;
		double fPtHardMax;
		double fDiJetAsymMin;
		Bool_t TagThisEvent[kNESE];

		ClassDef(AliJHSInterplayTask, 1); // example of analysis
};

#endif
