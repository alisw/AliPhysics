#ifndef AliJHSInterplay_cxx
#define AliJHSInterplay_cxx

#include <TVector.h>
#include <TRandom.h>
#include <TString.h>
#include <TPRegexp.h>
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisUtils.h"
#include "AliJEfficiency.h"
#include "AliJFFlucAnalysis.h"
#include "AliMCParticle.h"

//#include <fstream>

class AliESDEvent;
class AliMCEvent;
class AliESDtrackCuts;
class AliESDVertex;
class TClonesArray;
class AliJCard;
class AliJHistos;
class AliJBaseTrack;
class AliJEfficiency;
class AliJRunTable;
class AliJFFlucAnalysis;

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

		bool IsGoodEvent(AliVEvent *event);
		void SetDebugMode( int debug) { fDebugMode = debug; };
		AliJCard *GetCard() { return fCard; }
		void SetCard( AliJCard *c ) { fCard = c; }
		void SetIsMC( Bool_t ismc){ IsMC = ismc; cout << "Settint IsMC = " << ismc << endl; };
		double GetCentralityFromImpactPar(double ip);
		double GetImpactParFromCentrality(double cent);
		void ReadKineTracks( AliMCEvent *mcEvent, TClonesArray *TrackList);
		Bool_t IsThisAWeakDecayingParticle(AliMCParticle *thisGuy);
		void SetKineOnly (Bool_t iskineonly) {IsKinematicOnly = iskineonly;};
		void RegisterList(TClonesArray* listToFill, TClonesArray* listFromToFill,double lpt, double hpt);
		void SetPtHardMin( double pthardmin ){fPtHardMin = pthardmin; };
		void SetPtHardMax( double pthardmax ){fPtHardMax = pthardmax; };

	private:
		TDirectory           *fOutput;     // Output

		AliAnalysisUtils *fAnaUtils;
		AliJFFlucAnalysis *fFFlucAna;

		AliJCard * fCard;
		AliJHistos *fHistos;

		AliJEfficiency *fEfficiency;
		int fHadronSelectionCut;

		TClonesArray * fInputList;
		TClonesArray * fInputListSpectra;
		int fVnMethod; // 0; RunFlow 1 : JFluc
		int fESMethod; // 0; all pt 1 : leading pt
		TClonesArray * fInputListFlow;

		Bool_t fFirstEvent; //
		int cBin;
		int zBin;
		double zVert;
		Int_t fevt; // event number
		int fDebugMode;
		Int_t trkfilterBit;
		AliJRunTable *fRunTable; // 

		TRandom *fRandom;
		Bool_t IsMC;
		Bool_t IsKinematicOnly;
		double fPtHardMin;
		double fPtHardMax;
		Bool_t TagThisEvent;

		ClassDef(AliJHSInterplayTask, 1); // example of analysis
};

#endif
