#ifndef ALIANALYSISMULTIPLICITYTASK_H
#define ALIANALYSISMULTIPLICITYTASK_H

#include <AliAnalysisTaskSE.h>
#include <AliESDtrackCuts.h>
#include <AliESDEvent.h>
#include <AliMCEvent.h>
#include <AliGenCocktail.h>
#include <AliGenCocktailEventHeader.h>
#include "AliMultiplicityHelper.h"
#include "AliMultiplicityEventSelection.h"
#include <map>
#include <TObjString.h>

const Double_t ITSetaLimit = 1.3;
const Double_t TPCetaLimit = 0.9;

class AliAnalysisMultiplicityTask : public AliAnalysisTaskSE {
public:
	enum ESelection {
		kExclusiveSelection,
		kAnalysisSelection
	};
	
	void CreateOutput();														// output list initialization

	virtual void UserCreateOutputObjects();												// post the output data
	virtual void UserExec ( Option_t* );												// event processing function
	virtual void Terminate ( Option_t* );												// analysis finishing
	void SetUseMC ( Bool_t use );													// set MC flag
	Bool_t UseMC();															// check the MC flag

	void AddExclusiveSelection ( AliMultiplicityEventSelection* selection, const Int_t outputSlot = 1 );								//add selection to the list
	void AddAnalysisSelection ( AliMultiplicityEventSelection* selection, const Int_t outputSlot = 1 );								//add selection to the list

	void AddTrackCuts ( AliESDtrackCuts* externalTrackCuts = 0 );
	void AddTrackCuts ( TList* cutsList );

	AliAnalysisMultiplicityTask ( const char* name = "AliAnalysisMultiplicityTask", const Int_t nOutput = 1 );

	void DumpLists();

// 	void SetUseMultInV0 ( Bool_t use = kTRUE );
// 	Bool_t GetUseMultInV0 ( );

	void SetEventBits ( UInt_t mask ) {
		fEventClassification = mask;
	};
	UInt_t GetEventBits () {
		return fEventClassification;
	};
	void SetEventBit ( UInt_t bit ) {
		fEventClassification |= bit/* & ~((UInt_t) 0)*/;
	};
	Bool_t TestEventBit ( UInt_t bit ) const {
		return ( Bool_t ) ( ( fEventClassification & bit ) != 0 );
	};
	void InvertEventBit ( UInt_t bit ) {
		fEventClassification ^= bit/* & ~((UInt_t) 0)*/;
	};

	virtual void SetConsistencyCut ( Double_t cut ) {
		fConsistencyCut = cut;
	};
	virtual Double_t GetConsistencyCut() const {
		return fConsistencyCut;
	};
// 	virtual void SetRequireConsistency ( Bool_t req = kTRUE ) {
// 		fRequireConsistency = req;
// 	};
	
	void SetSPDPileupThreshold ( Double_t threshold = 0.8 ) {
		fSPDPileupThreshold = threshold;
	};
	
protected:
	TList* fExclusiveSelectionList;													// list of exclusive selections
	TList* fAnalysisSelectionList;													// list of analysis selections
	TList* fOutput; 														// Output list
	TList* fTrackCuts;														// ESD track cuts list

	TList* fCollectedEvents;													// event list output
private:
	UInt_t GetEventClassification () {
		return fEventClassification;
	};
	
	Bool_t fUseMC;											// MC access flag
	Int_t fNOutput; // number of output slots to be used by the task
	std::map<TString, Int_t> fOutputMap; // map selections to output slots
	Bool_t fIsCocktail;

	AliESDEvent* ESDEvent();													// get the pointer to the current ESD

	void Multiplicity ( Bool_t excludeESDEvent, Bool_t excludeMCEvent );								 // multiplicity loop
	void AddSelection ( ESelection selType, AliMultiplicityEventSelection* selection, const Int_t outputSlot );
	void PostDataLists();

	virtual ~AliAnalysisMultiplicityTask();
	void CreateDefaultTrackCuts ( Int_t year = 2010 );

	TH1D* nEvents;

	UInt_t fEventClassification;

	Double_t fConsistencyCut;
	
	Double_t fSPDPileupThreshold;
	
	virtual void CollectEvent ( const char* origin = "" );
	
	TH2D* hVertex;

	ClassDef ( AliAnalysisMultiplicityTask, 9 );
};

#endif // ALIANALYSISMULTIPLICITYTASK_H




