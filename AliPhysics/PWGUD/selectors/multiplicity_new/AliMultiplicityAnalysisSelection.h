#ifndef ALIMULTIPLICITYANALYSISSELECTION_H
#define ALIMULTIPLICITYANALYSISSELECTION_H

#include "AliMultiplicityEventSelection.h"
#include "AliPhiCut.h"
#include <AliStack.h>


class AliMultiplicityAnalysisSelection : public AliMultiplicityEventSelection {

public:
	virtual Bool_t AcceptESDTrack ( AliESDtrack*, AliMultiplicityEventSelection::mEstimators, Bool_t, Bool_t count = kTRUE ) {
		return count;
	};
	virtual Bool_t AcceptSecondaryESDTrack ( AliESDtrack*, AliMultiplicityEventSelection::mEstimators, Bool_t, Bool_t count = kTRUE ) {
		return count;
	};
	virtual Bool_t AcceptSecondaryESDTrackMC( AliESDtrack* , AliStack* , AliMultiplicityEventSelection::mEstimators, Bool_t, Bool_t count = kTRUE ) {
		return count;
	};
	virtual Bool_t AcceptTracklet ( const AliMultiplicity*, Int_t, AliMultiplicityEventSelection::mEstimators ) {
		return kTRUE;
	};
	virtual Bool_t AcceptParticle ( TParticle*, const Int_t, Bool_t count = kTRUE ) {
		return count;
	};

	virtual Bool_t AcceptESDTrackMC ( AliESDtrack*, AliStack*, AliMultiplicityEventSelection::mEstimators, Bool_t complement = kFALSE, Bool_t count = kTRUE ) {
		return complement && count;
	};
	virtual Bool_t AcceptTrackletMC ( const AliMultiplicity*, Int_t, AliStack*, AliMultiplicityEventSelection::mEstimators, Bool_t count = kTRUE ) {
		return count;
	};

	AliMultiplicityAnalysisSelection ( const char* name, const char* title );
	virtual ~AliMultiplicityAnalysisSelection() {};

	virtual void Conclude() = 0;
	virtual void ResetAccepted();
	virtual void ResetCounters() {};
	virtual void Add0bin ( const AliESDEvent* ) {};

	Bool_t IsVertexInCuts ( const AliESDEvent* esd, UInt_t mask );
	Bool_t IsVertexInCuts ( const AliMCEvent* mc );

	Double_t GetEtaCut();
	Double_t GetZvCut();
	Double_t GetZvCutTS( Int_t side = 0 );

	void SetEtaCut ( Double_t etaCut );
	void SetZvCut ( Double_t ZvCut );
	void SetZvCutTS ( Double_t zlow, Double_t zhi );
	
	virtual void SetQueryPhysSelection ( Bool_t query );
	virtual Bool_t GetQueryPhysSelection();


	void SetCorrelateTracks ( Bool_t corr );
	Bool_t GetCorrelateTracks();
	
	virtual void SetPhiCut( Int_t energy = 0 ) 
	{
		fPhiCut = new AliPhiCut(energy);
	};
	
	Bool_t IsPreselected() {
		return fPreSelected;
	}

protected:
	Double_t fEtaCut;
	Double_t fZvCut;
	Double_t fZvCutTS[2];
	Bool_t useZcutTS;
	
	Bool_t fCorrelateTracks;

	Bool_t checkIfPhysSelected;
	
	AliPhiCut* fPhiCut;
	
	Bool_t fSPDInCuts;
	Bool_t fGlobalInCuts;
	Bool_t fMCVInCuts;
	
	Bool_t fPreSelected;	//!

	ClassDef ( AliMultiplicityAnalysisSelection, 12 );
};

#endif // ALIMULTIPLICITYANALYSISSELECTION_H
