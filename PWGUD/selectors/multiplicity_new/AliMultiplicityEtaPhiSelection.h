#ifndef ALIMULTIPLICITYETAPHISELECTION_H
#define ALIMULTIPLICITYETAPHISELECTION_H
#include "AliMultiplicityAnalysisSelection.h"
#include "AliPhiCut.h"
#include <TH2.h>

class AliMultiplicityEtaPhiSelection : public AliMultiplicityAnalysisSelection {

public:
	virtual Bool_t AcceptESDTrack ( const AliESDtrack* track, AliMultiplicityEventSelection::mEstimators tType, Bool_t complement );
	virtual Bool_t AcceptTracklet ( const AliMultiplicity* spdmultiplicity, Int_t iTracklet, AliMultiplicityEventSelection::mEstimators tType );
	
// 	virtual Bool_t AcceptESDTrackMC ( const AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t complement );
	virtual Bool_t AcceptTrackletMC ( const AliMultiplicity* spdmultiplicity, Int_t iTracklet, AliStack* , AliMultiplicityEventSelection::mEstimators tType );

	virtual Bool_t AcceptESDEvent ( const AliESDEvent* esd );
	virtual Bool_t AcceptESDandMCEvent ( const AliESDEvent* esd, AliMCEvent* , Bool_t );	
	
	virtual void CreateSelectionHistograms();
	virtual void Result();
	virtual Long64_t Merge ( const TCollection* mergelist );
	virtual void Conclude();
	virtual void SaveSelectionHistograms();

	AliMultiplicityEtaPhiSelection ( const char* name = "EtaPhiSelection", const char* title = "#eta vs #phi" );
	virtual ~AliMultiplicityEtaPhiSelection();
	TH2D* GetHistogram ( AliMultiplicityEventSelection::mEstimators type );
	
	TH2D* GetEtaChi2 ( Bool_t phiCut = kFALSE ) {
		Int_t i = 0; if ( phiCut ) i = 1;
		return hEtaChi2[i];
	};
	
	TH2D* GetBkg ( Bool_t phiCut = kFALSE ) {
		Int_t i = 0; if ( phiCut ) i = 1;
		return hBkg[i];
	}
	
	TH1D* GetBkgFraction ( Double_t chi2cut, Bool_t phiCut = kTRUE, Bool_t rebinned = kTRUE );

private:
	TH2D* hEtaPhi[kSPD];
	TH2D* hVzEta[9];
	TH2D* hVzNch[9];
	TH2D* hEtaChi2[2];

	Double_t fVz;
	Bool_t isClass[kSPD];
	TH1D* nClass;
	Int_t nTotal[kSPD];
	
	TH2D* hBkg[2];

	ClassDef ( AliMultiplicityEtaPhiSelection, 6 );
};

#endif // ALIMULTIPLICITYETAPHISELECTION_H
