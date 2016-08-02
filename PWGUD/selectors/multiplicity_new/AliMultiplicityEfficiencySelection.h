#ifndef ALIMULTIPLICITYEFFICIENCYSELECTION_H
#define ALIMULTIPLICITYEFFICIENCYSELECTION_H
#include "AliMultiplicitySelection.h"
#include "AliAnalysisMultiplicityTask.h"
#include <TProfile.h>

class AliMultiplicityEfficiencySelection : public AliMultiplicitySelection {

public:
	virtual Bool_t AcceptESDandMCEvent ( const AliESDEvent* esd, AliMCEvent* mc, Bool_t excludeESD, UInt_t mask );
	virtual Long64_t Merge ( const TCollection* mergelist );

	virtual Bool_t AcceptESDTrackMC ( AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t, Bool_t );
	virtual Bool_t AcceptTrackletMC ( const AliMultiplicity* m, Int_t iTracklet, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t );

	virtual void Conclude();
	virtual void Result();

	virtual void CreateSelectionHistograms();
	virtual void SaveSelectionHistograms();

	AliMultiplicityEfficiencySelection ( const char* name = "eficiency", const char* title = "efficiency selection" );
	virtual ~AliMultiplicityEfficiencySelection();

	TH1D* GetEfficiency ( Int_t type );
	Double_t GetGlobalTriggeringEficiency ( Int_t type );
	TH1D* GetGeneratedMultiplicity ( Int_t type );
	TH1D* GetSelectedMultiplicity ( Int_t type );

	TH1D* GetAlphaSD ( ); //SD events fraction per multiplicity
	TH1D* GetAlphaDD ( ); //DD events fraction per multiplicity

	virtual void ResetCounters();

	void SetWeightDiffraction ( TH2D* hwSDs, TH2D* hwNDDDs );
	void SetWeightDiffraction ( Double_t diffMassCut ) {
		if ( diffMassCut > 0 ) {
			bWeightDiffraction = kTRUE;
			fDiffMassCut = diffMassCut;
		} else {
			AliWarning ( "Negative mass cut, ignoring." );
		}
	}
	Bool_t GetWeightDiffraction () {
		return bWeightDiffraction;
	};

	Double_t GetDiffCut() {
		return fDiffMassCut;
	}

	TCanvas* DrawC();

// 	TH1D* Get0bin ( Int_t itype = 0 ) {
// 		switch (itype) {
// 			case 0:
// 				return h0bin;
// 			case 1:
// 				return h0binZ;
// 			case 2:
// 				return h0binE;
// 			default:
// 				AliWarning ( "Wrong type" );
// 				return 0;
// 		}
// 	};

private:
	TH1D* hEventsGenerated[3];
	TH1D* hEventsSelected[3];

// 	TH2D* hEventsGeneratedN[3];

	Long64_t nSelected[3];
	Long64_t nGenerated[3];

	Int_t currentProcessType;	//!

	TProfile* pStrangeSurvival[3 * 2];

	Int_t nSurvived[3 * 2];
	Int_t nGeneratedS[2];

	TH1D* hSurvivedEta[3];
	TH1D* hGeneratedEta;

	Bool_t bIsGenerated;

	Bool_t bWeightDiffraction;
	Double_t fDiffMassCut;

	TH2D* hwSD;
	TH2D* hwNDDD;

// 	TH1D* hMZZ[2];
// 	TH1D* h0bin;
// 	TH1D* h0binZ;
// 	TH1D* h0binE;
//
// 	TH1D* hNZc[2];

	ClassDef ( AliMultiplicityEfficiencySelection, 12 );
};
#endif // ALIMULTIPLICITYEFFICIENCYSELECTION_H
