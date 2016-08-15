#ifndef ALIDATAMULTIPLICITYSELECTION_H
#define ALIDATAMULTIPLICITYSELECTION_H

#include "AliMultiplicityEventSelection.h"
#include "AliMultiplicityAnalysisSelection.h"
#include <TH2.h>
#include <TProfile.h>

class AliMultiplicitySelection : public AliMultiplicityAnalysisSelection {

public:
	AliMultiplicitySelection ( const char* name = "multiplicity", const char* title = "Raw multiplicity" );

	virtual Bool_t AcceptTracklet ( const AliMultiplicity* spdmultiplicity, Int_t iTracklet, AliMultiplicityEventSelection::mEstimators tType );
	virtual Bool_t AcceptESDTrack ( AliESDtrack* track, AliMultiplicityEventSelection::mEstimators tType, Bool_t complement, Bool_t count = kTRUE );
	virtual Bool_t AcceptParticle ( TParticle* particle, const Int_t, Bool_t count = kTRUE );
	
	virtual Bool_t AcceptESDTrackMC ( AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t, Bool_t count = kTRUE );
	virtual Bool_t AcceptTrackletMC ( const AliMultiplicity* spdm, Int_t iTracklet, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t count = kTRUE );
	
	virtual Bool_t AcceptSecondaryESDTrack ( AliESDtrack* track, AliMultiplicityEventSelection::mEstimators tType, Bool_t complement, Bool_t count = kTRUE );

	virtual Bool_t AcceptESDEvent ( const AliESDEvent* esd, UInt_t mask );
	virtual Bool_t AcceptMCEvent ( const AliMCEvent* mc, UInt_t );

	virtual void CreateSelectionHistograms();
	virtual void Result();

	virtual Long64_t Merge ( const TCollection* mergelist );

	virtual void Conclude();
	virtual void ResetCounters();

	virtual void SaveSelectionHistograms();
	TH1D* GetMultiplicityHistogram ( mEstimators estimator );
	
	TH1D* GetMultiplicityHistogramMC ( ){
		return hMultiplicityMC;
	};
	TH1D* GetTrackMultiplicity();
	TH1D* GetPhiDist(Int_t i){ if ( i > -1 && i < 3 ) return hPhi[i]; else return 0; };

	virtual ~AliMultiplicitySelection() {};

	void SetTrackletChi2Cut ( Double_t cut = 1.6 );
	Double_t GetTrackletChi2Cut ( );
	Bool_t UseTrackletChi2Cut () { return useTrackletChi2Cut; };
	
	Double_t GetPileupFraction ( AliMultiplicityEventSelection::mEstimators estimator );
	TProfile* GetMeanPtvsMult ( );
	TH1D* GetTrackPt ( Int_t Ngen );
	
	virtual Bool_t CollectCondition(  ); 
	
// 	TH2D* GetZZ ( Int_t vtx = 0 ) {
// 		if ( ( vtx > -1 ) && ( vtx < 3 ) ) return hZrg[vtx];
// 		return 0;
// 	}
// 	
// 	TH2D* GetDZ ( Int_t vtx = 2 ) {
// 		if ( ( vtx > -1 ) && ( vtx < 3 ) ) return hDZg[vtx];
// 		return 0;
// 	}
	
protected:
	Int_t nGlobalTracks;
	Int_t nGlobalITScomplements;
	Int_t nGlobalTrackletComplements;
	Int_t nITSTrackletComplements;
	Int_t nITSTracks;
	Int_t nTracklets;

	Int_t nMCparticles;

	TH1D* hMultiplicity[kSPD];
	TH1D* hMultiplicityTracks;
	Bool_t wITSComplements;
	TH1D* hMultiplicityMC;
// 	Long64_t nNoVertexEvents;

	Bool_t isPileup;
	TH1D* hPileupvsMult[kSPD];
// 	TH1D* hTrackletChi2;

	Double_t trackletChi2Cut;
	Bool_t useTrackletChi2Cut;
	
	TH2D* fhPtvsNch;
	TList* tracks;                                             // !
	TH1D* hPileUpDist[2];
	
	TH2D* fhPtvsNchGen;
	TH2D* fhPtvsNchGenS[kSPD];
	TList* particles;                                          // !
	TList* particlesS[kSPD];                                   // !
	
	Int_t nGlobalSecondaryTracks;                              // !
	Int_t nGlobalSecondaryITScomplements;                      // !
	Int_t nITSSecondaryTracks;                                 // !
	
	TList* secondaryTracks;                                    // !
	
	TH1D* hEta[kSPD];
	TH2D* hTT;
	TH2D* hZSPD;
	
	TH1D* hPhi[kSPD];
	
private:
	
// 	TH2D* hZrg[3];
// 	TH2D* hDZg[3];
// 	TH1D* hEtarg[2];
// 	TH2D* hZetarg[2];
	
	ClassDef ( AliMultiplicitySelection, 30 );
};

#endif // ALIDATAMULTIPLICITYSELECTION_H
