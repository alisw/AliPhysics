#ifndef ALIMULTIPLICITYRESPONSESELECTION_H
#define ALIMULTIPLICITYRESPONSESELECTION_H

#include "AliMultiplicitySelection.h"
#include "AliMultiplicityHelper.h"
#include <TH2.h>


class AliMultiplicityResponseSelection : public AliMultiplicitySelection {

public:
	virtual void Conclude();
	virtual void SaveSelectionHistograms();
	virtual void Result();
	virtual void CreateSelectionHistograms();
	virtual Long64_t Merge ( const TCollection* mergelist );
	virtual Bool_t AcceptESDandMCEvent ( const AliESDEvent* esd, AliMCEvent* mc, Bool_t excludeESD, UInt_t mask );
	
	virtual Bool_t AcceptESDTrackMC ( AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t complement = kFALSE, Bool_t count = kTRUE );
	virtual Bool_t AcceptTrackletMC ( const AliMultiplicity* spdm, Int_t iTracklet, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t );
	
	virtual Bool_t AcceptTracklet ( const AliMultiplicity* , Int_t , AliMultiplicityEventSelection::mEstimators ){
		return kFALSE;
	};
	virtual Bool_t AcceptESDTrack ( AliESDtrack* , AliMultiplicityEventSelection::mEstimators , Bool_t, Bool_t ){
		return kFALSE;
	};
	
	virtual Bool_t AcceptParticle ( TParticle* particle, const Int_t, Bool_t );

	AliMultiplicityResponseSelection ( const char* name = "response", const char* title = "Detector response" );
	virtual ~AliMultiplicityResponseSelection();

	static TH2D* NormalizeResponse ( TH2D* response, TH1D* efficiency = 0 );

	TH2D* GetResponse ( mEstimators estimator );
	TH2D* GetResponseNormalized ( mEstimators estimator );

	TH1D* GetSeedDistribution();
	
	void SetPtWeights( TH1D* hweights, Bool_t RWSec = kFALSE );
	
	virtual Bool_t CollectCondition(  ); 
	
private:
	TH2D* response[kSPD];
// 	TH1D* hMultiplicityMC;
	TH1D* hMinus1binRedist;

	Bool_t fNoVertexESD;
	
	TH1D* hPtWeights;
	Bool_t fRWStrangeness;
	
	Int_t GetTrackWeight ( AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators, Bool_t );
	Int_t GetTrackWeight ( const AliMultiplicity* spdmultiplicity, Int_t iTracklet, AliStack* stack, AliMultiplicityEventSelection::mEstimators );
	Int_t GetParticleWeight ( TParticle* particle, const Int_t );
	
	ClassDef ( AliMultiplicityResponseSelection, 16 );
};

#endif // ALIMULTIPLICITYRESPONSESELECTION_H
