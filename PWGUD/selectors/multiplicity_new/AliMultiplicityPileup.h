#ifndef ALIMULTIPLICITYPILEUP_H
#define ALIMULTIPLICITYPILEUP_H
#include "AliMultiplicitySelection.h"
#include "AliAnalysisMultiplicityTask.h"
#include <TProfile.h>

class AliMultiplicityPileup : public AliMultiplicitySelection {

public:
	virtual Bool_t AcceptESDandMCEvent ( const AliESDEvent* esd, AliMCEvent* mc, Bool_t excludeESD, UInt_t mask );
	virtual Long64_t Merge ( const TCollection* mergelist );

	virtual Bool_t AcceptESDTrackMC ( const AliESDtrack* track, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType, Bool_t complement );
	virtual Bool_t AcceptTrackletMC ( const AliMultiplicity* spdm, Int_t iTracklet, AliStack* stack, AliMultiplicityEventSelection::mEstimators tType );
	virtual Bool_t AcceptParticle ( const TParticle* particle, const Int_t iParticle = -1 );

	virtual void Conclude();
	virtual void Result();

	virtual void CreateSelectionHistograms();
	virtual void SaveSelectionHistograms();

	AliMultiplicityPileup ( const char* name = "plup", const char* title = "pileup" );
	virtual ~AliMultiplicityPileup();

	virtual void ResetCounters();

	virtual Int_t FindHeader ( const Int_t iParticle );
	
	//getters
	TH2D* GetResponse(Int_t iMethod, Bool_t bTrue);
	TH2D* GetExtraTvsZ(Int_t iMethod);
	TH2D* GetExtraTvsZwo0(Int_t iMethod);
	TH1D* GetGenerated1();
	
	Double_t GetAffectedFraction(Int_t iMethod);

private:

	// multiplicity histograms
	TH1D* hGenerated;	//
	TH1D* hSelected[3 * 2];	//
	
// 	TH1D* hSampleBias[4];	//

	TH2D* hResponse[3];	//
	TH2D* hResponse_true[3];	//

	//other histograms
	TH2D* hdZvsExtra[3];	//
	TH2D* hdZvsExtraNC[3];	//
	TH1D* hNHeaders;		//
	
	TH2D* hVertexCorrelation[2];		//
	
// 	TH1D* hSPDPileupEfficiency[3];		//
	
// 	TH2D* hFakeRate;					//
	

	//counters
	// 0 - coinciding vertex
	// 1 - others
	Int_t nParticlesPrimary;				//
	Int_t nParticlesSecondary;			//
	Int_t nSelectedGlobal[3];			//
	Int_t nSelectedITSSA[3];			//
	Int_t nSelectedTracklets[3];			//
	Int_t nSelectedGlobalITSSAcomplements[3];	//
	Int_t nSelectedGlobalTrackletcomplements[3];	//
	Int_t nSelectedITSSATrackletComplements[3];	//

	Int_t nHeaders;					//

	//internal cache
	TList* lHeaders;	//
	Int_t nCoincidence;	//
	Int_t nClosest;		//
	Double_t ZESD;		//
	Double_t dZESD;		//
	Double_t dZ;		//

	ClassDef ( AliMultiplicityPileup, 5 );
};
#endif // ALIMULTIPLICITYPILEUP_H
