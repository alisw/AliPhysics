#ifndef ALIMULTIPLICITYHELPER_H
#define ALIMULTIPLICITYHELPER_H

#include <Rtypes.h>
#include <TObject.h>
#include <AliVertex.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliMCEvent.h>
#include <AliTriggerAnalysis.h>
#include <TMath.h>
#include <AliGenEventHeader.h>
#include <AliStack.h>

// struct stat;
// enum MCProcessType {													// event generator process types
// 	kND = 0,															// generic inellastic process
// 	kSingleDiffractive = 1,												// single-diffractive process
// 	kDoubleDiffractive = 2,												// double-diffractive process
// 	kNonDiffractive = -1,												// non-diffractive process
// 	kUndefined = -2,													// undefined process
// 	kNonSingleDiffractive = 3,											// non-single difractive
// };

// inel = knd + ksingle + kdouble

class AliMultiplicityHelper : public TObject {

public:

	enum VertexType {													// vertex types
	    kVertexTracks, 													// track vertex
	    kVertexSPD, 													// pixel vertex (default)
	    kVertexMC,														// MC vertex
	    kPrimaryVertex													// best available data vertex
	};

	enum EClassification {
	    kVGlobal = BIT ( 0 ),
	    kVSPD = BIT ( 1 ),
	    kVSPD1D = BIT ( 2 ),
	    kVSPD3D = BIT ( 3 ),
	    kVConsistent = BIT ( 4 ),
	    kTriggered = BIT ( 5 ),
	    kMBOR = BIT ( 6 ),
	    kV0AND = BIT ( 7 ),
	    kNonDiffractive = BIT ( 8 ),											// non-diffractive process
	    kSingleDiffractive = BIT ( 9 ),										// single-diffractive process
	    kDoubleDiffractive = BIT ( 10 ),										// double-diffractive process
	    kSPDPileUp = BIT ( 11 ),
	    kVInconsistent = BIT ( 12 ),
	    kVNoSPD = BIT ( 13 ),
	    kVNoSPDDC = BIT ( 14 ),
	    kVNoGlobal = BIT ( 15 ),
	    kSPDNoPileup = BIT ( 16 ),
	    kHighMult = BIT ( 17 ),
	    kNotClusterCutBG = BIT ( 18 )
	};
	// inel = knd + ksingle + kdouble


	AliMultiplicityHelper();
	virtual ~AliMultiplicityHelper();

	static const AliVertex* GetVertex ( VertexType vType, const AliVEvent* event );		// get pointer to event vertex
	static AliGenEventHeader* GetGenEventHeader ( const AliVEvent* mc );				// get mc event header

	static Bool_t IsQualityVertex ( const AliVertex* v, Bool_t& disp );					// check vertex quality
	static Bool_t IsQualityVertex ( const AliVEvent* event );							// check vertex quality

	static Double_t GetTrackletChi2 ( const AliMultiplicity* m, Int_t iTracklet );		// computes tracklet chi2
	static Int_t FindPrimaryMother ( AliStack* stack, Int_t label );
	static TParticle* GetMother ( AliStack* stack, Int_t label );

	static Bool_t CheckDiff ( AliStack* stack, Int_t& gen_type, Int_t& weighted_type, Double_t energyCMS, TH2D* hwSD, TH2D* hwNDDD );	//diffraction re-weighting
	static Bool_t CheckDiff ( AliStack* stack, Int_t& gen_type, Int_t& weighted_type, Double_t energyCMS, Double_t MMax );
	static Double_t GetDiffractiveMass ( AliStack* stack, Double_t energyCMS );
	static Bool_t IsCocktail ( AliGenEventHeader* h );

	static UInt_t Classify ( AliESDEvent* esd, const AliMCEvent* mc, Double_t consistencyCut = 0.5, Double_t SPDPileUpThreshold = 0.8 );
	static void Declassify ( UInt_t mask );

private:
	static Bool_t GetWeights ( Double_t mass, Double_t& minMass, Double_t& maxMass, Double_t& wSD, Double_t& wDD, Double_t& wND, Int_t energyCMS, TH2D* hwSD, TH2D* hwNDDD );	//get weights for diffraction from supplied histogram
	ClassDef ( AliMultiplicityHelper, 8 );

};

#endif // ALIMULTIPLICITYHELPER_H
