//-*- Mode: C++ -*-
#ifndef ALIGAMMACONVERSIONAODBGHANDLER_H
#define ALIGAMMACONVERSIONAODBGHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class  for handling of background calculation
//---------------------------------------------
////////////////////////////////////////////////

#include <vector>


// --- ROOT system ---
#include <TObject.h> 
#include "AliKFParticle.h"
#include "AliAODConversionPhoton.h"
#include "AliAODConversionMother.h"
#include "TClonesArray.h"
#include "AliESDVertex.h"

typedef std::vector<AliAODConversionPhoton*> AliGammaConversionAODVector;
typedef std::vector<AliAODConversionMother*> AliGammaConversionMotherAODVector;

class AliGammaConversionAODBGHandler : public TObject {

	public: 
	struct GammaConversionVertex{
		Double_t fX;
		Double_t fY;
		Double_t fZ;
		Double_t fEP;
	};
	
	typedef struct GammaConversionVertex GammaConversionVertex; 																//!

	typedef std::vector<AliGammaConversionAODVector> AliGammaConversionBGEventVector;
	typedef std::vector<AliGammaConversionBGEventVector> AliGammaConversionMultipicityVector;
	typedef std::vector<AliGammaConversionMultipicityVector> AliGammaConversionBGVector;

	typedef std::vector<AliGammaConversionMotherAODVector> AliGammaConversionMotherBGEventVector;
	typedef std::vector<AliGammaConversionMotherBGEventVector> AliGammaConversionMotherMultipicityVector;
	typedef std::vector<AliGammaConversionMotherMultipicityVector> AliGammaConversionMotherBGVector;
	

	AliGammaConversionAODBGHandler();																							//constructor
    AliGammaConversionAODBGHandler(Int_t binsZ,Int_t binsMultiplicity,Int_t nEvents);										// constructor
    AliGammaConversionAODBGHandler(Int_t collisionSystem,Int_t centMin,Int_t centMax,Int_t nEvents, Bool_t useTrackMult, Int_t mode,Int_t binsZ, Int_t binsMultiplicity);
	AliGammaConversionAODBGHandler(const AliGammaConversionAODBGHandler & g);													//copy constructor
	AliGammaConversionAODBGHandler & operator = (const AliGammaConversionAODBGHandler & g);										//assignment operator
	virtual ~AliGammaConversionAODBGHandler();																					//virtual destructor

	void Initialize(Double_t * const zBinLimitsArray, Double_t * const multiplicityBinLimitsArray);

	Int_t GetZBinIndex(Double_t z) const;

	Int_t GetMultiplicityBinIndex(Int_t mult) const;

	Int_t GetNBackgroundEventsInBuffer(Int_t binz, int binMult) const;

	void AddEvent(TList* const eventGammas, Double_t xvalue,Double_t yvalue,Double_t zvalue, Int_t multiplicity, Double_t epvalue = -100);
	void AddMesonEvent(TList* const eventMothers, Double_t xvalue,Double_t yvalue,Double_t zvalue, Int_t multiplicity, Double_t epvalue = -100);
	void AddMesonEvent(const std::vector<AliAODConversionMother> &eventMother, Double_t xvalue, Double_t yvalue, Double_t zvalue, Int_t multiplicity, Double_t epvalue = -100);
	void AddElectronEvent(TClonesArray* const eventENeg, Double_t zvalue, Int_t multiplicity);

	Int_t GetNBGEvents()const {return fNEvents;}

	// Get BG photons
	AliGammaConversionAODVector* GetBGGoodV0s(Int_t zbin, Int_t mbin, Int_t event);
	
	// Get BG mesons
	AliGammaConversionMotherAODVector* GetBGGoodMesons(Int_t zbin, Int_t mbin, Int_t event);
	
	// Get BG electron
	AliGammaConversionAODVector* GetBGGoodENeg(Int_t event, Double_t zvalue, Int_t multiplicity);
	
	void PrintBGArray();

	GammaConversionVertex * GetBGEventVertex(Int_t zbin, Int_t mbin, Int_t event){return &fBGEventVertex[zbin][mbin][event];}

	Double_t GetBGProb(Int_t z, Int_t m){return fBGProbability[z][m];}

	private:

		Int_t 								fNEvents; 						// number of events
		Int_t ** 							fBGEventCounter;				//! bg counter
		Int_t ** 							fBGEventENegCounter;			//! bg electron counter
		Int_t ** 							fBGEventMesonCounter;			//! bg counter
		Double_t ** 						fBGProbability; 				//! prob per bin
		GammaConversionVertex *** 			fBGEventVertex;					//! array of event vertex
		Int_t 								fNBinsZ;	 					//n z bins
		Int_t 								fNBinsMultiplicity; 			//n bins multiplicity
		Double_t *							fBinLimitsArrayZ;				//! bin limits z array
		Double_t *							fBinLimitsArrayMultiplicity;	//! bin limit multiplicity array
		AliGammaConversionBGVector 			fBGEvents; 						// photon background events
		AliGammaConversionBGVector 			fBGEventsENeg; 					// electron background electron events
		AliGammaConversionMotherBGVector 	fBGEventsMeson; 				// neutral meson background events
		
	ClassDef(AliGammaConversionAODBGHandler,5)
};
#endif
