#ifndef ALIEMCALJETFINDERINPUT_H
#define ALIEMCALJETFINDERINPUT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  *  * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Initial input object for jetfinder
//
//*-- Author: Mark Horner (LBL/UCT)
//
//


#include "TObject.h"
#include "TParticle.h"
#include "TMCParticle6.h"
#include "TClonesArray.h"
#include "AliEMCALParton.h"
#include "AliEMCALDigit.h"
#include "AliEMCALJetFinderTypes.h"

class AliEMCALJetFinderInput : public TObject
{
	public:
		AliEMCALJetFinderInput();
		~AliEMCALJetFinderInput();
		void Reset(AliEMCALJetFinderResetType_t resettype);
		void SetDebug(Int_t debug=0){fDebug = debug;}
		void AddEnergyToDigit(Int_t digitID,Int_t denergy); 
		void AddTrack(TParticle track);
		void AddTrack(TMCParticle *track);
		void AddParton(AliEMCALParton *parton);
		void AddParticle(TParticle *particle);
		void AddParticle(TMCParticle *particle);
		AliEMCALDigit* GetDigit(Int_t digitID);
		Int_t GetNDigits() const {return fNDigits;}
		TParticle* GetTrack(Int_t trackID);
		Int_t GetNTracks() const {return fNTracks;}
		AliEMCALParton* GetParton(Int_t partonID);
		Int_t GetNPartons() const {return fNPartons;}
		TParticle* GetParticle(Int_t particleID);
		Int_t GetNParticles() const {return fNParticles;}

		AliEMCALJetFinderInput (const AliEMCALJetFinderInput&);
		AliEMCALJetFinderInput & operator = (const AliEMCALJetFinderInput &) {
		  Fatal("operator =", "not implemented") ;
		  return *this ;
		}

	private:
		void InitArrays();
		TClonesArray*	fDigitsArray;	//-> This is the digits array for the EMCAL
		Int_t		fNDigits;     	// This is the number of digits
		Int_t		fNMaxDigits;  	// This is the max number of digits
		TClonesArray*	fTracksArray; 	//-> This is the track array 
		Int_t		fNTracks;     	// This stores the number of tracks	
		Int_t		fNMaxTracks;	// This stores the maximum number of tracks
		TClonesArray*	fPartonsArray;  //->  This is the partons array
		Int_t		fNPartons;	// This stores the number of partons
		Int_t		fNMaxPartons;	// This stores the maximum number of partons
		TClonesArray*	fParticlesArray;//-> This stores the particles	
		Int_t		fNParticles;	// This stores the number of particles
	       	Int_t		fNMaxParticles; // This stroes the maximum number of particles
		Int_t		fDebug;		// This is the debug value 
		Bool_t		fInitialised;	// Stores whether or not the arrays have been initialised 
		
	ClassDef(AliEMCALJetFinderInput,5)
};
#endif
