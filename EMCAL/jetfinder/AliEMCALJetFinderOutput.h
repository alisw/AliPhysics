#ifndef ALIEMCALJETFINDEROUTPUT_H
#define ALIEMCALJETFINDEROUTPUT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  *  * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Output object for jetfinder
//
//*-- Author: 	Renan Cabrera (LBL)
//		Mark Horner (LBL/UCT)
//

class TClonesArray;
#include "TObject.h"
#include "TParticle.h"
#include "AliEMCALParton.h"
#include "AliEMCALJet.h"
#include "AliEMCALJetFinderTypes.h"

class AliEMCALJetFinderOutput : public TObject
{
		
	public:
		AliEMCALJetFinderOutput();
		~AliEMCALJetFinderOutput();
		void Reset(AliEMCALJetFinderResetType_t resettype);
		void AddJet(AliEMCALJet *jet); 
		void AddParton(AliEMCALParton *parton);
		void AddParticle(TParticle *particle);
		void SetDebug(Int_t debug){fDebug = debug;}
		AliEMCALJet* GetJet(Int_t jetID);
		Int_t GetNJets() const {return fNJets;}
		TClonesArray *GetJets() {return fJetsArray; }
		AliEMCALParton* GetParton(Int_t partonID);
		Int_t GetNPartons() const {return fNPartons;}
		TParticle* GetParticle(Int_t particleID);
		TClonesArray *GetParticles() {return fParticlesArray; }
		Int_t GetNParticles() const {return fNParticles;}

		AliEMCALJetFinderOutput (const AliEMCALJetFinderOutput&);
		AliEMCALJetFinderOutput & operator = (const AliEMCALJetFinderOutput & ) {
		  Fatal("operator =", "not implemented") ;
		  return *this ;
		}

	private:
		void InitArrays();
		TClonesArray	*fJetsArray;     	// Array of jet objects
		TClonesArray	*fPartonsArray;  	// Array of parton objects
		Int_t		fNPartons;		// Number of Partons actually stored
		Int_t		fNJets; 		// Number of jets actually stored
		TClonesArray    *fParticlesArray;	// Array of particles
		Int_t		fNParticles;		// Number of particles actually stored
                Int_t           fNMaxJets;      	// Maximum number of jets 
                Int_t           fNMaxParticles; 	// Maximum number of primary particles
                Int_t           fNMaxPartons;   	// Maximum number of primary particles
                Int_t           fDebug;			// Debug level
		Bool_t 		fInitialised;		// stores whether or not the arrays have been initialised

	ClassDef(AliEMCALJetFinderOutput,5)
		
};
#endif
