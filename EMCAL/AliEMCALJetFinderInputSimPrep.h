#ifndef ALIEMCALJETFINDERINPUTSIMPREP_H
#define ALIEMCALJETFINDERINPUTSIMPREP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  *  * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Class for JetFinder Input preparation from simulated data 
//*-- Author: Mark Horner (LBL/UCT)
//
//


#include "AliEMCALJetFinderInput.h"
#include "AliEMCALJetFinderInputPrep.h"
#include "TTask.h"
#include "AliEMCALJetFinderTypes.h"

class AliEMCALJetFinderInputSimPrep : public AliEMCALJetFinderInputPrep
{
	public:
	AliEMCALJetFinderInputSimPrep();
	~AliEMCALJetFinderInputSimPrep();
	void Reset(AliEMCALJetFinderResetType_t resettype);
        void SetEMCALType(AliEMCALJetFinderEMCALType_t emcaltype )  {fEMCALType = emcaltype;}
	//void SetDebug(Int_t debug = 0)  {fDebug = debug;}
        void SetSmearingType(AliEMCALJetFinderSmearingType_t smeartype ) {fSmearType = smeartype;}
	void SetTrackType(AliEMCALJetFinderTrackType_t tracktype){fTrackType = tracktype;}  
	void SetEfficiency(Float_t efficiency)  {fEfficiency = efficiency;}
	void SetTimeCut(Float_t timecut)  {fTimeCut = timecut; fEMCALType = kTimeCut;}
	Int_t FillFromFile(TString *filename, AliEMCALJetFinderFileType_t filetype,Int_t EventNumber,TString data);
	AliEMCALJetFinderInput* GetJetFinderInput(){return &fInputObject;}
	private:
	void FillHits();		// Fill from the hits to input object from simulation
	void FillSDigits();		// Fill from the hits to input object from simulation
	void FillTracks();		// Fill from particles simulating a TPC to input object from simulation
	void Smear(TParticle *particle);
	Bool_t Efficiency();
	void FillPartons();		// Fill partons to input object from simulation
	void FillPartonTracks(AliEMCALParton *parton);	// Fill partons to input object from simulation
	void FillParticles();		// Fill particles to input object from simulation
	void FillDigits();		// Fill digits to input object  

        AliEMCALJetFinderEMCALType_t 	fEMCALType;	// The EMCAL type set by the user	
        AliEMCALJetFinderSmearingType_t fSmearType;	// The efficiency and smearing for TPC	
        AliEMCALJetFinderTrackType_t 	fTrackType;	// The Track type set by the user	
	AliEMCALJetFinderFileType_t	fFileType;	//! The type of file being processed
	Float_t				fEfficiency;	// The TPC efficiency
	Float_t				fTimeCut;	// User specified time cut
	Float_t				fEtaMax;	// User specified time cut
	Float_t				fEtaMin;	// User specified time cut
	Float_t				fPhiMax;	// User specified time cut
	Float_t				fPhiMin;	// User specified time cut
	
	ClassDef(AliEMCALJetFinderInputSimPrep,1)
	
};
#endif
