#ifndef ALIEMCALJETFINDERINPUTPREP_H
#define ALIEMCALJETFINDERINPUTPREP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  *  * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Base class for jetfinder input preparation 
//
//*-- Author: Mark Horner (LBL/UCT)
//



#include "TTask.h"
#include "AliEMCALJetFinderInput.h"
#include "AliEMCALJetFinderTypes.h"

class AliEMCALJetFinderInputPrep : public TTask
{
		
	public:
	AliEMCALJetFinderInputPrep();
	~AliEMCALJetFinderInputPrep();
	void Reset(AliEMCALJetFinderResetType_t resettype);
	void SetDebug(Int_t debug = 0){fDebug = debug; fInputObject.SetDebug(debug-2); }
	Int_t FillFromFile(TString *filename, AliEMCALJetFinderFileType_t filetype, Int_t EventNumber);
	AliEMCALJetFinderInput* GetJetFinderInput()  {return &fInputObject;}
	void SetPythiaComparison(Bool_t value) {fPythiaComparison = value;}
	Bool_t GetPythiaComparison() {return fPythiaComparison;}
	protected:
	Int_t 		fDebug; // The debug flag to be used for messages
	AliEMCALJetFinderInput 	fInputObject;	// The JetFinder input object to be filled
	Bool_t		fPythiaComparison;	// Special flag for pyclus comparison
	private:

	ClassDef(AliEMCALJetFinderInputPrep,3)
	
};
#endif
