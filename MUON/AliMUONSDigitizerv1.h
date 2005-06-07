#ifndef ALIMUONSDIGITIZERV1_H
#define ALIMUONSDIGITIZERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

// The AliMUONSDigitizer produces
// SDigits from Hits 
// J.P Cussonneau Subatech Feb 2004

#include "AliMUONDigitizer.h"

class AliMUONHit;

class AliMUONSDigitizerv1 : public AliMUONDigitizer
{
  public:    
	AliMUONSDigitizerv1();
	virtual ~AliMUONSDigitizerv1();
	
	// Preferred constructor to call which sets the manager.
	AliMUONSDigitizerv1(AliRunDigitizer * manager);

	// methods from old AliMUONDigitizerv1
	void MakeTransientDigitsFromHit(Int_t itrack, Int_t ihit, AliMUONHit * mHit);
	void GenerateTransientDigits();

	void AddDigit(Int_t chamber, Int_t tracks[kMAXTRACKS], Int_t charges[kMAXTRACKS], Int_t digits[7]);
	Int_t GetSignalFrom(AliMUONTransientDigit* td);
	Bool_t InitOutputData(AliMUONLoader* muonloader);
	void FillOutputData();
	void CleanupOutputData(AliMUONLoader* muonloader);

	// methods from old AliMUONDigitizerv1
	virtual Bool_t InitInputData(AliMUONLoader* muonloader);
	virtual void CleanupInputData(AliMUONLoader* muonloader);

	// to disable trigger in SDigitizer
	void CreateTrigger(){return;}
	Bool_t FetchTriggerPointer(AliMUONLoader* /*loader*/ ){return kTRUE;}
	void CleanupTriggerArrays(){return;}
	void FillTriggerOutput(){return;}
	void AddDigitTrigger(Int_t /*chamber*/, Int_t* /*tracks[kMAXTRACKS]*/, 
			     Int_t* /*charges[kMAXTRACKS]*/, Int_t* /*digits[6]*/,
                             const Int_t /*digitindex*/
		) {return;}

	ClassDef(AliMUONSDigitizerv1, 0)
};    
#endif

