#ifndef ALIMUONDIGITIZERV1_H
#define ALIMUONDIGITIZERV1_H
/* Copyright(c) 1998-2001, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// The AliMUONDigitizer procees :
// - Addition of hits from different tracks
// - Merging of hits from different files
// - The response function of the chamber.
// - Simulation of the electronic noise, threshold and saturation
// 
// Gines MARTINEZ Subatech Feb 2003 

#include "AliMUONDigitizer.h"

class AliMUONDigitizerv1 : public AliMUONDigitizer 
{
  public:
	AliMUONDigitizerv1();
	virtual ~AliMUONDigitizerv1();
	
	// Preferred constructor which assigns the manager object.
	AliMUONDigitizerv1(AliRunDigitizer * manager);
    
  protected:
	// Generation of a TransientDigits from a hit object.
	void MakeTransientDigitsFromHit(Int_t itrack, Int_t ihit, AliMUONHit * mHit);
	
	// The following methods are all derived from AliMUONDigitizer
	virtual void GenerateTransientDigits();
	virtual void AddDigit(Int_t chamber, Int_t tracks[kMAXTRACKS], Int_t charges[kMAXTRACKS], Int_t digits[6]);
	virtual Int_t GetSignalFrom(AliMUONTransientDigit* td);
	virtual Bool_t InitOutputData(AliMUONLoader* muonloader);
	virtual void FillOutputData();
	virtual void CleanupOutputData(AliMUONLoader* muonloader);
	virtual Bool_t InitInputData(AliMUONLoader* muonloader);
	virtual void CleanupInputData(AliMUONLoader* muonloader);
   
	ClassDef(AliMUONDigitizerv1, 2)
};    
#endif

