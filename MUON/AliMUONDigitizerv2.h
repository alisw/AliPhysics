#ifndef ALIMUONDIGITIZERV2_H
#define ALIMUONDIGITIZERV2_H
/* Copyright(c) 1998-2001, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliMUONDigitizerv1.h"

class AliMUONDigitizerv2 : public AliMUONDigitizerv1
{
  public:    
	AliMUONDigitizerv2();
	virtual ~AliMUONDigitizerv2();
	
	// Preferred constructor which assigns the manager object.
	AliMUONDigitizerv2(AliRunDigitizer * manager);

	// Makes a transient digit from the specified s-digit from the specified chamber.
	void MakeTransientDigitFromSDigit(Int_t iChamber, AliMUONDigit* sDigit);

	// The following methods are inherited from AliMUONDigitizerv1 and overridden.
	virtual void GenerateTransientDigits();
	virtual void AddDigit(Int_t chamber, Int_t tracks[kMAXTRACKS], Int_t charges[kMAXTRACKS], Int_t digits[6]);
	virtual Bool_t InitInputData(AliMUONLoader* muonloader);
	virtual void CleanupInputData(AliMUONLoader* muonloader);

	ClassDef(AliMUONDigitizerv2, 0) 
};    
#endif

