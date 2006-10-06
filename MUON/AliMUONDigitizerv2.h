#ifndef ALIMUONDIGITIZERV2_H
#define ALIMUONDIGITIZERV2_H
/* Copyright(c) 1998-2001, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup sim
/// \class AliMUONDigitizerv2
/// \brief Class produces Digits from SDigits.
///
/// Do the Digitization (Digit) from summable Digits (SDigit)
/// Allow the merging of signal file with background file(s).

#include "AliMUONDigitizer.h"
#include "AliMUONDigit.h"

class AliMUONLoader;

class AliMUONDigitizerv2 : public AliMUONDigitizer
{
  public:    
	AliMUONDigitizerv2();
	virtual ~AliMUONDigitizerv2();
	
	// Preferred constructor which assigns the manager object.
	AliMUONDigitizerv2(AliRunDigitizer * manager);

	// Makes a transient digit from the specified s-digit from the specified chamber.
	void MakeTransientDigitFromSDigit(Int_t iChamber, AliMUONDigit* sDigit);

	// The following methods are inherited from AliMUONDigitizerv1 and overridden.
	void GenerateTransientDigits();
	void AddDigit(Int_t chamber, Int_t tracks[kMAXTRACKS], Int_t charges[kMAXTRACKS], Int_t digits[7]);
	Bool_t InitInputData(AliMUONLoader* muonloader);
	void CleanupInputData(AliMUONLoader* muonloader);

	// methods from old AliMUONDigitizerv1
	Int_t GetSignalFrom(AliMUONTransientDigit* td);
	Bool_t InitOutputData(AliMUONLoader* muonloader);
	void CleanupOutputData(AliMUONLoader* muonloader);
	void FillOutputData();

	ClassDef(AliMUONDigitizerv2, 0) 
};    
#endif

