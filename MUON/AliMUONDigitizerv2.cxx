/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//  Do the Digitization (Digit) from summable Digits (SDigit)
//  Allow the merging of signal file with background file(s).

/////////////////////////////////////////////////////////////////////////////////
//
// AliMUONDigitizerv2 digitizes digits from s-digits. 
// Merging is allowed from multiple input streams. The first input stream is 
// assumed to be the signal and all other input streams as assumed to be 
// background.
// Chamber response is applied to the digits identically to AliMUONDigitizerv1.
//
/////////////////////////////////////////////////////////////////////////////////

#include "AliMUONData.h"
#include "AliMUONLoader.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONDigitizerv2.h"
#include "AliMUONTransientDigit.h"
#include "AliLog.h"

ClassImp(AliMUONDigitizerv2)

//___________________________________________
AliMUONDigitizerv2::AliMUONDigitizerv2() : AliMUONDigitizerv1()
{
	// Default ctor - don't use it
}

//___________________________________________
AliMUONDigitizerv2::AliMUONDigitizerv2(AliRunDigitizer* manager) : AliMUONDigitizerv1(manager)
{
	// ctor which should be used
}

//___________________________________________
AliMUONDigitizerv2::~AliMUONDigitizerv2()
{
	// Destructor
}

//-----------------------------------------------------------------------
void AliMUONDigitizerv2::GenerateTransientDigits()
{
// Loop over all chambers and s-digits in the input stream and create 
// AliMUONTransientDigit objects from them. These are then added to fTDList.

	AliDebug(2,"Generating transient digits using treeH = 0x%X");
	//
	// Loop over SDigits
	Int_t ndig, k;
	AliMUONDigit* sDigit;
	TClonesArray* muonSDigits;
	for (Int_t ich = 0; ich < AliMUONConstants::NCh(); ich++)  // loop over chamber
	{
		fMUONData->ResetSDigits();
		fMUONData->GetCathodeS(0);
		muonSDigits = fMUONData->SDigits(ich); 
		ndig = muonSDigits->GetEntriesFast();
		for (k = 0; k < ndig; k++)
		{
			sDigit = (AliMUONDigit*) muonSDigits->UncheckedAt(k);
			MakeTransientDigitFromSDigit(ich,sDigit);
		}
		fMUONData->ResetSDigits();
		fMUONData->GetCathodeS(1);
		muonSDigits = fMUONData->SDigits(ich); 
		ndig=muonSDigits->GetEntriesFast();
		for (k = 0; k < ndig; k++)
		{
			sDigit = (AliMUONDigit*) muonSDigits->UncheckedAt(k);
			MakeTransientDigitFromSDigit(ich,sDigit);
		};
	} // SDigits loop, end loop over chamber
};

//------------------------------------------------------------------------
void AliMUONDigitizerv2::MakeTransientDigitFromSDigit(Int_t iChamber, AliMUONDigit* sDigit)
{
// Makes a transient digit from the specified s-digit from the specified chamber. 
// Once the digit is created it is added to the fTDList.

	AliDebug(4,Form("Making transient digit from s-digit for chamber %d.", iChamber));
	Int_t digits[6];
	//
	// Creating a new TransientDigits from SDigit
	digits[0] = sDigit->PadX();  // Padx of the Digit
	digits[1] = sDigit->PadY();  // Pady of the Digit
	digits[2] = sDigit->Cathode()+1;  // Cathode plane
	digits[3] = sDigit->Signal();  // Induced charge in the Pad
	if (fSignal)
		digits[4] = sDigit->Signal();
	else
		digits[4] = 0;
		
	digits[5] = sDigit->Hit();    // Hit number in the list

	AliDebug(5,Form("Made digit from sDigit 0x%X: PadX %d\tPadY %d\tPlane %d\tCharge %d\tHit %d",
			(void*)sDigit, digits[0], digits[1], digits[2], digits[3], digits[5]));
	
	AliMUONTransientDigit* mTD = new AliMUONTransientDigit(iChamber, digits);
	// Copy list of tracks and trackcharge
	for(Int_t itrack = 0; itrack < kMAXTRACKS; itrack++)
	{
		Int_t track = sDigit->Track(itrack);
		if (track < 0) break;  // Check if we reached the end of the track list.
		mTD->AddToTrackList( track + fMask, sDigit->TrackCharge(itrack) );
	};

	OnCreateTransientDigit(mTD, sDigit);
	AddOrUpdateTransientDigit(mTD);
};

//------------------------------------------------------------------------
void AliMUONDigitizerv2::AddDigit(Int_t chamber, Int_t tracks[kMAXTRACKS], Int_t charges[kMAXTRACKS], Int_t digits[6])
{
// Override to add new digits to the digits tree TreeD.
	fMUONData->AddDigit(chamber, tracks, charges, digits);   
};

//------------------------------------------------------------------------
Bool_t AliMUONDigitizerv2::InitInputData(AliMUONLoader* muonloader)
{
// Overridden to initialize fMUONData to read from the s-digits tree TreeS. 
// If the s-digits are not loaded then the muon loader is used to load the
// s-digits into memory.

	AliDebug(3,"Loading s-digits in READ mode and setting the tree address.");
	fMUONData->SetLoader(muonloader);

	if (muonloader->TreeS() == NULL)
	{
		muonloader->LoadSDigits("READ");
		if (muonloader->TreeS() == NULL)
		{
			AliError("Can not load the s-digits tree.");
			return kFALSE;
		};
	};

	fMUONData->SetTreeAddress("S");
	return kTRUE;
};

//------------------------------------------------------------------------
void AliMUONDigitizerv2::CleanupInputData(AliMUONLoader* muonloader)
{
// Overridden to release and unload s-digits from memory.

	AliDebug(3,"Releasing loaded s-digits.");
	fMUONData->ResetSDigits();
	muonloader->UnloadSDigits();
};
