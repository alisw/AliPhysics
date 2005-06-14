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

/////////////////////////////////////////////////////////////////////////////////
//
// AliMUONSDigitizerv1 digitizes hits from the input stream into s-digits. 
// The s-digits are written to the s-digit tree TreeS.
// The same algorithm is implemented as for AliMUONDigitizerv1 however the
// chamber response is not applied to the s-digits and we write to TreeS and 
// not TreeD.
//
/////////////////////////////////////////////////////////////////////////////////

#include "AliMUONSDigitizerv1.h"
#include "AliMUON.h"
#include "AliMUONLoader.h"
#include "AliMUONConstants.h"
#include "AliMUONChamber.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONHit.h"
#include "AliMUONTransientDigit.h"
#include "AliLog.h"

ClassImp(AliMUONSDigitizerv1)

//___________________________________________
AliMUONSDigitizerv1::AliMUONSDigitizerv1() 
  : AliMUONDigitizer()
{
	// Default ctor - don't use it
}

//___________________________________________
AliMUONSDigitizerv1::AliMUONSDigitizerv1(AliRunDigitizer* manager) 
  : AliMUONDigitizer(manager)
{
	// ctor which should be used
}

//___________________________________________
AliMUONSDigitizerv1::~AliMUONSDigitizerv1()
{
	// Destructor
}

//-----------------------------------------------------------------------
void AliMUONSDigitizerv1::GenerateTransientDigits()
{
// Loops over all tracks and hits in the current selected event and calls 
// MakeTransientDigitsFromHit for each hit. 
// Note: Charge correlation is applied to the tracking chambers. 

	TTree* treeH = fGime->TreeH();
	AliDebug(2, Form("Generating transient digits using treeH = 0x%X"
			, (void*)treeH));
	//
	// Loop over tracks
	Int_t ntracks = (Int_t) treeH->GetEntries();
	for (Int_t itrack = 0; itrack < ntracks; itrack++) 
	{
		AliDebug(3, Form("Processing track %d...", itrack));
		fMUONData->ResetHits();
		treeH->GetEvent(itrack);
		//
		//  Loop over hits
		TClonesArray* hits = fMUONData->Hits();
		for (Int_t ihit = 0; ihit < hits->GetEntriesFast(); ihit++) 
		{
			AliMUONHit* mHit = static_cast<AliMUONHit*>( hits->At(ihit) );
			Int_t ichamber = mHit->Chamber()-1;  // chamber number
			if (ichamber > AliMUONConstants::NCh()-1) 
			{
				AliError(Form("Hit 0x%X has a invalid chamber number: %d", ichamber));
				continue;
			}
			//
			//Dumping Hit content:
			AliDebug(3,Form("Hit %d: chamber = %d\tX = %f\tY = %f\tZ = %f\teloss = %f",
					ihit, mHit->Chamber(), mHit->X(), mHit->Y(), mHit->Z(), mHit->Eloss()
				    ));
			// 
			// Inititializing Correlation
			AliMUONChamber& chamber = fMUON->Chamber(ichamber);
			chamber.ChargeCorrelationInit();
			if (ichamber < AliMUONConstants::NTrackingCh()) 
			{
				// Tracking Chamber
				// Initialize hit position (cursor) in the segmentation model 
			  
			  chamber.SigGenInit(mHit);

			} // else do nothing for Trigger Chambers
			
			MakeTransientDigitsFromHit(itrack, ihit, mHit);
		} // hit loop
	} // track loop      
}

//--------------------------------------------------------------------------
void AliMUONSDigitizerv1::MakeTransientDigitsFromHit(Int_t track, Int_t iHit, AliMUONHit * mHit)
{  
// This method is called for every hit in an event to generate AliMUONTransientDigits 
// from the hit and add these to fTDList.
// The AliMUONChamber::DisIntegration method us used to figure out which pads are 
// fired for a given hit. We then loop over the fired pads and add an AliMUONTransientDigit
// for each pad.

	AliDebug(4,Form("Making transient digit for hit number %d.", iHit));
		
	//
	// Calls the charge disintegration method of the current chamber 
	AliDebug(5,"Calling AliMUONChamber::DisIngtegration...");

	Float_t newdigit[6][500];  // Pad information
	Int_t nnew=0;              // Number of touched Pads per hit
	Int_t ichamber = mHit->Chamber()-1;
	AliMUONChamber& chamber = fMUON->Chamber(ichamber);

	chamber.DisIntegration(mHit, nnew, newdigit);

	// Creating new TransientDigits from hit
	for(Int_t iTD = 0; iTD < nnew; iTD++) 
	{
		Int_t charge;   
		Int_t digits[7];
		
		digits[0] = Int_t(newdigit[1][iTD]);  // Padx of the Digit
		digits[1] = Int_t(newdigit[2][iTD]);  // Pady of the Digit
		digits[2] = Int_t(newdigit[5][iTD]);  // Cathode plane
		digits[3] = Int_t(newdigit[3][iTD]);  // Induced charge in the Pad
		if (fSignal)
		{ 
			charge = digits[3];
			digits[4] = Int_t(newdigit[3][iTD]);  // Signal due to physics
		}
		else
		{
			charge = digits[3] + fMask;
			digits[4] = 0;    // No signal due to physics since this is now background.
		}
		digits[5] = iHit+fMask;    // Hit number in the list
		digits[6] =  mHit->DetElemId();

		AliDebug(5,Form("MakeTransientDigitsFromHit", 
				"DisIntegration result %d: PadX %d\tPadY %d\tPlane %d\tCharge %d\tHit %d\tidDE %d",
				iTD, digits[0], digits[1], digits[2], digits[3], digits[5], digits[6]));

		AliMUONTransientDigit* mTD = new AliMUONTransientDigit(ichamber, digits);
		mTD->AddToTrackList(track + fMask, charge);

		OnCreateTransientDigit(mTD, mHit);
		AddOrUpdateTransientDigit(mTD);
	}
}

//------------------------------------------------------------------------
void AliMUONSDigitizerv1::AddDigit(Int_t chamber, Int_t tracks[kMAXTRACKS], Int_t charges[kMAXTRACKS], Int_t digits[7])
{
// Derived to write to the s-digit tree TreeS.

	fMUONData->AddSDigit(chamber, tracks, charges, digits);   
}

//------------------------------------------------------------------------
Int_t AliMUONSDigitizerv1::GetSignalFrom(AliMUONTransientDigit* td)
{
// Returns the transient digit signal as is without applying the chamber response.
	AliDebug(4,"Returning TransientDigit signal.");
	return td->Signal(); 
}

//------------------------------------------------------------------------
Bool_t AliMUONSDigitizerv1::InitOutputData(AliMUONLoader* muonloader)
{
// Overridden to initialise the output tree to be TreeS rather than TreeD.
	AliDebug(3,"Creating s-digits branch and setting the tree address.");

	fMUONData->SetLoader(muonloader);

	// New branch per chamber for MUON digit in the tree of digits
	if (muonloader->TreeS() == NULL)
	{
		muonloader->MakeSDigitsContainer();
		if (muonloader->TreeS() == NULL)
		{
			AliError("Could not create TreeS.");
			return kFALSE;
		}
	}

	fMUONData->MakeBranch("S");
	fMUONData->SetTreeAddress("S");
	
	return kTRUE;
}

//------------------------------------------------------------------------
void AliMUONSDigitizerv1::FillOutputData()
{
// Overridden to fill TreeS rather than TreeD.

	AliDebug(3,"Filling trees with s-digits.");
	fMUONData->Fill("S");
	fMUONData->ResetSDigits();
}

//------------------------------------------------------------------------
void AliMUONSDigitizerv1::CleanupOutputData(AliMUONLoader* muonloader)
{
// Overridden to write and then cleanup TreeS that was initialised in InitOutputData.
	AliDebug(3,"Writing s-digits and releasing pointers.");
	muonloader->WriteSDigits("OVERWRITE");
	fMUONData->ResetSDigits();
	muonloader->UnloadSDigits();
}

//------------------------------------------------------------------------
Bool_t AliMUONSDigitizerv1::InitInputData(AliMUONLoader* muonloader)
{
// Derived to initialise the input to read from TreeH the hits tree. 
// If the hits are not loaded then we load the hits using the muon loader.

	AliDebug(3, "Loading hits in READ mode and setting the tree address.");

	fMUONData->SetLoader(muonloader);

	if (muonloader->TreeH() == NULL)
	{
		muonloader->LoadHits("READ");
		if (muonloader->TreeH() == NULL)
		{
			AliError("Can not load the hits tree.");
			return kFALSE;
		}
	}

	fMUONData->SetTreeAddress("H");
	return kTRUE;
}

//------------------------------------------------------------------------
void AliMUONSDigitizerv1::CleanupInputData(AliMUONLoader* muonloader)
{
// Derived to release the loaded hits and unload them.

	AliDebug(3, "Releasing loaded hits.");
	fMUONData->ResetHits();
	muonloader->UnloadHits();
}

