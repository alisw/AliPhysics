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
#include "AliMUONLoader.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONTransientDigit.h"

ClassImp(AliMUONSDigitizerv1)

//___________________________________________
AliMUONSDigitizerv1::AliMUONSDigitizerv1() 
  : AliMUONDigitizerv1()
{
	// Default ctor - don't use it
}

//___________________________________________
AliMUONSDigitizerv1::AliMUONSDigitizerv1(AliRunDigitizer* manager) 
  : AliMUONDigitizerv1(manager)
{
	// ctor which should be used
}

//___________________________________________
AliMUONSDigitizerv1::~AliMUONSDigitizerv1()
{
	// Destructor
}

//------------------------------------------------------------------------
void AliMUONSDigitizerv1::AddDigit(Int_t chamber, Int_t tracks[kMAXTRACKS], Int_t charges[kMAXTRACKS], Int_t digits[6])
{
// Derived to write to the s-digit tree TreeS.

	fMUONData->AddSDigit(chamber, tracks, charges, digits);   
};

//------------------------------------------------------------------------
Int_t AliMUONSDigitizerv1::GetSignalFrom(AliMUONTransientDigit* td)
{
// Returns the transient digit signal as is without applying the chamber response.

	if (GetDebug() > 3) Info("GetSignalFrom", "Returning TransientDigit signal.");
	return td->Signal(); 
};

//------------------------------------------------------------------------
Bool_t AliMUONSDigitizerv1::InitOutputData(AliMUONLoader* muonloader)
{
// Overridden to initialise the output tree to be TreeS rather than TreeD.

	if (GetDebug() > 2)
		Info("InitOutputData", "Creating s-digits branch and setting the tree address.");

	fMUONData->SetLoader(muonloader);

	// New branch per chamber for MUON digit in the tree of digits
	if (muonloader->TreeS() == NULL)
	{
		muonloader->MakeSDigitsContainer();
		if (muonloader->TreeS() == NULL)
		{
			Error("InitOutputData", "Could not create TreeS.");
			return kFALSE;
		};
	};

	fMUONData->MakeBranch("S");
	fMUONData->SetTreeAddress("S");
	
	return kTRUE;
};

//------------------------------------------------------------------------
void AliMUONSDigitizerv1::FillOutputData()
{
// Overridden to fill TreeS rather than TreeD.

	if (GetDebug() > 2) Info("FillOutputData", "Filling trees with s-digits.");
	fMUONData->Fill("S");
	fMUONData->ResetSDigits();
};

//------------------------------------------------------------------------
void AliMUONSDigitizerv1::CleanupOutputData(AliMUONLoader* muonloader)
{
// Overridden to write and then cleanup TreeS that was initialised in InitOutputData.

	if (GetDebug() > 2) Info("CleanupOutputData", "Writing s-digits and releasing pointers.");
	muonloader->WriteSDigits("OVERWRITE");
	fMUONData->ResetSDigits();
	muonloader->UnloadSDigits();
};
