
#include <Riostream.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TTree.h> 
#include <TMath.h>

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONSDigitizerv1.h"
#include "AliMUONHit.h"
#include "AliMUONHitMapA1.h"
#include "AliMUONPadHit.h"
#include "AliMUONTransientDigit.h"

/////////////////////////////////////////////////////////////////////////////////
//
// AliMUONSDigitizerv1 digitizes hits from the input stream into s-digits. 
// The s-digits are written to the s-digit tree TreeS.
// The same algorithm is implemented as for AliMUONDigitizerv1 however the
// chamber response is not applied to the s-digits and we write to TreeS and 
// not TreeD.
//
/////////////////////////////////////////////////////////////////////////////////

ClassImp(AliMUONSDigitizerv1)

//___________________________________________
AliMUONSDigitizerv1::AliMUONSDigitizerv1() : AliMUONDigitizerv1()
{
	// Default ctor - don't use it
}

//___________________________________________
AliMUONSDigitizerv1::AliMUONSDigitizerv1(AliRunDigitizer* manager) : AliMUONDigitizerv1(manager)
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

	muondata->AddSDigit(chamber, tracks, charges, digits);   
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

	muondata->SetLoader(muonloader);

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

	muondata->MakeBranch("S");
	muondata->SetTreeAddress("S");
	
	return kTRUE;
};

//------------------------------------------------------------------------
void AliMUONSDigitizerv1::FillOutputData()
{
// Overridden to fill TreeS rather than TreeD.

	if (GetDebug() > 2) Info("FillOutputData", "Filling trees with s-digits.");
	muondata->Fill("S");
	muondata->ResetSDigits();
};

//------------------------------------------------------------------------
void AliMUONSDigitizerv1::CleanupOutputData(AliMUONLoader* muonloader)
{
// Overridden to write and then cleanup TreeS that was initialised in InitOutputData.

	if (GetDebug() > 2) Info("CleanupOutputData", "Writing s-digits and releasing pointers.");
	muonloader->WriteSDigits("OVERWRITE");
	muondata->ResetSDigits();
	muonloader->UnloadSDigits();
};
