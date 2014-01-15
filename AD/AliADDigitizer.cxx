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
 
/* $Id: AliADDigitizer.cxx  $ */

///_________________________________________________________________________
///
/// This class constructs Digits out of Hits
///
///

// --- Standard library ---

// --- ROOT system ---
#include <TMath.h>
#include <TTree.h>
#include <TMap.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <AliGeomManager.h>
#include <TRandom.h>
#include <TF1.h>
#include <TH1F.h>

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliAD.h"
#include "AliADhit.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliGRPObject.h"
#include "AliDigitizationInput.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
//#include "AliADCalibData.h"
#include "AliCTPTimeParams.h"
#include "AliLHCClockPhase.h"
#include "AliADdigit.h"
#include "AliADDigitizer.h"
//#include "AliADSDigit.h"

ClassImp(AliADDigitizer)

//____________________________________________________________________________ 
AliADDigitizer::AliADDigitizer()
  :AliDigitizer(),
   fNdigits(0),
   fDigits(0)
{
  // default constructor
}

//____________________________________________________________________________ 
AliADDigitizer::AliADDigitizer(AliDigitizationInput* digInput)
  :AliDigitizer(digInput),
   fNdigits(0),
   fDigits(0)
{
  // constructor
}
           
//____________________________________________________________________________ 
AliADDigitizer::~AliADDigitizer()
{
  // destructor
  
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    fDigits=0; 
  }

}

//____________________________________________________________________________ 
Bool_t AliADDigitizer::Init()
{
	fDigits = new TClonesArray ("AliADdigit",1000);
	return kTRUE;
}
//____________________________________________________________________________ 

void AliADDigitizer::Digitize(Option_t* /*option*/) 
{   
   // Creates digits from hits
 
	// tem. variables
  Int_t modules[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t moduls[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t mods;
 

	// loaders

	AliRunLoader* outRunLoader = AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());
	if (!outRunLoader)
	{
		Error("Exec","Can not get output Run Loader");
		return;
	}
	AliLoader* outLoader = outRunLoader->GetLoader("ADLoader");
	if (!outLoader)
	{
		Error("Exec","Can not get output AD Loader");
		return;
	}
  
	outLoader->LoadDigits("update");
	if (!outLoader->TreeD()) outLoader->MakeTree("D");
	outLoader->MakeDigitsContainer();
	TTree* treeD = outLoader->TreeD();
	Int_t bufsize = 16000;
	treeD->Branch("ADdigit",&fDigits, bufsize);

	const Float_t eMin = 0;//1.52; //! MeVs, minimum energy
	// loop over inputs

	for (Int_t iInput=0; iInput < fDigInput->GetNinputs();iInput++)
	{
		AliRunLoader* runLoader = AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(iInput));
		AliLoader* loader = runLoader->GetLoader("ADLoader");
		if (!loader)
		{
			Error("Exec","Can not get AD Loader for input %d",iInput);
			continue;
		}
		if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
		AliAD* ad = (AliAD*) runLoader->GetAliRun()->GetDetector("AD");
		if (!ad)
		{
			Error("Exec","No AD detector for input %d",iInput);
			continue;
		}
		loader->LoadHits();
		TTree* treeH = loader->TreeH();
		if (!treeH)
		{
			Error("Exec","Cannot get TreeH for input %d",iInput);
			continue;
		}
		TClonesArray* hits = ad->Hits();

		// here I loop over tracks
		Int_t nTracks = (Int_t) treeH->GetEntries();
		for (Int_t iTrack=0; iTrack < nTracks; iTrack++)
		{
			ad->ResetHits();
			treeH->GetEvent(iTrack);
			Int_t nHits = hits->GetEntriesFast();
			// here comes the loop over AD hits
			for (Int_t iHit=0; iHit < nHits; iHit++)
			{
				AliADhit* hit = (AliADhit *)hits->UncheckedAt(iHit);
				Float_t eloss_mev = hit->GetEloss()*1000.0;
				Int_t module = hit->GetModule();
				//cout << "Module AD!!! " << module << endl;
				// at some point we shoukd have some calib objects to set real theresholds 
				// simple checking on hit, minimum energy at scintillator pad should be > 1.52 MeV's
				if (eloss_mev > eMin)
				{
				  modules[module] = 1;//cout << "energy: " << eloss_mev << endl;
				}else modules[module] = 0;
			}// end loop over hits
		} // endo loop over tracks
		for (Int_t i=0; i<16; i++)
		{
			moduls[i] = modules[i];
			//cout << "iModule: " << i <<  " AD hits: " << moduls[i] << endl;
		}

		loader->UnloadHits();

	} // end loop over inputs

	// here I add the hits to the TreeD

	Int_t tracks[3] = {-1,-1,-1};
	for (Int_t i=0; i<16; i++)
	{
		if (moduls[i]==1)
		{
			mods = i;
			AddDigit(tracks,mods,mods);
		}
	}
	treeD->Fill();
	outLoader->WriteDigits("OVERWRITE");
	outLoader->UnloadDigits();
	ResetDigit();
}

//____________________________________________________________________________
void AliADDigitizer::AddDigit(Int_t* track, Int_t module, Float_t cell) 
{  
   // Adds Digit 
   TClonesArray &ldigits = *fDigits;  
   new(ldigits[fNdigits++]) AliADdigit(track,module,cell);
}
//____________________________________________________________________________
void AliADDigitizer::AddDigit(Int_t* modul,Float_t cell) 
{  
   // Adds Digit 
   TClonesArray &ldigits = *fDigits;  
   new(ldigits[fNdigits++]) AliADdigit(modul,cell);
}

//____________________________________________________________________________
void AliADDigitizer::ResetDigit()
{
   // Clears Digits
   fNdigits = 0;
   if (fDigits) fDigits->Delete();
}

