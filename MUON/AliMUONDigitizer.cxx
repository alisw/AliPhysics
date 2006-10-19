/**************************************************************************
 * Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliMUONDigitizer.h"
#include "AliMUONConstants.h"
#include "AliMUONSegmentation.h"
#include "AliMUONHitMapA1.h"
#include "AliMUON.h"
#include "AliMUONLoader.h"
#include "AliMUONConstants.h"
#include "AliMUONTransientDigit.h"

#include "AliMpDEIterator.h"

#include "AliLog.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"

/////////////////////////////////////////////////////////////////////////////////////
//
// AliMUONDigitizer should be base abstract class of all digitisers in the MUON 
// module. It implements the common functionality of looping over input streams
// filling the fTDList and writing the fTDList to the output stream. 
// Inheriting digitizers need to override certain methods to choose and initialize
// the correct input and output trees, apply the correct detector response if any
// and implement how the transient digits are generated from the input stream.
//
/////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliMUONDigitizer)

//___________________________________________
AliMUONDigitizer::AliMUONDigitizer() : 
	AliDigitizer(),
	fRunLoader(0),
	fGime(0),
	fMUON(0),
	fMUONData(0),
	fHitMap(0),
	fTDList(0),
	fTDCounter(0),
	fMask(0),
	fSignal(0),
	fSegmentation(0)
{
/// Default constructor.
/// Initializes all pointers to NULL.
}

//___________________________________________
AliMUONDigitizer::AliMUONDigitizer(AliRunDigitizer* manager) : 
	AliDigitizer(manager),
	fRunLoader(0),
	fGime(0),
	fMUON(0),
	fMUONData(0),
	fHitMap(0),
	fTDList(0),
	fTDCounter(0),
	fMask(0),
	fSignal(0),
	fSegmentation(0)
{
/// Constructor which should be used rather than the default constructor.
/// Initializes all pointers to NULL.
}

//___________________________________________
AliMUONDigitizer::~AliMUONDigitizer() 
{
/// Destructor

  if (fMUONData)
    delete fMUONData;

}

//------------------------------------------------------------------------
Bool_t AliMUONDigitizer::Init()
{
/// Initialize - Does nothing. 

  return kTRUE;
}

//------------------------------------------------------------------------
void AliMUONDigitizer::Exec(Option_t* /*option*/)
{
/// The main work loop starts here. 
/// The digitization process is broken up into two steps: 
/// 1) Loop over input streams and create transient digits from the input.
///    Done in GenerateTransientDigits()
/// 2) Loop over the generated transient digits and write them to the output
///    stream. Done in CreateDigits()

	AliDebug(1, "Running digitiser.");

	if (fManager->GetNinputs() == 0)
	{
		AliWarning("No inputs set, nothing to do.");
		return;
	}

	if (!FetchLoaders(fManager->GetInputFolderName(0), fRunLoader, fGime) ) return;
	if (! FetchGlobalPointers(fRunLoader) ) return;

	InitArrays();
	
	AliDebug(2, Form("Event Number is %d.", fManager->GetOutputEventNr()));

	// Loop over files to merge and to digitize
	fSignal = kTRUE;
	for (Int_t inputFile = 0; inputFile < fManager->GetNinputs(); inputFile++)
	{
		fMask = fManager->GetMask(inputFile);
		AliDebug(2, Form("Digitising folder %d, with fMask = %d: %s", inputFile, fMask,
			     (const char*)fManager->GetInputFolderName(inputFile)));

		if (inputFile != 0)
			// If this is the first file then we already have the loaders loaded.
			if (! FetchLoaders(fManager->GetInputFolderName(inputFile), fRunLoader, fGime) )
				continue;
		else
			// If this is not the first file then it is assumed to be background.
			fSignal = kFALSE;

		if (! InitInputData(fGime) ) continue;
		GenerateTransientDigits();
		CleanupInputData(fGime);
	}

	Bool_t ok = FetchLoaders(fManager->GetOutputFolderName(), fRunLoader, fGime);
	if (ok) ok = InitOutputData(fGime);
	if (ok) CreateDigits();
	if (ok) CreateTrigger();
	if (ok) CleanupOutputData(fGime);

	CleanupArrays();
	CleanupTriggerArrays();
}

//--------------------------------------------------------------------------
void AliMUONDigitizer::AddOrUpdateTransientDigit(AliMUONTransientDigit* mTD)
{
/// Checks to see if the transient digit exists in the corresponding fHitMap.
/// If it does then the digit is updated otherwise it is added. 

	if (ExistTransientDigit(mTD)) 
	{
		UpdateTransientDigit(mTD);
		delete mTD;  // The new digit can be deleted because it was not added.
	}
	else 
		AddTransientDigit(mTD);
}

//------------------------------------------------------------------------
void AliMUONDigitizer::UpdateTransientDigit(AliMUONTransientDigit* mTD)
{
/// Update the transient digit that is already in the fTDList by adding the new
/// transient digits charges and track lists to the existing one.

	AliDebug(4,Form( "Updating transient digit 0x%X", (void*)mTD));
	// Choosing the maping of the cathode plane of the chamber:
	Int_t detElemId =  mTD->DetElemId();

	Int_t iNchCpl= fNDetElemId[detElemId] + (mTD->Cathode()-1) * AliMUONConstants::NDetElem();

	AliMUONTransientDigit *pdigit = 
		static_cast<AliMUONTransientDigit*>( fHitMap[iNchCpl]->GetHit(mTD->PadX(),mTD->PadY()) );

	// update charge
	pdigit->AddSignal( mTD->Signal() );
	pdigit->AddPhysicsSignal( mTD->Physics() );  
	
	// update list of tracks
	Int_t ntracks = mTD->GetNTracks();
	if (ntracks > kMAXTRACKS)  // Truncate the number of tracks to kMAXTRACKS if we have to.
	{
		AliDebug(1,Form( 
			"TransientDigit returned the number of tracks to be %d, which is bigger than kMAXTRACKS.",
			ntracks));
		AliDebug(1,Form( "Reseting the number of tracks to be %d.", kMAXTRACKS));
		ntracks = kMAXTRACKS;
	}
	
	for (Int_t i = 0; i < ntracks; i++)
	{
		pdigit->UpdateTrackList( mTD->GetTrack(i), mTD->GetCharge(i) );
	}
}

//------------------------------------------------------------------------
void AliMUONDigitizer::AddTransientDigit(AliMUONTransientDigit* mTD)
{
/// Adds the transient digit to the fTDList and sets the appropriate entry
/// in the fHitMap arrays. 

	AliDebug(4,Form( "Adding transient digit 0x%X", (void*)mTD));
	// Choosing the maping of the cathode plane of the chamber:

	Int_t detElemId =  mTD->DetElemId();
	Int_t iNchCpl= fNDetElemId[detElemId] + (mTD->Cathode()-1) * AliMUONConstants::NDetElem();

	fTDList->AddAtAndExpand(mTD, fTDCounter);
	if (iNchCpl>-1 && iNchCpl<2*AliMUONConstants::NDetElem()) {
	  fHitMap[iNchCpl]->SetHit( mTD->PadX(), mTD->PadY(), fTDCounter);
	  fTDCounter++;
	}
}

//------------------------------------------------------------------------
Bool_t AliMUONDigitizer::ExistTransientDigit(AliMUONTransientDigit* mTD)
{
/// Checks if the transient digit already exists on the corresponding fHitMap.
/// i.e. is there a transient digit on the same chamber, cathode and pad position 
/// as mTD. If yes then kTRUE is returned else kFASLE is returned.

	// Choosing the maping of the cathode plane of the chamber:
	Int_t detElemId =  mTD->DetElemId();

	Int_t iNchCpl= fNDetElemId[detElemId] + (mTD->Cathode()-1) *AliMUONConstants::NDetElem() ;

	//	Int_t iNchCpl= mTD->Chamber() + (mTD->Cathode()-1) * AliMUONConstants::NCh();
	if (iNchCpl>-1 && iNchCpl<2*AliMUONConstants::NDetElem())
	  return( fHitMap[iNchCpl]->TestHit(mTD->PadX(), mTD->PadY()) );
	else return kFALSE;
}

//-----------------------------------------------------------------------
void AliMUONDigitizer::CreateDigits()
{
/// Loops over the fTDList for each cathode, gets the correct signal for the 
/// digit and adds the new digit to the output stream.

	fTDList->Sort(); // sort by idDE
        AliDebug(2, "Creating digits...");
	//	for (Int_t icat = 0; icat < 2; icat++) {

          Int_t digitindex[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	  //
	  // Filling Digit List
	  Int_t nentries = fTDList->GetEntriesFast();
	  for (Int_t nent = 0; nent < nentries; nent++) {

	    AliMUONTransientDigit* td = (AliMUONTransientDigit*)fTDList->At(nent);
	    if (td == NULL) continue; 
			
	    // Must be the same cathode, otherwise we will fill a mixture
	    // of digits from both cathodes.
	    //if (icat != td->Cathode() - 1) continue;
			
	    AliDebug(3,Form( "Creating digit from transient digit 0x%X", (void*)td));

	    Int_t q = GetSignalFrom(td);
	    if (q > 0) {
	      Int_t chamber = td->Chamber();
	      if (0 <= chamber && chamber <= 13 )
		AddDigit(td, q, digitindex[chamber]++);
	      else
		AliError(Form("Invalid chamber %d\n",chamber));
	    }
	  }
	  FillOutputData();
	  //	}
	  fTDCounter = 0;
}

//------------------------------------------------------------------------
void AliMUONDigitizer::AddDigit(
		AliMUONTransientDigit* td, Int_t responseCharge,
		Int_t digitindex
	)
{
/// Prepares the digits, track and charge arrays in preparation for a call to
/// AddDigit(Int_t, Int_t[kMAXTRACKS], Int_t[kMAXTRACKS], Int_t[6])
/// This method is called by CreateDigits() whenever a new digit needs to be added
/// to the output stream trees.
/// The responseCharge value is used as the Signal of the new digit. 
/// The OnWriteTransientDigit method is also called just before the adding the 
/// digit to allow inheriting digitizers to be able to do some specific processing 
/// at this point. 

	Int_t tracks[kMAXTRACKS];
	Int_t charges[kMAXTRACKS];
	Int_t digits[7];
      
	digits[0] = td->PadX();
	digits[1] = td->PadY();
	digits[2] = td->Cathode() - 1;
	digits[3] = responseCharge;
	digits[4] = td->Physics();
	digits[5] = td->Hit();
	digits[6] =  td->DetElemId();

	Int_t nptracks = td->GetNTracks();
	if (nptracks > kMAXTRACKS) {

	  AliDebug(1, Form(
			   "TransientDigit returned the number of tracks to be %d, which is bigger than kMAXTRACKS.",
			   nptracks));
	  AliDebug(1, Form("Reseting the number of tracks to be %d.", kMAXTRACKS));
	  nptracks = kMAXTRACKS;
	}
	
	for (Int_t i = 0; i < nptracks; i++) {

	  tracks[i]   = td->GetTrack(i);
	  charges[i]  = td->GetCharge(i);
	}

	// Sort list of tracks according to charge
	SortTracks(tracks,charges,nptracks);

	if (nptracks < kMAXTRACKS ) {

	  for (Int_t i = nptracks; i < kMAXTRACKS; i++) {
	    tracks[i]  = -1;
	    charges[i] = 0;
	  }
	}

	AliDebug(4,Form( "Adding digit with charge %d.", responseCharge));

	OnWriteTransientDigit(td);
	AddDigit(td->Chamber(), tracks, charges, digits);
	AddDigitTrigger(td->Chamber(), tracks, charges, digits, digitindex);
}

//------------------------------------------------------------------------
void AliMUONDigitizer::OnCreateTransientDigit(AliMUONTransientDigit* /*digit*/, TObject* /*source_object*/)
{
/// Does nothing.
///
/// This is derived by Digitisers that want to trace which digits were made from
/// which hits.
}

//------------------------------------------------------------------------
void AliMUONDigitizer::OnWriteTransientDigit(AliMUONTransientDigit* /*digit*/)
{
/// Does nothing.
///
/// This is derived by Digitisers that want to trace which digits were made from
/// which hits.
}

//------------------------------------------------------------------------
Bool_t AliMUONDigitizer::FetchLoaders(const char* foldername, AliRunLoader*& runloader, AliMUONLoader*& muonloader)
{
/// Fetches the run loader from the current folder, specified by 'foldername'. 
/// The muon loader is then loaded from the fetched run loader. 
/// kTRUE is returned if no error occurred otherwise kFALSE is returned. 

	AliDebug(3, Form("Fetching run loader and muon loader from folder: %s", foldername));

	runloader = AliRunLoader::GetRunLoader(foldername);
	if (runloader == NULL)
	{
		AliError(Form("RunLoader not found in folder: %s", foldername));
		return kFALSE; 
	}                                                                                       
	muonloader = (AliMUONLoader*) runloader->GetLoader("MUONLoader");
	if (muonloader == NULL) 
	{
		AliError(Form("MUONLoader not found in folder: %s", foldername));
		return kFALSE; 
	}
	return kTRUE;

}

//------------------------------------------------------------------------
Bool_t AliMUONDigitizer::FetchGlobalPointers(AliRunLoader* runloader)
{
/// Fetches the AliRun object into the global gAlice pointer from the specified
/// run loader. The AliRun object is loaded into memory using the run loader if 
/// not yet loaded. The MUON module object is then loaded from gAlice and 
/// AliMUONData fetched from the MUON module. 
/// kTRUE is returned if no error occurred otherwise kFALSE is returned. 

        AliDebug(3, Form("Fetching gAlice, MUON module and AliMUONData from runloader 0x%X.",
			(void*)runloader
		    ));

	if (runloader->GetAliRun() == NULL) runloader->LoadgAlice();
	gAlice = runloader->GetAliRun();
	if (gAlice == NULL)
	{
		AliError(Form("Could not find the AliRun object in runloader 0x%X.", (void*)runloader));
		return kFALSE;
	}
	fMUON = (AliMUON*) gAlice->GetDetector("MUON");
	if (fMUON == NULL)
	{
		AliError(Form("Could not find the MUON module in runloader 0x%X.", (void*)runloader));
		return kFALSE;
	}

	AliMUONLoader *muonloader = (AliMUONLoader*) runloader->GetLoader("MUONLoader");
	if (muonloader == NULL) 
	{
		AliError( "MUONLoader not found ");
		return kFALSE; 
	}


	if (fMUONData == NULL) fMUONData = new AliMUONData(muonloader,"MUON","MUON");
	if (fMUONData == NULL)
	{
		AliError(Form("Could not find AliMUONData object in runloader 0x%X.", (void*)runloader));
		return kFALSE;
	}

	return kTRUE;
}

//------------------------------------------------------------------------
void AliMUONDigitizer::InitArrays()
{
/// Creates a new fTDList object. 
/// Also creates an array of 2 * chamber_number AliMUONHitMapA1 objects
/// in the fHitMaps array. Each one is set to a chamber and cathode
/// specific segmentation model. 
///
/// Note: the fTDList and fHitMap arrays must be NULL before calling this method.

    AliDebug(2, "Initialising internal arrays.");
    AliDebug(4, "Creating transient digits list.");
    fTDList = new TObjArray;
	
    // Array of pointer of the AliMUONHitMapA1:
    //  two HitMaps per chamber, or one HitMap per cahtode plane
    fHitMap = new AliMUONHitMapA1* [2*AliMUONConstants::NDetElem()];
    for (Int_t i=0; i<2*AliMUONConstants::NDetElem(); i++) fHitMap[i] = 0x0;

    Int_t k = 0;
    Int_t idDE;

    for (Int_t i = 0; i < AliMUONConstants::NCh(); i++) {


      AliDebug(4,Form( "Creating hit map for chamber %d, cathode 1.", i+1));

      AliMUONSegmentation* segmentation = fMUON->GetSegmentation();

      // Loop over detection elements
      AliMpDEIterator it;
      for ( it.First(i+1); !it.IsDone(); it.Next() ) { 
        
	idDE = it.CurrentDE();
	fNDetElemId[idDE] = k;

        AliMUONGeometrySegmentation* c1Segmentation 
          = segmentation->GetModuleSegmentationByDEId(idDE, 0); // Cathode plane 1
        Int_t npx1 = c1Segmentation->Npx(idDE)+1;
        Int_t npy1 = c1Segmentation->Npy(idDE)+1;
	fHitMap[k] = new AliMUONHitMapA1(npx1, npy1, fTDList); 
     
        AliMUONGeometrySegmentation* c2Segmentation 
          = segmentation->GetModuleSegmentationByDEId(idDE, 1); // Cathode plane 1
        Int_t npx2 = c2Segmentation->Npx(idDE)+1;
        Int_t npy2 = c2Segmentation->Npy(idDE)+1;
	fHitMap[k+AliMUONConstants::NDetElem()] = new AliMUONHitMapA1(npx2, npy2, fTDList);
	k++;
      }
    }
}
//------------------------------------------------------------------------
void AliMUONDigitizer::CleanupArrays()
{
/// The arrays fTDList and fHitMap are deleted and the pointers set to NULL.

	AliDebug(2, "Deleting internal arrays.");
	for(Int_t i = 0; i < 2*AliMUONConstants::NDetElem(); i++) {
		delete fHitMap[i];
	}
	delete [] fHitMap;
	fHitMap = NULL;
	
	AliDebug(4, "Deleting transient digits list.");
	fTDList->Delete();
	delete fTDList;
	fTDList = NULL;

}

//------------------------------------------------------------------------
void AliMUONDigitizer::SortTracks(Int_t *tracks, Int_t *charges, Int_t ntr) const
{
/// Sort the list of tracks contributing to a given digit
/// Only the 3 most significant tracks are actually sorted

       if (ntr <= 1) return;

       //
       //  Loop over signals, only 3 times
       //

       Int_t qmax;
       Int_t jmax;
       Int_t idx[3] = {-2,-2,-2};
       Int_t jch[3] = {-2,-2,-2};
       Int_t jtr[3] = {-2,-2,-2};
       Int_t i, j, imax;

       if (ntr < 3) imax = ntr;
       else imax=3;
	
       for(i = 0; i < imax; i++)
	 {
	   qmax=0;
	   jmax=0;

	   for(j = 0; j < ntr; j++)
	     {
	       if (	(i == 1 && j == idx[i-1]) || 
			(i == 2 && (j == idx[i-1] || j == idx[i-2]))
			) 
		 continue;

	       if(charges[j] > qmax) 
		 {
		   qmax = charges[j];
		   jmax = j;
		 }       
	     } 

	   if(qmax > 0)
	     {
	       idx[i] = jmax;
	       jch[i] = charges[jmax]; 
	       jtr[i] = tracks[jmax]; 
	     }

	 } 

       for(i = 0; i < 3; i++)
	 {
	   if (jtr[i] == -2) 
	     {
	       charges[i] = 0;
	       tracks[i] = 0;
	     } 
	   else 
	     {
	       charges[i] = jch[i];
	       tracks[i] = jtr[i];
	     }
	 }
}
