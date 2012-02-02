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

//====================================================================================================================================================
//
//      Digitizer class for the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliDigitizationInput.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "AliMFTDigitizer.h"
#include "AliMFTDigit.h"
#include "AliMFT.h"
#include "AliMFTSegmentation.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "AliDigitizer.h"

ClassImp(AliMFTDigitizer)

//====================================================================================================================================================
    
AliMFTDigitizer::AliMFTDigitizer():
  AliDigitizer(),
  fNPlanes(0),
  fSegmentation(0)
{

  // default constructor 

}

//====================================================================================================================================================

AliMFTDigitizer::AliMFTDigitizer(AliDigitizationInput *digInp) :
  AliDigitizer(digInp),
  fNPlanes(0),
  fSegmentation(0)
{

}

//====================================================================================================================================================

void AliMFTDigitizer::Digitize(Option_t*) {

  // This method is responsible for merging sdigits to a list of digits

  AliDebug(1, "************************************************************************");
  AliDebug(1, "************************ AliMFTDigitizer::Exec *************************");
  AliDebug(1, "************************************************************************");
  
  if (!fSegmentation) {
    fSegmentation = new AliMFTSegmentation("AliMFTGeometry.root");
    fNPlanes = fSegmentation -> GetNPlanes();
  }

  AliDebug(1, Form("nPlanes = %d",fNPlanes));

  AliDebug(1,Form("Start with %i input(s) for event %i", fDigInput->GetNinputs(), fDigInput->GetOutputEventNr()));
    
  AliRunLoader *pInRunLoader=0;
  AliLoader    *pInMFTLoader=0;
  
  TClonesArray sDigits[fNMaxPlanes];
  for (Int_t iPlane=0; iPlane<fNMaxPlanes; iPlane++) sDigits[iPlane].SetClass("AliMFTDigit");  // tmp storage for sdigits sum up from all input files
  
  // filling the arrays of sdigits...

  for (Int_t iFile=0; iFile<fDigInput->GetNinputs(); iFile++) {
    
    pInRunLoader = AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(iFile));
    pInMFTLoader = pInRunLoader->GetLoader("MFTLoader"); 
    if (!pInMFTLoader) {
      AliDebug(1,"no MFT lodader, checking in the other input \n"); 
      continue;
    }
    
    if (!pInRunLoader->GetAliRun()) pInRunLoader->LoadgAlice();
    AliMFT *pInMFT = (AliMFT*) pInRunLoader->GetAliRun()->GetDetector("MFT"); 
    
    AliDebug(1, "pInMFTLoader->LoadSDigits()...");
    pInMFTLoader->LoadSDigits();
    AliDebug(1, "... done!");
    AliDebug(1, "    pInMFTLoader->TreeS()->GetEntry(0);");
    pInMFTLoader->TreeS()->GetEntry(0);
    AliDebug(1, "... done!");
    
    for (Int_t iPlane=0; iPlane<pInMFT->GetSDigitsList()->GetEntries(); iPlane++) {      
      for(Int_t iSDig=0; iSDig<((TClonesArray*)pInMFT->GetSDigitsList()->At(iPlane))->GetEntries(); iSDig++) {
	AliDebug(2, Form("Reading digit %03d of plane %02d (A)", iSDig, iPlane));
	AliMFTDigit *pSDig = (AliMFTDigit*) ((TClonesArray*)pInMFT->GetSDigitsList()->At(iPlane))->At(iSDig);
	AliDebug(2, Form("Reading digit %03d of plane %02d (B)", iSDig, iPlane));
	pSDig->AddOffset2TrackID(fDigInput->GetMask(iFile));   // -> To be introduced for merging (since all inputs count tracks independently from 0)
	AliDebug(2, Form("Reading digit %03d of plane %02d (C)", iSDig, iPlane));
	new ((sDigits[iPlane])[sDigits[iPlane].GetEntries()]) AliMFTDigit(*pSDig);  
	AliDebug(2, Form("Reading digit %03d of plane %02d (D)", iSDig, iPlane));
      }
    }
    
    pInMFTLoader->UnloadSDigits();   
    pInMFT->ResetSDigits();

  }
  
  AliRunLoader *pOutRunLoader = AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());  // open output stream (only 1 possible)
  AliLoader    *pOutMFTLoader = pOutRunLoader->GetLoader("MFTLoader");                        
  AliRun       *pAliRun       = pOutRunLoader->GetAliRun();
  AliMFT       *pOutMFT       = (AliMFT*) pAliRun->GetDetector("MFT");      
  pOutMFTLoader->MakeTree("D");   
  pOutMFT->MakeBranch("D");
  pOutMFT->SetTreeAddress();
  
  SDigits2Digits(sDigits, pOutMFT->GetDigitsList());   // here the sdigits are merged into digits
  
  pOutMFTLoader->TreeD()->Fill();              // fill the output tree with the list of digits
  pOutMFTLoader->WriteDigits("OVERWRITE");     // serialize them to file
  
  for (Int_t iPlane=0; iPlane<fNPlanes; iPlane++) sDigits[iPlane].Clear();   // remove all tmp sdigits
  pOutMFTLoader->UnloadDigits();   
  pOutMFT->ResetDigits(); 

}

//====================================================================================================================================================

void AliMFTDigitizer::SDigits2Digits(TClonesArray *pSDigitList, TObjArray *pDigitList) {   

  TClonesArray *myDigitList[fNMaxPlanes] = {0};
 
  for (Int_t iPlane=0; iPlane<fNPlanes; iPlane++) { 
    myDigitList[iPlane] = (TClonesArray*)(*pDigitList)[iPlane];
    if (myDigitList[iPlane]->GetEntries()!=0) AliErrorClass("Some of digits lists is not empty");   // in principle those lists should be empty 
  }
   
  AliDebug(1,"starting loop over planes");
   
  for (Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {
     
    AliMFTDigit *newDig=NULL;
    AliMFTDigit *oldDig=NULL;
    pSDigitList[iPlane].Sort();

    AliDebug(1,"starting loop over sdigits to create digits");
    AliDebug(1, Form("MFT plane #%02d has %d SDigits", iPlane, Int_t(pSDigitList[iPlane].GetEntries())));

    for (Int_t iSDig=0; iSDig<pSDigitList[iPlane].GetEntries(); iSDig++) {

      newDig = (AliMFTDigit*) (pSDigitList[iPlane].At(iSDig));
      Bool_t digitExists = kFALSE;
      Int_t nDigits = myDigitList[iPlane]->GetEntries();
      
      for (Int_t iDig=0; iDig<nDigits; iDig++) {
	oldDig = (AliMFTDigit*) (myDigitList[iPlane]->At(iDig));
	if (newDig->GetDetElemID()==oldDig->GetDetElemID() &&
	    newDig->GetPixelX()==oldDig->GetPixelX() &&
	    newDig->GetPixelY()==oldDig->GetPixelY() &&
	    newDig->GetPixelZ()==oldDig->GetPixelZ()) {
	  digitExists = kTRUE;
	  MergeDigits(oldDig, newDig);
	  break;
	}
      }

      if (!digitExists) new ((*myDigitList[iPlane])[myDigitList[iPlane]->GetEntries()]) AliMFTDigit(*newDig);

    }
    
    AliDebug(1, Form("MFT plane #%02d has %d Digits", iPlane, Int_t(myDigitList[iPlane]->GetEntries())));
    
//     for (Int_t iDigit=0; iDigit<myDigitList[iPlane]->GetEntries(); iDigit++) {
//       AliDebug(1, Form("Digit %03d of MFT plane %02d has pixel coordinates (%05d, %05d)", 
// 		       iDigit, iPlane, ((AliMFTDigit*) myDigitList[iPlane]->At(iDigit))->GetPixelX(), ((AliMFTDigit*) myDigitList[iPlane]->At(iDigit))->GetPixelY()) );
//     }

    AliDebug(1, "ending loop over sdigits to create digits");

  }

  AliDebug(1,"ending loop over layers");  

}

//====================================================================================================================================================

void AliMFTDigitizer::MergeDigits(AliMFTDigit *mainDig, AliMFTDigit *digToSum) {
  
  mainDig -> SetEloss(mainDig->GetEloss() + digToSum->GetEloss());
  
  Bool_t trackExists = kFALSE;
  for (Int_t iTrack=0; iTrack<mainDig->GetNMCTracks(); iTrack++) {
    if (digToSum->GetMCLabel(0) == mainDig->GetMCLabel(iTrack)) {
      trackExists = kTRUE;
      break;
    }
  }
  
  if (!trackExists) mainDig->AddMCLabel(digToSum->GetMCLabel(0));
  
}

//====================================================================================================================================================

