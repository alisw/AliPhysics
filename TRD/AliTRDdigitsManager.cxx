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

/*
$Log$
Revision 1.1.2.1  2000/05/08 14:44:01  cblume
Add new class AliTRDdigitsManager

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Manages the digits and the track dictionary in the form of               //
//  AliTRDdataArray objects.                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliRun.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDconst.h"

ClassImp(AliTRDdigitsManager)

//_____________________________________________________________________________
AliTRDdigitsManager::AliTRDdigitsManager():TObject()
{
  //
  // Default constructor
  //

  fIsRaw = kFALSE;

  fDigits = new AliTRDsegmentArray("AliTRDdataArrayI",kNdet);

  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    fDictionary[iDict] = new AliTRDsegmentArray("AliTRDdataArrayI",kNdet);
  }

}

//_____________________________________________________________________________
AliTRDdigitsManager::~AliTRDdigitsManager()
{

  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
  }

  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    fDictionary[iDict]->Delete();
    delete fDictionary[iDict];
  }

}

//_____________________________________________________________________________
void AliTRDdigitsManager::SetRaw()
{

  fIsRaw = kTRUE;

  fDigits->SetBit(kRawDigit);
  
}

//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::MakeBranch()
{
  //
  // Creates the branches for the digits and the dictionary in the digits tree
  //

  Int_t buffersize = 64000;

  Bool_t status = kTRUE;

  if (gAlice->TreeD()) {

    // Make the branch for the digits
    if (fDigits) {
      const AliTRDdataArrayI *Digits = 
           (AliTRDdataArrayI *) fDigits->At(0);
      if (Digits) {
        gAlice->TreeD()->Branch("TRDdigits",Digits->IsA()->GetName()
                                           ,&Digits,buffersize,1);
        printf("AliTRDdigitsManager::MakeBranch -- ");
        printf("Making branch TRDdigits\n");
      }
      else {
        status = kFALSE;
      }
    }
    else {
      status = kFALSE;
    }

    // Make the branches for the dictionaries
    for (Int_t iDict = 0; iDict < kNDict; iDict++) {

      Char_t branchname[15];
      sprintf(branchname,"TRDdictionary%d",iDict);
      if (fDictionary[iDict]) {
        const AliTRDdataArrayI *Dictionary = 
             (AliTRDdataArrayI *) fDictionary[iDict]->At(0);
        if (Dictionary) {
          gAlice->TreeD()->Branch(branchname,Dictionary->IsA()->GetName()
                                            ,&Dictionary,buffersize,1);
          printf("AliTRDdigitsManager::MakeBranch -- ");
          printf("Making branch %s\n",branchname);
	}
        else {
          status = kFALSE;
	}
      }
      else {
        status = kFALSE;
      }
    }

  }
  else {
    status = kFALSE;
  }

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::ReadDigits()
{

  Bool_t status = kTRUE;

  status = fDigits->LoadArray("TRDdigits");

  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    Char_t branchname[15];
    sprintf(branchname,"TRDdictionary%d",iDict);
    status = fDictionary[iDict]->LoadArray(branchname);
  }  

  if (fDigits->TestBit(kRawDigit)) {
    fIsRaw = kTRUE;
  }
  else {
    fIsRaw = kFALSE;
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::WriteDigits()
{
  //
  // Writes out the TRD-digits and the dictionaries
  //

  // Create the branches
  if (!(gAlice->TreeD()->GetBranch("TRDdigits"))) { 
    if (!MakeBranch()) return kFALSE;
  }

  // Store the contents of the segment array in the tree
  if (!fDigits->StoreArray("TRDdigits")) {
    printf("AliTRDdigitsManager::WriteDigits -- ");
    printf("Error while storing digits in branch TRDdigits\n");
    return kFALSE;
  }
  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    Char_t branchname[15];
    sprintf(branchname,"TRDdictionary%d",iDict);
    if (!fDictionary[iDict]->StoreArray(branchname)) {
      printf("AliTRDdigitsManager::WriteDigits -- ");
      printf("Error while storing dictionary in branch %s\n",branchname);
      return kFALSE;
    }
  }

  return kTRUE;

}
