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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Manages the digits and the track dictionary in the form of               //
//  AliTRDdataArray objects.                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
 
#include <TROOT.h>
#include <TTree.h>                                                              
#include <TFile.h>

#include "AliRun.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDsegmentArray.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdigit.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDdigitsManager)

//_____________________________________________________________________________

  // Number of track dictionary arrays
  const Int_t AliTRDdigitsManager::fgkNDict = kNDict;

//_____________________________________________________________________________
AliTRDdigitsManager::AliTRDdigitsManager():TObject()
{
  //
  // Default constructor
  //

  fIsRaw   = kFALSE;
  fEvent   = 0;
  fDebug   = 0;
  fSDigits = 0;

  fTree    = NULL;
  fDigits  = NULL;
  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    fDictionary[iDict] = NULL;
  }

}

//_____________________________________________________________________________
AliTRDdigitsManager::AliTRDdigitsManager(const AliTRDdigitsManager &m)
:TObject(m)
{
  //
  // AliTRDdigitsManager copy constructor
  //

  ((AliTRDdigitsManager &) m).Copy(*this);

}

//_____________________________________________________________________________
AliTRDdigitsManager::~AliTRDdigitsManager()
{
  //
  // AliTRDdigitsManager destructor
  //

  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    fDigits = NULL;
  }

  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    fDictionary[iDict]->Delete();
    delete fDictionary[iDict];
    fDictionary[iDict] = NULL;
  }

}

//_____________________________________________________________________________
void AliTRDdigitsManager::Copy(TObject &m) const
{
  //
  // Copy function
  //

  ((AliTRDdigitsManager &) m).fIsRaw   = fIsRaw;
  ((AliTRDdigitsManager &) m).fEvent   = fEvent;
  ((AliTRDdigitsManager &) m).fDebug   = fDebug;
  ((AliTRDdigitsManager &) m).fSDigits = fSDigits;

  TObject::Copy(m);

}

//_____________________________________________________________________________
void AliTRDdigitsManager::CreateArrays()
{
  //
  // Create the data arrays
  //

  fDigits = new AliTRDsegmentArray("AliTRDdataArrayI",AliTRDgeometry::Ndet());

  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    fDictionary[iDict] = new AliTRDsegmentArray("AliTRDdataArrayI"
                                               ,AliTRDgeometry::Ndet());
  }

}
//_____________________________________________________________________________
void AliTRDdigitsManager::ResetArrays()
{
  //
  // Reset the data arrays
  //

  if (fDigits) {
    delete fDigits;
  }
  fDigits = new AliTRDsegmentArray("AliTRDdataArrayI",AliTRDgeometry::Ndet());

  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    if (fDictionary[iDict]) {  
      delete fDictionary[iDict];
    }
    fDictionary[iDict] = new AliTRDsegmentArray("AliTRDdataArrayI"
                                               ,AliTRDgeometry::Ndet());
  }

}

//_____________________________________________________________________________
void AliTRDdigitsManager::SetRaw()
{
  //
  // Switch on the raw digits flag
  //

  fIsRaw = kTRUE;

  fDigits->SetBit(AliTRDdigit::RawDigit());
  
}

//_____________________________________________________________________________
Short_t AliTRDdigitsManager::GetDigitAmp(Int_t row, Int_t col,Int_t time
                                       , Int_t det) const
{
  //
  // Returns the amplitude of a digit
  //

  return ((Short_t) GetDigits(det)->GetData(row,col,time));

}
 
//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::MakeBranch(TTree *tree)
{
  //
  // Creates the tree and branches for the digits and the dictionary
  //

  Int_t buffersize = 64000;

  Bool_t status = kTRUE;

  if (tree) {
    fTree = tree;
  }

  // Make the branch for the digits
  if (fDigits) {
    const AliTRDdataArray *kDigits = (AliTRDdataArray *) fDigits->At(0);
    if (kDigits) {
      if (!fTree) return kFALSE;
      TBranch* branch = fTree->GetBranch("TRDdigits");
      if (!branch) fTree->Branch("TRDdigits",kDigits->IsA()->GetName(),
                                 &kDigits,buffersize,99);
      if (fDebug > 0) {
        printf("<AliTRDdigitsManager::MakeBranch> ");
        printf("Making branch TRDdigits\n");
      }
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
      const AliTRDdataArray *kDictionary = 
              (AliTRDdataArray *) fDictionary[iDict]->At(0);
      if (kDictionary) {
	if (!fTree) return kFALSE;
	TBranch* branch = fTree->GetBranch(branchname);
	if (!branch) fTree->Branch(branchname,kDictionary->IsA()->GetName(),
				   &kDictionary,buffersize,99);
        if (fDebug > 0) {
          printf("<AliTRDdigitsManager::MakeBranch> ");
          printf("Making branch %s\n",branchname);
	}
      }
      else {
        status = kFALSE;
      }
    }
    else {
      status = kFALSE;
    }
  }

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::ReadDigits(TTree *tree)
{
  //
  // Reads the digit information from the input file
  //

  Bool_t status = kTRUE;

  if (tree) {

    fTree = tree;

  }

  if (!fDigits) {
    if (fDebug > 0) {
      printf("<AliTRDdigitsManager::ReadDigits> ");
      printf("Create the data arrays.\n");
    }
    CreateArrays();
  }

  status = fDigits->LoadArray("TRDdigits",fTree);

  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    Char_t branchname[15];
    sprintf(branchname,"TRDdictionary%d",iDict);
    status = fDictionary[iDict]->LoadArray(branchname,fTree);
  }  

  if (fDigits->TestBit(AliTRDdigit::RawDigit())) {
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

  // Store the contents of the segment array in the tree
  if (!fDigits->StoreArray("TRDdigits",fTree)) {
    printf("<AliTRDdigitsManager::WriteDigits> ");
    printf("Error while storing digits in branch TRDdigits\n");
    return kFALSE;
  }
  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    Char_t branchname[15];
    sprintf(branchname,"TRDdictionary%d",iDict);
    if (!fDictionary[iDict]->StoreArray(branchname,fTree)) {
      printf("<AliTRDdigitsManager::WriteDigits> ");
      printf("Error while storing dictionary in branch %s\n",branchname);
      return kFALSE;
    }
  }

  // Write the new tree to the output file
  //fTree->Write();
  fTree->AutoSave();  // Modification by Jiri

  return kTRUE;

}

//_____________________________________________________________________________
AliTRDdigit *AliTRDdigitsManager::GetDigit(Int_t row, Int_t col
                                         , Int_t time, Int_t det) const
{
  // 
  // Creates a single digit object 
  //

  Int_t digits[4];
  Int_t amp[1];

  digits[0] = det;
  digits[1] = row;
  digits[2] = col;
  digits[3] = time;

  amp[0]    = GetDigits(det)->GetData(row,col,time);
  
  return (new AliTRDdigit(fIsRaw,digits,amp));

}

//_____________________________________________________________________________
Int_t AliTRDdigitsManager::GetTrack(Int_t track
                                  , Int_t row, Int_t col, Int_t time
                                  , Int_t det) const
{
  // 
  // Returns the MC-track numbers from the dictionary.
  //

  if ((track < 0) || (track >= kNDict)) {
    TObject::Error("GetTracks"
                  ,"track %d out of bounds (size: %d, this: 0x%08x)"
                  ,track,kNDict,this);
    return -1;
  }

  // Array contains index+1 to allow data compression
  return (GetDictionary(det,track)->GetData(row,col,time) - 1);

}

//_____________________________________________________________________________
AliTRDdataArrayI *AliTRDdigitsManager::GetDigits(Int_t det) const
{
  //
  // Returns the digits array for one detector
  //

  return (AliTRDdataArrayI *) fDigits->At(det);

}

//_____________________________________________________________________________
AliTRDdataArrayI *AliTRDdigitsManager::GetDictionary(Int_t det, Int_t i) const
{
  //
  // Returns the dictionary for one detector
  //

  return (AliTRDdataArrayI *) fDictionary[i]->At(det);

}

//_____________________________________________________________________________
Int_t AliTRDdigitsManager::GetTrack(Int_t track, AliTRDdigit *Digit) const
{
  // 
  // Returns the MC-track numbers from the dictionary for a given digit
  //

  Int_t row  = Digit->GetRow();
  Int_t col  = Digit->GetCol();
  Int_t time = Digit->GetTime();
  Int_t det  = Digit->GetDetector();

  return GetTrack(track,row,col,time,det);

}

//_____________________________________________________________________________
AliTRDdigitsManager &AliTRDdigitsManager::operator=(const AliTRDdigitsManager &m)
{
  //
  // Assignment operator
  //

  if (this != &m) ((AliTRDdigitsManager &) m).Copy(*this);
  return *this;

}
