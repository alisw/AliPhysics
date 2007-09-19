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
#include "AliLog.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDsegmentArray.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdigit.h"
#include "AliTRDgeometry.h"

#include "AliTRDSignalIndex.h"

ClassImp(AliTRDdigitsManager)

//_____________________________________________________________________________

  // Number of track dictionary arrays
  const Int_t AliTRDdigitsManager::fgkNDict = kNDict;

//_____________________________________________________________________________
AliTRDdigitsManager::AliTRDdigitsManager()
  :TObject()
  ,fEvent(0)
  ,fTree(0)
  ,fDigits(0)
  ,fIsRaw(0)
  ,fSDigits(0)
  ,fSignalIndexes(NULL)
  ,fUseDictionaries(kTRUE)
{
  //
  // Default constructor
  //

  for (Int_t iDict = 0; iDict < kNDict; iDict++) {
    fDictionary[iDict] = NULL;
  }
  
  //fSignalIndexes = new TList();
  fSignalIndexes = new TObjArray(AliTRDgeometry::Ndet());
  
}

//_____________________________________________________________________________
AliTRDdigitsManager::AliTRDdigitsManager(const AliTRDdigitsManager &m)
  :TObject(m)
  ,fEvent(m.fEvent)
  ,fTree(0)
  ,fDigits(0)
  ,fIsRaw(m.fIsRaw)
  ,fSDigits(m.fSDigits)
  ,fSignalIndexes(NULL)
  ,fUseDictionaries(kTRUE)
{
  //
  // AliTRDdigitsManager copy constructor
  //

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

  if (fSignalIndexes) {
    fSignalIndexes->Delete();
    delete fSignalIndexes;
  }
  fSignalIndexes = NULL;

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

//_____________________________________________________________________________
void AliTRDdigitsManager::Copy(TObject &m) const
{
  //
  // Copy function
  //

  ((AliTRDdigitsManager &) m).fIsRaw   = fIsRaw;
  ((AliTRDdigitsManager &) m).fEvent   = fEvent;
  ((AliTRDdigitsManager &) m).fSDigits = fSDigits;
  
  ((AliTRDdigitsManager &) m).fSignalIndexes = fSignalIndexes;
  ((AliTRDdigitsManager &) m).fUseDictionaries = fUseDictionaries;

  TObject::Copy(m);

}

//_____________________________________________________________________________
void AliTRDdigitsManager::CreateArrays()
{
  //
  // Create the data arrays
  //

  fDigits = new AliTRDsegmentArray("AliTRDdataArrayI",AliTRDgeometry::Ndet());

  if (fUseDictionaries)
    {
      for (Int_t iDict = 0; iDict < kNDict; iDict++) {
	fDictionary[iDict] = new AliTRDsegmentArray("AliTRDdataArrayI"
						    ,AliTRDgeometry::Ndet());
      }
    }

  for (Int_t i = 0; i < AliTRDgeometry::Ndet(); i++)
    {
      fSignalIndexes->AddLast(new AliTRDSignalIndex());
    }
}
//_____________________________________________________________________________
void AliTRDdigitsManager::ResetArrays()
{
  //
  // Reset the data arrays
  //

  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
  }
  fDigits = new AliTRDsegmentArray("AliTRDdataArrayI",AliTRDgeometry::Ndet());

  if (fUseDictionaries)
    {
      for (Int_t iDict = 0; iDict < kNDict; iDict++) {
	if (fDictionary[iDict]) { 
	  fDictionary[iDict]->Delete();
	  delete fDictionary[iDict];
	}
	fDictionary[iDict] = new AliTRDsegmentArray("AliTRDdataArrayI"
						    ,AliTRDgeometry::Ndet());
      }
    }

  for (Int_t i = 0; i < AliTRDgeometry::Ndet(); i++)
    {
      AliTRDSignalIndex *idx = (AliTRDSignalIndex *)fSignalIndexes->At(i);
      if (idx) idx->Reset();
    }
}

//_____________________________________________________________________________
void AliTRDdigitsManager::SetRaw()
{
  //
  // Switch on the raw digits flag
  //

  fIsRaw = kTRUE;
  if (fDigits)
    fDigits->SetBit(AliTRDdigit::RawDigit());
  
}

//_____________________________________________________________________________
Short_t AliTRDdigitsManager::GetDigitAmp(Int_t row, Int_t col,Int_t time
                                       , Int_t det) const
{
  //
  // Returns the amplitude of a digit
  //

  if (!GetDigits(det)) return 0;
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
      AliDebug(1,"Making branch TRDdigits\n");
    }
    else {
      status = kFALSE;
    }
  }
  else {
    status = kFALSE;
  }

  if (fUseDictionaries)
    {
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
	    AliDebug(1,Form("Making branch %s\n",branchname));
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
    AliDebug(1,"Create the data arrays.\n");
    CreateArrays();
  }

  status = fDigits->LoadArray("TRDdigits",fTree);

  if (fUseDictionaries)
    {
      for (Int_t iDict = 0; iDict < kNDict; iDict++) {
	Char_t branchname[15];
	sprintf(branchname,"TRDdictionary%d",iDict);
	status = fDictionary[iDict]->LoadArray(branchname,fTree);
	if (status == kFALSE)
	  {
	    fUseDictionaries = kFALSE;
	    AliWarning("Unable to load dict arrays. Will not use them.\n");
	    break;
	  }
      }  
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
    AliError("Error while storing digits in branch TRDdigits\n");
    return kFALSE;
  }

  if (fUseDictionaries)
    {
      for (Int_t iDict = 0; iDict < kNDict; iDict++) {
	Char_t branchname[15];
	sprintf(branchname,"TRDdictionary%d",iDict);
	if (!fDictionary[iDict]->StoreArray(branchname,fTree)) {
	  AliError(Form("Error while storing dictionary in branch %s\n",branchname));
	  return kFALSE;
	}
      }
    }

  // Write the new tree to the output file
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
    AliError(Form("track %d out of bounds (size: %d, this: 0x%08x)"
                 ,track,kNDict,this));
    return -1;
  }

  if (fUseDictionaries == kFALSE)
    {
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

  if (!fDigits) return 0x0;
  return (AliTRDdataArrayI *) fDigits->At(det);

}

//_____________________________________________________________________________
AliTRDdataArrayI *AliTRDdigitsManager::GetDictionary(Int_t det, Int_t i) const
{
  //
  // Returns the dictionary for one detector
  //
  if (fUseDictionaries == kFALSE)
    {
      return NULL;
    }

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
AliTRDSignalIndex *AliTRDdigitsManager::GetIndexes(Int_t det) 
{
  // 
  // Returns indexes of active pads
  //

  return (AliTRDSignalIndex*)fSignalIndexes->At(det);

}

//_____________________________________________________________________________
void AliTRDdigitsManager::RemoveDigits(Int_t det) 
{
  // 
  // Clear memory
  //

  fDigits->ClearSegment(det);

}

//_____________________________________________________________________________
void AliTRDdigitsManager::RemoveDictionaries(Int_t det) 
{
  // 
  // Clear memory
  //
  if (fUseDictionaries == kFALSE)
    {
      return;
    }

  for (Int_t i = 0; i < kNDict; i++) {
    fDictionary[i]->ClearSegment(det);
  }

}

//_____________________________________________________________________________
void AliTRDdigitsManager::ClearIndexes(Int_t det) 
{
  // 
  // Clear memory
  //
  fSignalIndexes->At(det)->Clear();  
}


//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::BuildIndexes(Int_t det)
{
  //
  // Build the list of indices
  //

  Int_t nRows = 0;
  Int_t nCols = 0;
  Int_t nTbins = 0;

  AliTRDgeometry    geom;
  AliTRDdataArrayI *digits = GetDigits(det);
  //digits should be expanded by now!!!
  if (digits->GetNtime() > 0) 
    {
      digits->Expand();
      nRows = digits->GetNrow();
      nCols = digits->GetNcol();
      nTbins = digits->GetNtime();

      //AliInfo(Form("rows %d cols %d tbins %d", nRows, nCols, nTbins));

      AliTRDSignalIndex* indexes = GetIndexes(det);
      indexes->SetSM(geom.GetSector(det));
      indexes->SetChamber(geom.GetChamber(det));
      indexes->SetPlane(geom.GetPlane(det));
      indexes->SetDetNumber(det);
      if (indexes->IsAllocated() == kFALSE)
	{
	  indexes->Allocate(nRows, nCols, nTbins);
	  //AliInfo(Form("Allocating 0x%x %d", indexes->GetArray(), indexes->GetArray()->GetSize()));
	}
      for (Int_t ir = 0; ir < nRows; ir++)
	{
	  for (Int_t ic = 0; ic < nCols; ic++)
	    {
	      for (Int_t it = 0; it < nTbins; it++)
		{	  
		  //AliInfo(Form("row %d col %d tbin %d", ir, ic, it));
		  
		  Int_t isig = digits->GetDataUnchecked(ir, ic, it);
  		  if (isig > 0)
		    {
		      //AliInfo(Form("row %d col %d tbin %d", ir, ic, it));
		      indexes->AddIndexTBin(ir, ic, it);	    
		    }
		} //tbins
	    } //cols
	} // rows
    } // if GetNtime
  else
    {
      return kFALSE;
    }
  return kTRUE;

}
