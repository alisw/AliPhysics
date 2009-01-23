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
//  TObjArray objects                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h> 
#include <TROOT.h>
#include <TTree.h>                                                              
#include <TFile.h>
#include "AliRun.h"
#include "AliLog.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDarrayDictionary.h"
#include "AliTRDarrayADC.h"
#include "AliTRDarraySignal.h"
#include "AliTRDdigit.h"
#include "AliTRDgeometry.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDReconstructor.h"

ClassImp(AliTRDdigitsManager)

//_____________________________________________________________________________

  // Number of track dictionary arrays
  const Int_t AliTRDdigitsManager::fgkNDict = kNDict;

//_____________________________________________________________________________
AliTRDdigitsManager::AliTRDdigitsManager(const AliTRDReconstructor *const rec)
  :TObject()
  ,fEvent(0)
  ,fTree(0)
  ,fDigits(0) 
  ,fHasSDigits(0)
  ,fSignalIndexes(NULL)
  ,fUseDictionaries(kTRUE)
  ,fTreeD(0)
  ,fBranch(0)
  ,fDets(AliTRDgeometry::Ndet())
  ,fRawRec(kFALSE)
{
  //
  // Default constructor
  //
  
  if(rec)
    {
      fDets=1;
      fRawRec=kTRUE;
    }

  for (Int_t iDict = 0; iDict < kNDict; iDict++) 
    {
      fDict[iDict] = NULL;
    }  
}

//_____________________________________________________________________________
AliTRDdigitsManager::AliTRDdigitsManager(const AliTRDdigitsManager &m)
  :TObject(m)
  ,fEvent(m.fEvent)
  ,fTree(0)
  ,fDigits(0) 
  ,fHasSDigits(m.fHasSDigits)
  ,fSignalIndexes(NULL)
  ,fUseDictionaries(kTRUE)
  ,fTreeD(m.fTree)
  ,fBranch(m.fBranch)
  ,fDets(m.fDets)
  ,fRawRec(m.fRawRec)
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


  if (fDigits) 
    {
      fDigits->Delete();
      delete fDigits;
      fDigits = NULL;
    }

  for (Int_t iDict = 0; iDict < kNDict; iDict++) 
    {
      if(fDict[iDict])
	{
	  fDict[iDict]->Delete();
	  delete fDict[iDict];
	  fDict[iDict] = NULL;
	}
    }

  if (fSignalIndexes) 
    {
      fSignalIndexes->Delete();
      delete fSignalIndexes;
      fSignalIndexes = NULL;
    }

}

//_____________________________________________________________________________
AliTRDdigitsManager &AliTRDdigitsManager::operator=(const AliTRDdigitsManager &m)
{
  //
  // Assignment operator
  //

  if (this != &m) 
    {
      ((AliTRDdigitsManager &) m).Copy(*this);
    }

  return *this;

}

//_____________________________________________________________________________
void AliTRDdigitsManager::Copy(TObject &m) const
{
  //
  // Copy function
  //

  ((AliTRDdigitsManager &) m).fEvent           = fEvent;
  ((AliTRDdigitsManager &) m).fHasSDigits      = fHasSDigits;
  ((AliTRDdigitsManager &) m).fDigits          = fDigits;
  for(Int_t i=0; i<kNDict; i++)
    {
      ((AliTRDdigitsManager &) m).fDict[i]  = fDict[i];
    }
  ((AliTRDdigitsManager &) m).fSignalIndexes   = fSignalIndexes;
  ((AliTRDdigitsManager &) m).fUseDictionaries = fUseDictionaries;
  ((AliTRDdigitsManager &) m).fDets            = fDets;
  ((AliTRDdigitsManager &) m).fRawRec           = fRawRec;

  TObject::Copy(m);

}

//_____________________________________________________________________________
void AliTRDdigitsManager::CreateArrays()
{
  //
  // Create the data arrays
  //

  if (fHasSDigits) 
    {
      if(fDigits)                                        
	{                                                   
	  fDigits->Delete();                                
	  delete fDigits;                                   
	}                                                    
      fDigits = new TObjArray(fDets);
      for (Int_t index = 0; index < fDets; index++) 
	fDigits->AddAt(new AliTRDarraySignal(),index);
    }
  else 
    {
      if(fDigits)                                          
	{                                                    
	  fDigits->Delete();                                
	  delete fDigits;                                   
	}                                                   
      fDigits = new TObjArray(fDets);    
      for (Int_t index = 0; index < fDets;index++) 
	fDigits->AddAt(new AliTRDarrayADC(),index);
    }

  if (fUseDictionaries) 
    {
      for(Int_t iDict = 0; iDict < kNDict; iDict++)
	if(fDict[iDict])                                           
	  {
	    fDict[iDict]->Delete();                                
	    delete fDict[iDict];                                    
	  }
      for(Int_t iDict = 0; iDict < kNDict; iDict++)
	fDict[iDict] = new TObjArray(fDets);

      for (Int_t iDict = 0; iDict < kNDict; iDict++)
	for (Int_t index = 0; index < fDets; index++) 
	  fDict[iDict]->AddAt(new AliTRDarrayDictionary(),index);
    }
  
  if(fSignalIndexes)
    {
      fSignalIndexes->Delete();
      delete fSignalIndexes;
    }
  fSignalIndexes = new TObjArray(fDets);
  for (Int_t i = 0; i < fDets; i++) 
    fSignalIndexes->AddLast(new AliTRDSignalIndex());

}

//_____________________________________________________________________________
void AliTRDdigitsManager::ResetArrays()
{
  //
  // Reset the data arrays
  //

  if (fDigits)
    {
      fDigits->Delete();
      delete fDigits;
    }
  if (fHasSDigits)
    {
      fDigits = new TObjArray(fDets);     
      for (Int_t index = 0; index < fDets; index++) 
	fDigits->AddAt(new AliTRDarraySignal(),index);
    }
  else
    {
      fDigits = new TObjArray(fDets);
      for (Int_t index = 0; index < fDets; index++)
	fDigits->AddAt(new AliTRDarrayADC(),index);
    }
  
  for (Int_t iDict = 0; iDict < kNDict; iDict++)
    {
      if (fDict[iDict])
	{
	  fDict[iDict]->Delete();
	  delete fDict[iDict];
	  fDict[iDict]=NULL;
	}
    }
  if (fUseDictionaries) 
    {
      for(Int_t iDict = 0; iDict < kNDict; iDict++)
	fDict[iDict] = new TObjArray(fDets);
      
      for (Int_t iDict = 0; iDict < kNDict; iDict++)
	for (Int_t index = 0; index < fDets; index++)
	  fDict[iDict]->AddAt(new AliTRDarrayDictionary(),index);
    }
  
  if(fSignalIndexes)
    {
      fSignalIndexes->Delete();
      delete fSignalIndexes;
    }
  fSignalIndexes = new TObjArray(fDets);
  for (Int_t i = 0; i < fDets; i++)
    fSignalIndexes->AddLast(new AliTRDSignalIndex());
}

//_____________________________________________________________________________
void AliTRDdigitsManager::ResetArrays(Int_t det)
{
  Int_t RecoDet = fRawRec ? 0 : det;

  RemoveDigits(RecoDet);
  RemoveDictionaries(RecoDet);
  RemoveIndexes(RecoDet);

  if (fHasSDigits)
    fDigits->AddAt(new AliTRDarraySignal(),RecoDet);
  else
    fDigits->AddAt(new AliTRDarrayADC(),RecoDet);

  if (fUseDictionaries) 
    {
      for (Int_t iDict = 0; iDict < kNDict; iDict++)
	fDict[iDict]->AddAt(new AliTRDarrayDictionary(),RecoDet);
    }
  
  fSignalIndexes->AddAt(new AliTRDSignalIndex(),RecoDet);
}

//_____________________________________________________________________________
Short_t AliTRDdigitsManager::GetDigitAmp(Int_t row, Int_t col,Int_t time, Int_t det) const
{
  //
  // Returns the amplitude of a digit
  //

  if (!GetDigits(det)) return 0;
  
  return ((Short_t) ((AliTRDarrayADC *) GetDigits(det))->GetDataB(row,col,time));

}

//_____________________________________________________________________________
UChar_t AliTRDdigitsManager::GetPadStatus(Int_t row, Int_t col, Int_t time, Int_t det) const
{
  //
  // Returns the pad status for the requested pad
  //
	
  if (!GetDigits(det)) return 0;

  return ((UChar_t) ((AliTRDarrayADC *) GetDigits(det))->GetPadStatus(row,col,time));
 
}

//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::MakeBranch(TTree *tree)  
{
  //
  // Creates the tree and branches for the digits and the dictionary
  //

  Int_t  buffersize = 64000;
  Bool_t status     = kTRUE;

  if (tree) 
    {
      fTree = tree;
    }

  // Make the branch for the digits
  if (fDigits) 
    {
      if(fHasSDigits)
	{
	  const AliTRDarraySignal *kDigits = (AliTRDarraySignal *) fDigits->At(0); 
	  if (kDigits) 
	    {
	      if (!fTree) return kFALSE;
	      AliDebug(1,"Making branch for SDigits!\n");
	      TBranch* branch = fTree->GetBranch("TRDdigits");
	      if (!branch) fTree->Branch("TRDdigits","AliTRDarraySignal",&kDigits,buffersize,99);
	      AliDebug(1,"Making branch TRDdigits\n");
	    }
	  else 
	    {
	      status = kFALSE;
	    }
	}

      if(!fHasSDigits)
	{
	  const AliTRDarrayADC *kDigits = (AliTRDarrayADC *) fDigits->At(0);
	  if (kDigits) 
	    {
	      if (!fTree) return kFALSE;
	      AliDebug(1,"Making branch for Digits!\n");
	      TBranch* branch = fTree->GetBranch("TRDdigits");
	      if (!branch) fTree->Branch("TRDdigits","AliTRDarrayADC",&kDigits,buffersize,99);
	      AliDebug(1,"Making branch TRDdigits\n");	      
	    }
	  else 
	    {
	      status = kFALSE;
	    }
	}

    }    
  else
    {
      status = kFALSE;
    }
  
  if (fUseDictionaries) 
    {
      // Make the branches for the dictionaries
      for (Int_t iDict = 0; iDict < kNDict; iDict++) 
	{
	  Char_t branchname[15];
	  sprintf(branchname,"TRDdictionary%d",iDict); 
	  if (fDict[iDict]) 
	    {
	      const AliTRDarrayDictionary *kDictionary = (AliTRDarrayDictionary *) fDict[iDict]->At(0);
	      if (kDictionary) 
		{
		  if (!fTree) return kFALSE;
		  AliDebug(2,"Making branch for dictionary!\n");
		  TBranch* branch = fTree->GetBranch(branchname);
		  if (!branch) fTree->Branch(branchname,"AliTRDarrayDictionary",&kDictionary,buffersize,99);
		  AliDebug(1,Form("Making branch %s\n",branchname));
		}
	      else 
		{
		  status = kFALSE;
		}
	    }
	  else 
	    {
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

  if (tree) 
    {
      fTree = tree;
    }

  if (!fDigits) 
    {
      AliDebug(1,"Create the data arrays.\n");
      CreateArrays();
    }

  status = LoadArray(fDigits,"TRDdigits",fTree);

  if (fUseDictionaries) 
    {
      for (Int_t iDict = 0; iDict < kNDict; iDict++) 
	{
	  Char_t branchname[15];
	  sprintf(branchname,"TRDdictionary%d",iDict);
	  status = LoadArrayDict(fDict[iDict],branchname,fTree);
	  if (status == kFALSE) 
	    {
	      fUseDictionaries = kFALSE;
	      AliWarning("Unable to load dict arrays. Will not use them.\n");
	      break;
	    }
	}  
    }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::WriteDigits()
{
  //
  // Writes out the TRD-digits and the dictionaries
  //

  // Store the contents of the detector array in the tree

  if (!StoreArray(fDigits,"TRDdigits",fTree))
    {
      AliError("Error while storing digits in branch TRDdigits\n");
      return kFALSE;
    }

  if (fUseDictionaries) 
    {
      for (Int_t iDict = 0; iDict < kNDict; iDict++)
	{
	  Char_t branchname[15];
	  sprintf(branchname,"TRDdictionary%d",iDict);
	  if (!StoreArrayDict(fDict[iDict],branchname,fTree)) 
	    {
	      AliError(Form("Error while storing dictionary in branch %s\n",branchname));
	      return kFALSE;
	    }
	}
    }
  
  // Write the new tree to the output file
  fTree->AutoSave();

  return kTRUE;

}

//_____________________________________________________________________________
AliTRDdigit *AliTRDdigitsManager::GetDigit(Int_t row
                                         , Int_t col
                                         , Int_t time
                                         , Int_t det) const
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

  amp[0]    = ((AliTRDarrayADC *) GetDigits(det))->GetData(row,col,time);
  
  return (new AliTRDdigit(digits,amp));

}

//_____________________________________________________________________________
Int_t AliTRDdigitsManager::GetTrack(Int_t track
                                  , Int_t row
                                  , Int_t col
                                  , Int_t time
                                  , Int_t det) const
{
  // 
  // Returns the MC-track numbers from the dictionary.
  //

  if ((track < 0) || (track >= kNDict)) 
    {
      AliError(Form("track %d out of bounds (size: %d, this: 0x%08x)",track,kNDict,this));
      return -1;
    }

  if (fUseDictionaries == kFALSE) 
    {
      return -1;
    }

  // Array contains index+1 to allow data compression--->Changed
  return (((AliTRDarrayDictionary *) GetDictionary(det,track))->GetData(row,col,time) );

}

//________________________________________________________________________________
AliTRDarrayADC *AliTRDdigitsManager::GetDigits(Int_t det) const
{
  //
  // Returns the digits array for one detector
  //

  Int_t RecoDet = fRawRec ? 0 : det;

  if (!fDigits)   
    {
      return 0x0;
    }

  if (!fHasSDigits)
    {
      ((AliTRDarrayADC *) fDigits->At(RecoDet))->SetNdet(det);
      return (AliTRDarrayADC *) fDigits->At(RecoDet); 
    }
  else
    {
      AliDebug(2,"ERROR IN DATA TYPE!!!!");
      return 0x0;
    }

}

//_____________________________________________________________________________
AliTRDarraySignal *AliTRDdigitsManager::GetSDigits(Int_t det) const
{
  //
  // Returns the sdigits array for one detector
  //

  Int_t RecoDet = fRawRec ? 0 : det;

  if (!fDigits)   
    {
      //      AliDebug(1,"NO FDIGITS!");	
      return 0x0;
    }

  if (fHasSDigits)
    {
      ((AliTRDarraySignal *) fDigits->At(RecoDet))->SetNdet(det);
      return (AliTRDarraySignal *) fDigits->At(RecoDet);
    }
  else
    {
      AliDebug(2,"ERROR IN DATA TYPE!!!!");
      return 0x0;
    }

}

//_____________________________________________________________________________
AliTRDarrayDictionary *AliTRDdigitsManager::GetDictionary(Int_t det
                                                        , Int_t i) const
{
  //
  // Returns the dictionary for one detector
  //

  Int_t RecoDet = fRawRec ? 0 : det;

  if (fUseDictionaries == kFALSE)
    {
      return 0x0;
    }

  ((AliTRDarrayDictionary *) fDigits->At(RecoDet))->SetNdet(det);
  return (AliTRDarrayDictionary *) fDict[i]->At(RecoDet);
  
}

//_____________________________________________________________________________
Int_t AliTRDdigitsManager::GetTrack(Int_t track, AliTRDdigit *digit) const
{
  // 
  // Returns the MC-track numbers from the dictionary for a given digit
  //

  Int_t row  = digit->GetRow();
  Int_t col  = digit->GetCol();
  Int_t time = digit->GetTime();
  Int_t det  = digit->GetDetector();

  return GetTrack(track,row,col,time,det);

}

//_____________________________________________________________________________
AliTRDSignalIndex *AliTRDdigitsManager::GetIndexes(Int_t det) 
{
  // 
  // Returns indexes of active pads
  //

  Int_t RecoDet = fRawRec ? 0 : det;

  return (AliTRDSignalIndex *) fSignalIndexes->At(RecoDet);

}

//_____________________________________________________________________________
void AliTRDdigitsManager::RemoveDigits(Int_t det) 
{
   // 
   // Clear memory at det for Digits
   //

  Int_t RecoDet = fRawRec ? 0 : det;

  if (fDigits->At(RecoDet))
    {
      if (fHasSDigits) 
        {
          AliTRDarraySignal *arr = (AliTRDarraySignal *) fDigits->RemoveAt(RecoDet);
          delete arr;
	}
      else 
        {
          AliTRDarrayADC    *arr = (AliTRDarrayADC *)    fDigits->RemoveAt(RecoDet);
          delete arr;
	}
    }

}

//_____________________________________________________________________________
void AliTRDdigitsManager::RemoveDictionaries(Int_t det) 
{
  // 
  // Clear memory
  //

  Int_t RecoDet = fRawRec ? 0 : det;

  if (fUseDictionaries == kFALSE) 
    {
      return;
    }

  for (Int_t i = 0; i < kNDict; i++) 
    {
      if (fDict[i]->At(RecoDet))
	{
	  AliTRDarrayDictionary *arr = (AliTRDarrayDictionary *) fDict[i]->RemoveAt(RecoDet);
          delete arr;
	}
    }

}

//_____________________________________________________________________________
void AliTRDdigitsManager::RemoveIndexes(Int_t det) 
{
   // 
   // Clear memory
   //

  Int_t RecoDet = fRawRec ? 0 : det;

  if (fSignalIndexes->At(RecoDet))
    {
      AliTRDSignalIndex *arr = (AliTRDSignalIndex *) fSignalIndexes->RemoveAt(RecoDet);
      delete arr;
    }

}


//_____________________________________________________________________________
void AliTRDdigitsManager::ClearIndexes(Int_t det) 
{
  // 
  // Clear memory
  //
  
  Int_t RecoDet = fRawRec ? 0 : det;

  ((AliTRDSignalIndex *) fSignalIndexes->At(RecoDet))->ClearAll();  

}

//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::BuildIndexes(Int_t det)
{
  //
  // Build the list of indices
  //

  Int_t nRows  = 0;
  Int_t nCols  = 0;
  Int_t nTbins = 0;

  AliTRDgeometry  geom;
  AliTRDarrayADC *digits = 0x0;

  if (fHasSDigits) 
    {
      return kFALSE;
    }
  else 
    {
      digits = (AliTRDarrayADC *) GetDigits(det);
    }

  //digits should be expanded by now!!!
  if (digits->GetNtime() > 0) 
    {      
      digits->Expand(); 
      nRows  = digits->GetNrow();
      nCols  = digits->GetNcol();
      nTbins = digits->GetNtime();
      
      AliTRDSignalIndex *indexes = GetIndexes(det);
      indexes->SetSM(geom.GetSector(det));
      indexes->SetStack(geom.GetStack(det));
      indexes->SetLayer(geom.GetLayer(det));
      indexes->SetDetNumber(det);

      if (indexes->IsAllocated() == kFALSE)
	{
	  indexes->Allocate(nRows,nCols,nTbins);
	}

      for (Int_t ir = 0; ir < nRows; ir++) 
	{
	  for (Int_t ic = 0; ic < nCols; ic++) 
	    {
	      for (Int_t it = 0; it < nTbins; it++)
		{	  
		  Int_t isig = digits->GetDataB(ir,ic,it);
		  if (isig > 0) 
		    {
		      indexes->AddIndexRC(ir,ic);	    
		    }
		} // tbins
	    } // cols
	} // rows

    } // if GetNtime
  else 
    {
      return kFALSE;
    }
  
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::LoadArray(TObjArray *object
                                    , const Char_t *branchname
                                    , TTree *tree)
{
  //
  // Loads all detectors of the array from the branch <branchname> of
  // the digits tree <tree>
  // Adapted from code of the class AliTRDsegmentArray
  //

  fTreeD = tree;

  if (!fTreeD) 
    {
      AliError("Digits tree is not defined\n");
      return kFALSE;
    }

  // Get the branch
  fBranch = fTreeD->GetBranch(branchname);
  if (!fBranch) 
    {
      AliError(Form("Branch %s is not defined\n",branchname));
      return kFALSE;
    }

  // Loop through all detectors and read them from the tree
  Bool_t status = kTRUE;
  for (Int_t iDet = 0; iDet < fDets; iDet++) 
    {
      if(fHasSDigits)
	{
	  AliTRDarraySignal *dataArray = (AliTRDarraySignal *) object->At(iDet);
	  if (!dataArray) 
	    {
	      status = kFALSE;
	      break;    
	    }

	  fBranch->SetAddress(&dataArray);
	  fBranch->GetEntry(iDet);
	}
      else
	{
	  AliTRDarrayADC *dataArray = (AliTRDarrayADC *) object->At(iDet);
	  if (!dataArray) 
	    {
	      status = kFALSE;
	      break;    
	    }
	  fBranch->SetAddress(&dataArray);
	  fBranch->GetEntry(iDet);
	}
    }

  return status;

}

//________________________________________________________________________________________________
Bool_t AliTRDdigitsManager::LoadArrayDict(TObjArray *object
                                        , const Char_t *branchname
                                        , TTree *tree)
{
  //
  // Loads all detectors of the array from the branch <branchname> of
  // the dictionary tree <tree>
  // Adapted from code of the class AliTRDsegmentArray
  //

  fTreeD = tree;

  if (!fTreeD) 
    {
      AliError("Digits tree is not defined\n");
      return kFALSE;
    }

  // Get the branch
  fBranch = fTreeD->GetBranch(branchname);
  if (!fBranch) 
    {
      AliError(Form("Branch %s is not defined\n",branchname));
      return kFALSE;
    }

  // Loop through all detectors and read them from the tree
  Bool_t status = kTRUE;
  for (Int_t iDet = 0; iDet < fDets; iDet++) 
    {
      AliTRDarrayDictionary *dataArray = (AliTRDarrayDictionary *) object->At(iDet);
      if (!dataArray) 
	{
	  status = kFALSE;
	  break;    
	}
      fBranch->SetAddress(&dataArray);
      fBranch->GetEntry(iDet);
    }

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::StoreArray(TObjArray *array1
                                     , const Char_t *branchname
                                     , TTree *tree)
{
  //
  // Stores all the detectors of the array in the branch <branchname> of 
  // the digits tree <tree>
  // Adapted from code of the class AliTRDsegmentArray
  //

  fTree = tree;

  if (!fTree) 
    {
      AliError("Digits tree is not defined\n");
      return kFALSE;
    }

  // Get the branch
  fBranch = fTree->GetBranch(branchname);
  if (!fBranch) 
    {
      AliError(Form("Branch %s is not defined\n",branchname));
      return kFALSE;
    }

  // Loop through all detectors and fill them into the tree
  Bool_t status = kTRUE;
  for (Int_t iDet = 0; iDet < fDets; iDet++) 
    {
      if (fHasSDigits)
	{
	  const AliTRDarraySignal *kDataArray = (AliTRDarraySignal *) array1->At(iDet);
	  if (!kDataArray) 
	    {
	      status = kFALSE;
	      break;
	    }
	  fBranch->SetAddress(&kDataArray);
	  fBranch->Fill();
	}
      else
	{
	  const AliTRDarrayADC *kDataArray = (AliTRDarrayADC *) array1->At(iDet); 
	  if (!kDataArray) 
	    {
	      status = kFALSE;
	      break;
	    }
	  fBranch->SetAddress(&kDataArray);
	  fBranch->Fill();
	}
    }

  return status;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitsManager::StoreArrayDict(TObjArray *array3
                                         , const Char_t *branchname
                                         , TTree *tree)
{
  //
  // Stores all the dictionary arrays of the detectors of the array in the branch <branchname> of 
  // the dictionary tree <tree>
  // Adapted from code of the class AliTRDsegmentArray
  //

  //  AliDebug(1,"Storing Arrays of Dictionary");
  fTree = tree;

  if (!fTree) 
    {
      AliError("Digits tree is not defined\n");
      return kFALSE;
    }

  // Get the branch
  fBranch = fTree->GetBranch(branchname);
  if (!fBranch) 
    {
      AliError(Form("Branch %s is not defined\n",branchname));
      return kFALSE;
    }

  // Loop through all detectors and fill them into the tree
  Bool_t status = kTRUE;
  for (Int_t iDet = 0; iDet < fDets; iDet++) 
    {
      const AliTRDarrayDictionary *kDataArray = (AliTRDarrayDictionary *) array3->At(iDet);
      if (!kDataArray) 
	{
	  status = kFALSE;
	  break;
	}
      fBranch->SetAddress(&kDataArray);
      fBranch->Fill();
    }

  return status;

}
