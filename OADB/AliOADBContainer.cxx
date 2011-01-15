/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     Offline Analysis Database Container and Service Class 
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------




#include "AliOADBContainer.h"
#include "AliLog.h"
#include <TObjArray.h>
#include <TArrayI.h>
#include <TFile.h>
#include <TList.h>

ClassImp(AliOADBContainer);

//______________________________________________________________________________
AliOADBContainer::AliOADBContainer() : 
  TNamed(),
  fArray(new TObjArray(100)),
  fDefaultList(new TList()),
  fLowerLimits(),
  fUpperLimits(),
  fEntries(0)
{
}

AliOADBContainer::AliOADBContainer(char* name) : 
  TNamed(name, "OADBContainer"),
  fArray(new TObjArray(100)),
  fDefaultList(new TList()),
  fLowerLimits(),
  fUpperLimits(),
  fEntries(0)
{
}


//______________________________________________________________________________
AliOADBContainer::~AliOADBContainer() 
{
  // destructor
}

//______________________________________________________________________________
AliOADBContainer::AliOADBContainer(const AliOADBContainer& cont) :
  TNamed(cont),
  fArray(cont.fArray),
  fDefaultList(cont.fDefaultList),
  fLowerLimits(cont.fLowerLimits),
  fUpperLimits(cont.fUpperLimits),
  fEntries(cont.fEntries)
{
  // Copy constructor.
}

//______________________________________________________________________________
AliOADBContainer& AliOADBContainer::operator=(const AliOADBContainer& cont)
{
  // Assignment operator
  if(this!=&cont) {
    TNamed::operator=(cont);
    fEntries = cont.fEntries;
    fLowerLimits.Set(fEntries);
    fUpperLimits.Set(fEntries);
    for (Int_t i = 0; i < fEntries; i++) {
	fLowerLimits[i] = cont.fLowerLimits[i]; 
	fUpperLimits[i] = cont.fUpperLimits[i];
	fArray->AddAt(cont.fArray->At(i), i);
    }
  }
  return *this;
}

void AliOADBContainer::AppendObject(TObject* obj, Int_t lower, Int_t upper)
{
  // Append a new object to the list 
  fEntries++;
  fLowerLimits.Set(fEntries);
  fUpperLimits.Set(fEntries);

  fLowerLimits[fEntries - 1] = lower;
  fUpperLimits[fEntries - 1] = upper;

  fArray->Add(obj);
}

void AliOADBContainer::RemoveObject(Int_t idx)
{
  // Remove object from the list 
  if (idx < 0 || idx >= fEntries) 
    {
      AliError(Form("Index out of Range %5d >= %5d", idx, fEntries));
      return;
    }
    
  TObject* obj = fArray->RemoveAt(idx);
  delete obj;

  for (Int_t i = idx; i < (fEntries-1); i++) {
    fLowerLimits[i] = fLowerLimits[i + 1]; 
    fUpperLimits[i] = fUpperLimits[i + 1];
    fArray->AddAt(fArray->At(i+1), i);
  }
  fArray->RemoveAt(fEntries - 1);
  fEntries--;
}

void AliOADBContainer::UpdateObject(Int_t idx, TObject* obj, Int_t lower, Int_t upper)
{
  // Append a new object to the list 
  if (idx < 0 || idx >= fEntries) 
    {
      AliError(Form("Index out of Range %5d >= %5d", idx, fEntries));
      return;
    }

  fLowerLimits[idx] = lower;
  fUpperLimits[idx] = upper;
  fArray->AddAt(obj, idx);
}

void  AliOADBContainer::AddDefaultObject(TNamed* obj)
{
  // Add a default object
  fDefaultList->Add(obj);
}

void  AliOADBContainer::CleanDefaultList()
{
  // Clean default list
  fDefaultList->Delete();
}

Int_t AliOADBContainer::GetIndexForRun(Int_t run) const
{
  //
  // Find the index for a given run 
  Int_t found = 0;
  Int_t index = -1;
  for (Int_t i = 0; i < fEntries; i++) 
    {
      if (run >= fLowerLimits[i] && run <= fUpperLimits[i])
	{
	  found++;
	  index = i;
	}
    }

  if (found > 1) {
    AliError(Form("More than one (%5d) object found; return last (%5d) !\n", found, index));
  } else if (index == -1) {
    AliWarning(Form("No object found for run %5d !\n", run));
  }
  
  return index;
}

TObject* AliOADBContainer::GetObject(Int_t run, char* def) const
{
  // Return object for given run or default if not found
  TObject* obj = 0;
  Int_t idx = GetIndexForRun(run);
  if (idx == -1) {
    // no object found, try default
    obj = fDefaultList->FindObject(def);
    if (!obj) {
      AliError("Default Object not found !\n");
      return (0);
    } else {
      return (obj);
    }
  } else {
    return (fArray->At(idx));
  }
}

TObject* AliOADBContainer::GetObjectByIndex(Int_t run) const
{
  // Return object for given index
  return (fArray->At(run));
}

void AliOADBContainer::WriteToFile(char* fname) const
{
  //
  // Write object to file
  TFile* f = new TFile(fname, "recreate");
  Write();
  f->Close();
}

Int_t AliOADBContainer::InitFromFile(char* fname, char* key)
{
    // 
    // Initialize object from file
    TFile* file = TFile::Open(fname);
    if (!file) return (1);
    AliOADBContainer* cont  = 0;
    file->GetObject(key, cont);
    if (!cont)
    {
	AliError("Object not found in file \n");	
	return 1;
    }
    
    fEntries = cont->GetNumberOfEntries();
    fLowerLimits.Set(fEntries);
    fUpperLimits.Set(fEntries);
    for (Int_t i = 0; i < fEntries; i++) {
	fLowerLimits[i] = cont->LowerLimit(i); 
	fUpperLimits[i] = cont->UpperLimit(i);
	fArray->AddAt(cont->GetObjectByIndex(i), i);
    }
    
    return 0;
    
}


void AliOADBContainer::List()
{
  //
  // List Objects
  for (Int_t i = 0; i < fEntries; i++) {
    printf("Lower %5d Upper %5d \n", fLowerLimits[i], fUpperLimits[i]);
    (fArray->At(i))->Dump();
  }
}
