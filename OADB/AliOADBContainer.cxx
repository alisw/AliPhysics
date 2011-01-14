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

ClassImp(AliOADBContainer);

//______________________________________________________________________________
AliOADBContainer::AliOADBContainer() : 
  TNamed(),
  fArray(new TObjArray(100)),
  fLowerLimits(),
  fUpperLimits(),
  fEntries(0)
{
}

AliOADBContainer::AliOADBContainer(char* name) : 
  TNamed(name, "OADBContainer"),
  fArray(new TObjArray(100)),
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

Int_t AliOADBContainer::GetIndexForRun(Int_t run) 
{
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
  return index;
}

void AliOADBContainer::WriteToFile(char* fname)
{
  // Write object to file
  TFile* f = new TFile(fname, "recreate");
  Write();
  f->Close();
}

void AliOADBContainer::List()
{
  // List Objects
  for (Int_t i = 0; i < fEntries; i++) {
    printf("Lower %5d Upper %5d \n", fLowerLimits[i], fUpperLimits[i]);
    (fArray->At(i))->Dump();
  }
}
