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

/* $Id: AliTHn.cxx 20164 2007-08-14 15:31:50Z morsch $ */

//
// this storage container is optimized for small memory usage
//   under/over flow bins do not exist
//   sumw2 structure is float only and only create when the weight != 1
// all histogram functionality (projections, axis, ranges, etc) are taken from THnSparse by propagating 
// the information up into the THnSparse structure (in the ananalysis *after* data processing and merging)
//
// the derivation from THnSparse is obviously against many OO rules. correct would be a common baseclass of THnSparse and THn.
// 
// Author: Jan Fiete Grosse-Oetringhaus

#include "AliTHn.h"
#include "TList.h"
#include "TCollection.h"
#include "AliLog.h"
#include "TArrayF.h"
#include "THnSparse.h"
#include "TMath.h"

ClassImp(AliTHn)

AliTHn::AliTHn() : 
  AliCFContainer(),
  fNBins(0),
  fNVars(0),
  fNSteps(0),
  fValues(0),
  fSumw2(0)
{
  // Constructor
}

AliTHn::AliTHn(const Char_t* name, const Char_t* title,const Int_t nSelStep, const Int_t nVarIn, const Int_t* nBinIn) : 
  AliCFContainer(name, title, nSelStep, nVarIn, nBinIn),
  fNBins(0),
  fNVars(nVarIn),
  fNSteps(nSelStep),
  fValues(0),
  fSumw2(0)
{
  // Constructor

  fNBins = 1;
  for (Int_t i=0; i<fNVars; i++)
    fNBins *= nBinIn[i];
  
  Init();
}

void AliTHn::Init()
{
  // initialize
  
  fValues = new TArrayF*[fNSteps];
  fSumw2 = new TArrayF*[fNSteps];
  
  for (Int_t i=0; i<fNSteps; i++)
  {
    fValues[i] = 0;
    fSumw2[i] = 0;
  }
} 

AliTHn::AliTHn(const AliTHn &c) :
  AliCFContainer(),
  fNBins(0),
  fNVars(0),
  fNSteps(0),
  fValues(0),
  fSumw2(0)
{
  //
  // AliTHn copy constructor
  //

  ((AliTHn &) c).Copy(*this);
}

AliTHn::~AliTHn()
{
  // Destructor
  
  DeleteContainers();
  
  if (fValues)
  {
    delete fValues;
    fValues = 0;
  }

  if (fSumw2)
  {
    delete fSumw2;
    fSumw2 = 0;
  }
}

void AliTHn::DeleteContainers()
{
  // delete data containers
  
  for (Int_t i=0; i<fNSteps; i++)
  {
    if (fValues && fValues[i])
    {
      delete fValues[i];
      fValues[i] = 0;
    }
    
    if (fSumw2 && fSumw2[i])
    {
      delete fSumw2[i];
      fSumw2[i] = 0;
    }
  }
}

//____________________________________________________________________
AliTHn &AliTHn::operator=(const AliTHn &c)
{
  // assigment operator

  if (this != &c)
    ((AliTHn &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
void AliTHn::Copy(TObject& c) const
{
  // copy function

  AliTHn& target = (AliTHn &) c;
  
  AliCFContainer::Copy(target);
  
  target.fNSteps = fNSteps;
  target.fNBins = fNBins;
  target.fNVars = fNVars;
  
  target.Init();

  for (Int_t i=0; i<fNSteps; i++)
  {
    if (fValues[i])
      target.fValues[i] = new TArrayF(*(fValues[i]));
    else
      target.fValues[i] = 0;
    
    if (fSumw2[i])
      target.fSumw2[i] = new TArrayF(*(fSumw2[i]));
    else
      target.fSumw2[i] = 0;
  }
}

//____________________________________________________________________
Long64_t AliTHn::Merge(TCollection* list)
{
  // Merge a list of AliTHn objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;
  
  AliCFContainer::Merge(list);

  TIterator* iter = list->MakeIterator();
  TObject* obj;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliTHn* entry = dynamic_cast<AliTHn*> (obj);
    if (entry == 0) 
      continue;

    for (Int_t i=0; i<fNSteps; i++)
    {
      if (entry->fValues[i])
      {
	if (!fValues[i])
	  fValues[i] = new TArrayF(fNBins);
      
	for (Long64_t l = 0; l<fNBins; l++)
	  fValues[i]->GetArray()[l] += entry->fValues[i]->GetArray()[l];
      }

      if (entry->fSumw2[i])
      {
	if (!fSumw2[i])
	  fSumw2[i] = new TArrayF(fNBins);
      
	for (Long64_t l = 0; l<fNBins; l++)
	  fSumw2[i]->GetArray()[l] += entry->fSumw2[i]->GetArray()[l];
      }
    }
    
    count++;
  }

  return count+1;
}

void AliTHn::Fill(const Double_t *var, Int_t istep, Double_t weight)
{
  // fills an entry

  // calculate global bin index
  Long64_t bin = 0;
  for (Int_t i=0; i<fNVars; i++)
  {
    bin *= GetAxis(i, 0)->GetNbins();
    
    Int_t tmpBin = GetAxis(i, 0)->FindBin(var[i]);
//     Printf("%d", tmpBin);
    // under/overflow not supported
    if (tmpBin < 1 || tmpBin > GetAxis(i, 0)->GetNbins())
      return;
    
    // bins start from 0 here
    bin += tmpBin - 1;
//     Printf("%lld", bin);
  }

  if (!fValues[istep])
  {
    fValues[istep] = new TArrayF(fNBins);
    AliInfo(Form("Created values container for step %d", istep));
  }

  if (weight != 1)
  {
    // initialize with already filled entries (which have been filled with weight == 1), in this case fSumw2 := fValues
    if (!fSumw2[istep])
    {
      fSumw2[istep] = new TArrayF(*fValues[istep]);
      AliInfo(Form("Created sumw2 container for step %d", istep));
    }
  }

  fValues[istep]->GetArray()[bin] += weight;
  if (fSumw2[istep])
    fSumw2[istep]->GetArray()[bin] += weight * weight;
  
//   Printf("%f", fValues[istep][bin]);
  
  // debug
//   AliCFContainer::Fill(var, istep, weight);
}

Long64_t AliTHn::GetGlobalBinIndex(const Int_t* binIdx)
{
  // calculates global bin index
  // binIdx contains TAxis bin indexes
  // here bin count starts at 0 because we do not have over/underflow bins
  
  Long64_t bin = 0;
  for (Int_t i=0; i<fNVars; i++)
  {
    bin *= GetAxis(i, 0)->GetNbins();
    bin += binIdx[i] - 1;
  }

  return bin;
}

void AliTHn::FillParent()
{
  // fills the information stored in the buffer in this class into the baseclass containers
  
  for (Int_t i=0; i<fNSteps; i++)
  {
    if (!fValues[i])
      continue;
      
    Float_t* source = fValues[i]->GetArray();
    // if fSumw2 is not stored, the sqrt of the number of bin entries in source is filled below; otherwise we use fSumw2
    Float_t* sourceSumw2 = source;
    if (fSumw2[i])
      sourceSumw2 = fSumw2[i]->GetArray();
    
    THnSparse* target = GetGrid(i)->GetGrid();
    
    Int_t* binIdx = new Int_t[fNVars];
    for (Int_t j=0; j<fNVars; j++)
      binIdx[j] = 1;
    
    Long64_t count = 0;
    
    while (1)
    {
//       for (Int_t j=0; j<fNVars; j++)
// 	printf("%d ", binIdx[j]);
      
      Long64_t globalBin = GetGlobalBinIndex(binIdx);
//       Printf(" --> %lld", globalBin);
      
      if (source[globalBin] != 0)
      {
	target->SetBinContent(binIdx, source[globalBin]);
	target->SetBinError(binIdx, TMath::Sqrt(sourceSumw2[globalBin]));
	
	count++;
      }
      
      binIdx[fNVars-1]++;
      
      for (Int_t j=fNVars-1; j>0; j--)
      {
	if (binIdx[j] > target->GetAxis(j)->GetNbins())
	{
	  binIdx[j] = 1;
	  binIdx[j-1]++;
	}
      }
      
      if (binIdx[0] > target->GetAxis(0)->GetNbins())
	break;
    }
    
    AliInfo(Form("Step %d: copied %lld entries out of %lld bins", i, count, GetGlobalBinIndex(binIdx)));
  }
}
