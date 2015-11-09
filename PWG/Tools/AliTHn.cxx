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

/* $Id: AliTHnT.cxx 20164 2007-08-14 15:31:50Z morsch $ */

// Use AliTHnT instead of AliCFContainer and your memory consumption will be drastically reduced
// As AliTHnT derives from AliCFContainer, you can just replace your current AliCFContainer object by AliTHnT
// Once you have the merged output, call FillParent() and you can use AliCFContainer as usual
//
// this storage container is optimized for small memory usage
//   under/over flow bins do not exist
//   sumw2 structure is float only and only create when the weight != 1
// all histogram functionality (projections, axis, ranges, etc) are taken from THnSparse by propagating 
// the information up into the THnSparse structure (in the ananalysis *after* data processing and merging)
//
// the derivation from THnSparse is obviously against many OO rules. correct would be a common baseclass of THnSparse and THn.
//
// Templated version allows also the use of double as storage container
// 
// Author: Jan Fiete Grosse-Oetringhaus

#include "AliTHn.h"
#include "TList.h"
#include "TCollection.h"
#include "AliLog.h"
#include "TArrayF.h"
#include "TArrayD.h"
#include "THnSparse.h"
#include "TMath.h"

templateClassImp(AliTHnT)

template <class TemplateArray, typename TemplateType>
AliTHnT<TemplateArray, TemplateType>::AliTHnT() : 
  AliTHnBase(),
  fNBins(0),
  fNVars(0),
  fNSteps(0),
  fValues(0),
  fSumw2(0),
  axisCache(0),
  fNbinsCache(0),
  fLastVars(0),
  fLastBins(0)
{
  // Constructor
}

template <class TemplateArray, typename TemplateType>
AliTHnT<TemplateArray, TemplateType>::AliTHnT(const Char_t* name, const Char_t* title,const Int_t nSelStep, const Int_t nVarIn, const Int_t* nBinIn) : 
  AliTHnBase(name, title, nSelStep, nVarIn, nBinIn),
  fNBins(0),
  fNVars(nVarIn),
  fNSteps(nSelStep),
  fValues(0),
  fSumw2(0),
  axisCache(0),
  fNbinsCache(0),
  fLastVars(0),
  fLastBins(0)
{
  // Constructor

  fNBins = 1;
  for (Int_t i=0; i<fNVars; i++)
    fNBins *= nBinIn[i];
  
  Init();
}

template <class TemplateArray, typename TemplateType>
void AliTHnT<TemplateArray, TemplateType>::AliTHnT::Init()
{
  // initialize
  
  fValues = new TemplateArray*[fNSteps];
  fSumw2 = new TemplateArray*[fNSteps];
  
  for (Int_t i=0; i<fNSteps; i++)
  {
    fValues[i] = 0;
    fSumw2[i] = 0;
  }
} 

template <class TemplateArray, typename TemplateType>
AliTHnT<TemplateArray, TemplateType>::AliTHnT(const AliTHnT &c) :
  AliTHnBase(c),
  fNBins(c.fNBins),
  fNVars(c.fNVars),
  fNSteps(c.fNSteps),
  fValues(new TemplateArray*[c.fNSteps]),
  fSumw2(new TemplateArray*[c.fNSteps]),
  axisCache(0),
  fNbinsCache(0),
  fLastVars(0),
  fLastBins(0)
{
  //
  // AliTHnT copy constructor
  //

  memset(fValues,0,fNSteps*sizeof(TemplateArray*));
  memset(fSumw2,0,fNSteps*sizeof(TemplateArray*));

  for (Int_t i=0; i<fNSteps; i++) {
    if (c.fValues[i]) fValues[i] = new TemplateArray(*(c.fValues[i]));
    if (c.fSumw2[i])  fSumw2[i]  = new TemplateArray(*(c.fSumw2[i]));
  }

}

template <class TemplateArray, typename TemplateType>
AliTHnT<TemplateArray, TemplateType>::~AliTHnT()
{
  // Destructor
  
  DeleteContainers();
  
  delete[] fValues;
  delete[] fSumw2;
  delete[] axisCache;
  delete[] fNbinsCache;
  delete[] fLastVars;
  delete[] fLastBins;
}

template <class TemplateArray, typename TemplateType>
void AliTHnT<TemplateArray, TemplateType>::DeleteContainers()
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
template <class TemplateArray, typename TemplateType>
AliTHnT<TemplateArray, TemplateType> &AliTHnT<TemplateArray, TemplateType>::operator=(const AliTHnT<TemplateArray, TemplateType> &c)
{
  // assigment operator

  if (this != &c) {
    AliCFContainer::operator=(c);
    fNBins=c.fNBins;
    fNVars=c.fNVars;
    if(fNSteps) {
      for(Int_t i=0; i< fNSteps; ++i) {
	delete fValues[i];
	delete fSumw2[i];
      }
      delete [] fValues;
      delete [] fSumw2;
    }
    fNSteps=c.fNSteps;
    if(fNSteps) {
      fValues=new TemplateArray*[fNSteps];
      fSumw2=new TemplateArray*[fNSteps];
      memset(fValues,0,fNSteps*sizeof(TemplateArray*));
      memset(fSumw2,0,fNSteps*sizeof(TemplateArray*));

      for (Int_t i=0; i<fNSteps; i++) {
	if (c.fValues[i]) fValues[i] = new TemplateArray(*(c.fValues[i]));
	if (c.fSumw2[i])  fSumw2[i]  = new TemplateArray(*(c.fSumw2[i]));
      }
    } else {
      fValues = 0;
      fSumw2 = 0;
    }
    delete [] axisCache;
    axisCache = new TAxis*[fNVars];
    memcpy(axisCache, c.axisCache, fNVars*sizeof(TAxis*));
  }
  return *this;
}

//____________________________________________________________________
template <class TemplateArray, typename TemplateType>
void AliTHnT<TemplateArray, TemplateType>::Copy(TObject& c) const
{
  // copy function

  AliTHnT& target = (AliTHnT &) c;
  
  AliCFContainer::Copy(target);
  
  target.fNSteps = fNSteps;
  target.fNBins = fNBins;
  target.fNVars = fNVars;
  
  target.Init();

  for (Int_t i=0; i<fNSteps; i++)
  {
    if (fValues[i])
      target.fValues[i] = new TemplateArray(*(fValues[i]));
    else
      target.fValues[i] = 0;
    
    if (fSumw2[i])
      target.fSumw2[i] = new TemplateArray(*(fSumw2[i]));
    else
      target.fSumw2[i] = 0;
  }
}

//____________________________________________________________________
template <class TemplateArray, typename TemplateType>
Long64_t AliTHnT<TemplateArray, TemplateType>::Merge(TCollection* list)
{
  // Merge a list of AliTHnT objects with this (needed for
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
    
    AliTHnT* entry = dynamic_cast<AliTHnT*> (obj);
    if (entry == 0) 
      continue;

    for (Int_t i=0; i<fNSteps; i++)
    {
      if (entry->fValues[i])
      {
	if (!fValues[i])
	  fValues[i] = new TemplateArray(fNBins);
      
	for (Long64_t l = 0; l<fNBins; l++)
	  fValues[i]->GetArray()[l] += entry->fValues[i]->GetArray()[l];
      }

      if (entry->fSumw2[i])
      {
	if (!fSumw2[i])
	  fSumw2[i] = new TemplateArray(fNBins);
      
	for (Long64_t l = 0; l<fNBins; l++)
	  fSumw2[i]->GetArray()[l] += entry->fSumw2[i]->GetArray()[l];
      }
    }
    
    count++;
  }

  return count+1;
}

template <class TemplateArray, typename TemplateType>
void AliTHnT<TemplateArray, TemplateType>::Fill(const Double_t *var, Int_t istep, Double_t weight)
{
  // fills an entry

  // fill axis cache
  if (!axisCache)
  {
    axisCache = new TAxis*[fNVars];
    fNbinsCache = new Int_t[fNVars];
    for (Int_t i=0; i<fNVars; i++)
    {
      axisCache[i] = GetAxis(i, 0);
      fNbinsCache[i] = axisCache[i]->GetNbins();
    }
    
    fLastVars = new Double_t[fNVars];
    fLastBins = new Int_t[fNVars];
    
    // initial values to prevent checking for 0 below
    for (Int_t i=0; i<fNVars; i++)
    {
      fLastBins[i] = axisCache[i]->FindBin(var[i]);
      fLastVars[i] = var[i];
    }
  }
  
  // calculate global bin index
  Long64_t bin = 0;
  for (Int_t i=0; i<fNVars; i++)
  {
    bin *= fNbinsCache[i];
    
    Int_t tmpBin = 0;
    if (fLastVars[i] == var[i])
      tmpBin = fLastBins[i];
    else
    {
      tmpBin = axisCache[i]->FindBin(var[i]);
      fLastBins[i] = tmpBin;
      fLastVars[i] = var[i];
    }
    //Printf("%d", tmpBin);

    // under/overflow not supported
    if (tmpBin < 1 || tmpBin > fNbinsCache[i])
      return;
    
    // bins start from 0 here
    bin += tmpBin - 1;
//     Printf("%lld", bin);
  }

  if (!fValues[istep])
  {
    fValues[istep] = new TemplateArray(fNBins);
    AliInfo(Form("Created values container for step %d", istep));
  }

  if (weight != 1)
  {
    // initialize with already filled entries (which have been filled with weight == 1), in this case fSumw2 := fValues
    if (!fSumw2[istep])
    {
      fSumw2[istep] = new TemplateArray(*fValues[istep]);
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

template <class TemplateArray, typename TemplateType>
Long64_t AliTHnT<TemplateArray, TemplateType>::GetGlobalBinIndex(const Int_t* binIdx)
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

template <class TemplateArray, typename TemplateType>
void AliTHnT<TemplateArray, TemplateType>::FillContainer(AliCFContainer* cont)
{
  // fills the information stored in the buffer in this class into the container <cont>
  
  for (Int_t i=0; i<fNSteps; i++)
  {
    if (!fValues[i])
      continue;
      
    TemplateType* source = fValues[i]->GetArray();
    // if fSumw2 is not stored, the sqrt of the number of bin entries in source is filled below; otherwise we use fSumw2
    TemplateType* sourceSumw2 = source;
    if (fSumw2[i])
      sourceSumw2 = fSumw2[i]->GetArray();
    
    THnSparse* target = cont->GetGrid(i)->GetGrid();
    
    Int_t* binIdx = new Int_t[fNVars];
    Int_t* nBins  = new Int_t[fNVars];
    for (Int_t j=0; j<fNVars; j++)
    {
      binIdx[j] = 1;
      nBins[j] = target->GetAxis(j)->GetNbins();
    }
    
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
	if (binIdx[j] > nBins[j])
	{
	  binIdx[j] = 1;
	  binIdx[j-1]++;
	}
      }
      
      if (binIdx[0] > nBins[0])
	break;
    }
    
    AliInfo(Form("Step %d: copied %lld entries out of %lld bins", i, count, GetGlobalBinIndex(binIdx)));

    delete[] binIdx;
    delete[] nBins;
  }  
}

template <class TemplateArray, typename TemplateType>
void AliTHnT<TemplateArray, TemplateType>::FillParent()
{
  // fills the information stored in the buffer in this class into the baseclass containers
  
  FillContainer(this);
}

template <class TemplateArray, typename TemplateType>
void AliTHnT<TemplateArray, TemplateType>::ReduceAxis()
{
  // "removes" one axis by summing over the axis and putting the entry to bin 1
  // TODO presently only implemented for the last axis
  
  Int_t axis = fNVars-1;
  
  for (Int_t i=0; i<fNSteps; i++)
  {
    if (!fValues[i])
      continue;
      
    TemplateType* source = fValues[i]->GetArray();
    TemplateType* sourceSumw2 = 0;
    if (fSumw2[i])
      sourceSumw2 = fSumw2[i]->GetArray();
    
    THnSparse* target = GetGrid(i)->GetGrid();
    
    Int_t* binIdx = new Int_t[fNVars];
    Int_t* nBins  = new Int_t[fNVars];
    for (Int_t j=0; j<fNVars; j++)
    {
      binIdx[j] = 1;
      nBins[j] = target->GetAxis(j)->GetNbins();
    }
    
    Long64_t count = 0;

    while (1)
    {
      // sum over axis <axis>
      TemplateType sumValues = 0;
      TemplateType sumSumw2 = 0;
      for (Int_t j=1; j<=nBins[axis]; j++)
      {
	binIdx[axis] = j;
	Long64_t globalBin = GetGlobalBinIndex(binIdx);
	sumValues += source[globalBin];
	source[globalBin] = 0;

	if (sourceSumw2)
	{
	  sumSumw2 += sourceSumw2[globalBin];
	  sourceSumw2[globalBin] = 0;
	}
      }
      binIdx[axis] = 1;
	
      Long64_t globalBin = GetGlobalBinIndex(binIdx);
      source[globalBin] = sumValues;
      if (sourceSumw2)
	sourceSumw2[globalBin] = sumSumw2;

      count++;

      // next bin
      binIdx[fNVars-2]++;
      
      for (Int_t j=fNVars-2; j>0; j--)
      {
	if (binIdx[j] > nBins[j])
	{
	  binIdx[j] = 1;
	  binIdx[j-1]++;
	}
      }
      
      if (binIdx[0] > nBins[0])
	break;
    }
    
    AliInfo(Form("Step %d: reduced %lld bins to %lld entries", i, GetGlobalBinIndex(binIdx), count));

    delete[] binIdx;
    delete[] nBins;
  }
}

template class AliTHnT<TArrayF, Float_t>;
template class AliTHnT<TArrayD, Double_t>;
