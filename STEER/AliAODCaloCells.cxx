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
//     AOD class to store calorimeter cell data
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliAODCaloCells.h"

ClassImp(AliAODCaloCells)

AliAODCaloCells::AliAODCaloCells() : TNamed(), fNCells(0), fCellNumber(0), fAmplitude(0), fIsSorted(kTRUE), fType(kUndef)
{
  // default constructor
}

AliAODCaloCells::AliAODCaloCells(const char* name, const char* title, AODTwrs_t ttype) : TNamed(name, title), fNCells(0), fCellNumber(0), fAmplitude(0), fIsSorted(kTRUE), fType(ttype)
{
  // TNamed constructor
}

AliAODCaloCells::~AliAODCaloCells()
{
  // destructor

  DeleteContainer();
}

void AliAODCaloCells::CreateContainer(Short_t nCells)
{
  // function that creates container to store calorimeter cell data

  DeleteContainer();
  
  if (nCells <= 0) {
    fNCells = 0;
    return;
  }

  fNCells = nCells;

  fCellNumber = new Short_t[fNCells];
  fAmplitude = new Double32_t[fNCells];
}

void AliAODCaloCells::DeleteContainer()
{
  // deletes allocated memory

  if (fCellNumber)
  {
    delete[] fCellNumber;
    fCellNumber = 0;
  }

  if (fAmplitude)
  {
    delete[] fAmplitude;
    fAmplitude = 0;
  }

  fNCells = 0;
  fIsSorted = kFALSE;
}

void AliAODCaloCells::Sort() 
{
  // sort the cell array by cell number
  
  Int_t *idxArray = new Int_t[fNCells];
  TMath::Sort(fNCells,fCellNumber,idxArray,kFALSE);
  
  Short_t *newIndex = new Short_t[fNCells];
  Double32_t *newAmplitude = new Double32_t[fNCells];
  for (Int_t i=0; i < fNCells; i++) {
    newIndex[i] = fCellNumber[idxArray[i]];
    newAmplitude[i] = fAmplitude[idxArray[i]];
  }
  delete [] fCellNumber;
  delete [] fAmplitude;
  fCellNumber = newIndex;
  fAmplitude = newAmplitude;
  
  delete [] idxArray;
  
  fIsSorted = kTRUE;
} 

Bool_t AliAODCaloCells::SetCell(Short_t pos, Short_t cellNumber, Double32_t amplitude)
{
  // Sets a cell at the given position

  if (pos>=0 && pos < fNCells) {
    fCellNumber[pos] = cellNumber;
    fAmplitude[pos] = amplitude;
    fIsSorted = kFALSE;
    return kTRUE;
  } else {
    return kFALSE;
  }
}
