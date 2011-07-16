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

AliAODCaloCells::AliAODCaloCells() : AliVCaloCells(), fNCells(0), fCellNumber(0), fAmplitude(0), fIsSorted(kTRUE), fType(kUndef)
{
  // default constructor
}

AliAODCaloCells::AliAODCaloCells(const char* name, const char* title, VCells_t ttype) :
    AliVCaloCells(name, title), fNCells(0), fCellNumber(0), fAmplitude(0), fIsSorted(kTRUE), fType(ttype)
{
  //constructor
}

AliAODCaloCells::AliAODCaloCells(const AliAODCaloCells& cells) :
    AliVCaloCells(cells),
    fNCells(cells.fNCells),
    fCellNumber(0),
    fAmplitude(0),
    fIsSorted(cells.fIsSorted),
    fType(cells.fType)
{
// Copy constructor
  fCellNumber = new Short_t[fNCells];
  fAmplitude  = new Double32_t[fNCells]; 
  
  for (Int_t i = 0; i < fNCells; i++) {
    fCellNumber[i]    = cells.fCellNumber[i];
    fAmplitude[i]     = cells.fAmplitude[i];
  }
}

AliAODCaloCells& AliAODCaloCells::operator=(const AliAODCaloCells& cells)
{
    // Assignment operator
  if(&cells == this) return *this;
  delete [] fCellNumber;
  delete [] fAmplitude;

  fNCells = cells.fNCells;

  fCellNumber = new Short_t[fNCells];
  fAmplitude  = new Double32_t[fNCells];

  for (Int_t i = 0; i < fNCells; i++) {
    fCellNumber[i]    = cells.fCellNumber[i];
    fAmplitude[i]     = cells.fAmplitude[i];
  }
  fIsSorted = cells.fIsSorted;
  fType = cells.fType;
  
  SetName(cells.GetName()) ; 
  SetTitle(cells.GetTitle()) ; 
  return *this;
}

AliAODCaloCells::~AliAODCaloCells()
{
  // destructor

  DeleteContainer();
}

void AliAODCaloCells::Clear(const Option_t*)
{
  // clear
  
  DeleteContainer();
}

void AliAODCaloCells::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent
  
  if(this==&obj)return;
  AliAODCaloCells *robj = dynamic_cast<AliAODCaloCells*>(&obj);
  if(!robj)return; // not an AliAODCaloCells
  *robj = *this;
  
}

AliVCaloCells *AliAODCaloCells::CopyCaloCells(Bool_t all = kTRUE) const {
  
  // copy the calo cells into a new object. If option all=FALSE, just the object type, 
  // for mixing
  
  AliVCaloCells *obj =  new AliAODCaloCells();
  
  if(all){
    obj->SetName (GetName()) ; 
    obj->SetTitle(GetTitle()) ; 
    obj->SetType (GetType()) ; 
  
    obj->SetNumberOfCells(fNCells);
    for (Short_t i = 0; i < fNCells; i++) 
      obj->SetCell(i,fCellNumber[i],fAmplitude[i],-1);
  }
  
  return obj;
  
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
  fAmplitude  = new Double32_t[fNCells];

  // set to zero
  for(int i = 0;i<fNCells;++i){
    fAmplitude[i] = fCellNumber[i] = 0 ;
  }
}

void AliAODCaloCells::DeleteContainer()
{
  // deletes allocated memory

  if (fCellNumber)
  {
    delete[] fCellNumber;
    fCellNumber = NULL;
  }

  if (fAmplitude)
  {
    delete[] fAmplitude;
    fAmplitude = NULL;
  }
  
  fNCells = 0;
  fIsSorted = kFALSE;
}

void AliAODCaloCells::Sort() 
{
  // sort the cell array by cell number
  
  Int_t *idxArray = new Int_t[fNCells];
  TMath::Sort(fNCells,fCellNumber,idxArray,kFALSE);
  
  Short_t    *newIndex     = new Short_t[fNCells];
  Double32_t *newAmplitude = new Double32_t[fNCells];
  for (Int_t i=0; i < fNCells; i++) {
    newIndex[i]     = fCellNumber[idxArray[i]];
    newAmplitude[i] = fAmplitude[idxArray[i]];
  }
  delete [] fCellNumber;
  delete [] fAmplitude;
  fCellNumber = newIndex;
  fAmplitude  = newAmplitude;
  
  delete [] idxArray;
  
  fIsSorted = kTRUE;
} 

Bool_t AliAODCaloCells::SetCell(Short_t pos, Short_t cellNumber, Double32_t amplitude, Double32_t /*time*/)
{
  // Sets a cell at the given position

  if (pos>=0 && pos < fNCells) {
    fCellNumber[pos] = cellNumber;
    fAmplitude[pos]  = amplitude;
    fIsSorted = kFALSE;
    return kTRUE;
  } else {
    return kFALSE;
  }
}
