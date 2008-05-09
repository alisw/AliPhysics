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
//     ESD class to store calorimeter cell data
//     Clone of AliAODCaloCells made by Markus Oldenburg, CERN
//     Author: Gustavo Conesa Balbastre INFN-LNF
//-------------------------------------------------------------------------

#include "AliESDCaloCells.h"

ClassImp(AliESDCaloCells)

//_______________________________________________________________________
AliESDCaloCells::AliESDCaloCells() : 
  TNamed(), fNCells(0), fCellNumber(0), fAmplitude(0), fTime(0), fIsSorted(kTRUE), fType(kUndef)
{
  // default constructor
}
//_______________________________________________________________________
 AliESDCaloCells::AliESDCaloCells(const char* name, const char* title, ESDCells_t ttype) : 
   TNamed(name, title), fNCells(0), fCellNumber(0), fAmplitude(0),  fTime(0), fIsSorted(kTRUE), fType(ttype)
 {
   // TNamed constructor
 }

//_______________________________________________________________________
AliESDCaloCells::AliESDCaloCells(const AliESDCaloCells& c) : 
  TNamed(c), fNCells(c.fNCells),  fCellNumber(), fAmplitude(), fTime(), fIsSorted(c.fIsSorted), fType(c.fType)
{
  // copy constructor

  for(Int_t i = 0; i < fNCells; i++){
    fCellNumber[i] = c.fCellNumber[i];
    fAmplitude[i] = c.fAmplitude[i];
    fTime[i] = c.fTime[i];
  }
}

//_______________________________________________________________________
AliESDCaloCells & AliESDCaloCells::operator =(const AliESDCaloCells& source)  
{
  // assignment operator

  if(&source == this) return *this;
  TNamed::operator=(source);

  if(fNCells != source.fNCells){
    DeleteContainer();
    CreateContainer(source.fNCells);
  }

  fNCells = source.fNCells; 
  fIsSorted = source.fIsSorted;
  fType = source.fType;



  for(Int_t i = 0; i < fNCells; i++){
    fCellNumber[i] = source.fCellNumber[i];
    fAmplitude[i] = source.fAmplitude[i];
    fTime[i] = source.fTime[i];
  }

  return *this;

}


void AliESDCaloCells::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliESDCaloCells *robj = dynamic_cast<AliESDCaloCells*>(&obj);
  if(!robj)return; // not an AliESDCaloCells
  *robj = *this;

}

//_______________________________________________________________________
AliESDCaloCells::~AliESDCaloCells()
{
  // destructor

  DeleteContainer();
}

//_______________________________________________________________________
void AliESDCaloCells::CreateContainer(Short_t nCells)
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
  fTime = new Double32_t[fNCells];
}

//_______________________________________________________________________
void AliESDCaloCells::DeleteContainer()
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

  if (fTime)
  {
    delete[] fTime;
    fTime = 0;
  }

  fNCells = 0;
  fIsSorted = kFALSE;
}

//_______________________________________________________________________
void AliESDCaloCells::Sort() 
{
  // sort the cell array by cell number
  
  Int_t *idxArray = new Int_t[fNCells];
  TMath::Sort(fNCells,fCellNumber,idxArray,kFALSE);
  
  Short_t *newIndex = new Short_t[fNCells];
  Double32_t *newAmplitude = new Double32_t[fNCells];
  Double32_t *newTime = new Double32_t[fNCells];
  for (Int_t i=0; i < fNCells; i++) {
    newIndex[i] = fCellNumber[idxArray[i]];
    newAmplitude[i] = fAmplitude[idxArray[i]];
    newTime[i] = fTime[idxArray[i]];
  }
  delete [] fCellNumber;
  delete [] fAmplitude;
  delete [] fTime;
  fCellNumber = newIndex;
  fAmplitude = newAmplitude;
  fTime = newTime;

  delete [] idxArray;
  
  fIsSorted = kTRUE;
} 

//_______________________________________________________________________
Bool_t AliESDCaloCells::SetCell(Short_t pos, Short_t cellNumber, Double32_t amplitude, Double32_t  time)
{
  // Sets a cell at the given position

  if (pos>=0 && pos < fNCells) {
    fCellNumber[pos] = cellNumber;
    fAmplitude[pos] = amplitude;
    fTime[pos] = time;
    fIsSorted = kFALSE;
    return kTRUE;
  } else {
    return kFALSE;
  }
}
