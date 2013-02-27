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

//_____________________________________________________
AliESDCaloCells::AliESDCaloCells() : 
  AliVCaloCells(), fNCells(0), fCellNumber(0), 
  fAmplitude(0), fTime(0), fEFraction(0), fMCLabel(0),
  fIsSorted(kTRUE), fType(kUndef)
{
  // default constructor
}
//______________________________________________________________________________________
 AliESDCaloCells::AliESDCaloCells(const char* name, const char* title, VCells_t ttype) : 
   AliVCaloCells(name, title), fNCells(0), fCellNumber(0), 
   fAmplitude(0),  fTime(0), fEFraction(0), fMCLabel(0),
   fIsSorted(kTRUE), fType(ttype)
{
   // AliVCaloCells constructor
 }

//__________________________________________________________
AliESDCaloCells::AliESDCaloCells(const AliESDCaloCells& c) : 
  AliVCaloCells(c), fNCells(c.fNCells),  fCellNumber(0), 
  fAmplitude(0), fTime(0), fEFraction(0), fMCLabel(0),
  fIsSorted(c.fIsSorted), fType(c.fType)
{
  // copy constructor

  fCellNumber = new Short_t[fNCells];
  fAmplitude  = new Double32_t[fNCells];
  fTime       = new Double32_t[fNCells];
  fMCLabel    = new Int_t[fNCells];
  fEFraction  = new Double32_t[fNCells];

  for(Int_t i = 0; i < fNCells; i++)
  {
    fCellNumber[i]    = c.fCellNumber[i];
    fAmplitude[i]     = c.fAmplitude[i];
    fTime[i]          = c.fTime[i];
    if(c.fMCLabel)  fMCLabel[i]   = c.fMCLabel[i];
    if(c.fEFraction)fEFraction[i] = c.fEFraction[i];
  }
}

//__________________________________________________________________________
AliESDCaloCells & AliESDCaloCells::operator =(const AliESDCaloCells& source)  
{
  // assignment operator

  if(this != &source)
  {
    AliVCaloCells::operator=(source);
    
    if(fNCells != source.fNCells) 
    {
      delete [] fCellNumber;
      delete [] fAmplitude;
      delete [] fTime;
      delete [] fMCLabel;
      delete [] fEFraction;
      
      fNCells = source.fNCells;
      
      fCellNumber = new Short_t[fNCells];
      fAmplitude  = new Double32_t[fNCells];
      fTime       = new Double32_t[fNCells];
      fMCLabel    = new Int_t[fNCells];
      fEFraction  = new Double32_t[fNCells];
    }
    
    memcpy(fCellNumber,source.fCellNumber,fNCells*sizeof(Short_t));
    memcpy(fAmplitude, source.fAmplitude, fNCells*sizeof(Double32_t));
    memcpy(fTime,      source.fTime,      fNCells*sizeof(Double32_t));
    if(source.fMCLabel   && fMCLabel)   memcpy(fMCLabel,   source.fMCLabel,   fNCells*sizeof(Int_t));
    if(source.fEFraction && fEFraction) memcpy(fEFraction, source.fEFraction, fNCells*sizeof(Double32_t));

    fIsSorted = source.fIsSorted;
    fType = source.fType;
    
  }

  return *this;
}

//____________________________________________
void AliESDCaloCells::Copy(TObject &obj) const 
{
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliESDCaloCells *robj = dynamic_cast<AliESDCaloCells*>(&obj);
  if(!robj)return; // not an AliESDCaloCells
  *robj = *this;

}

//______________________________________________________________________
AliVCaloCells* AliESDCaloCells::CopyCaloCells(Bool_t all = kTRUE) const
{  
  // copy the calo cells into a new object. If option all=FALSE, just the object type, 
  // for mixing
  
  AliVCaloCells *obj = new AliESDCaloCells();
  
  if(all){
    obj->SetName (GetName()) ; 
    obj->SetTitle(GetTitle()) ; 
    obj->SetType (GetType()) ; 
    
    obj->SetNumberOfCells(fNCells);
    for (Short_t i = 0; i < fNCells; i++) 
    {
      Int_t mclabel = -1;
      if(fMCLabel) mclabel = fMCLabel[i];
      
      Float_t efrac = 0.;
      if(fEFraction) efrac = fEFraction[i];
      
      obj->SetCell(i,fCellNumber[i],fAmplitude[i],fTime[i],mclabel,efrac);
    }
  }

  return obj;
}


//_________________________________
AliESDCaloCells::~AliESDCaloCells()
{
  // destructor

  DeleteContainer();
}

//__________________________________________
void AliESDCaloCells::Clear(const Option_t*)
{
  // clear
  
  DeleteContainer();
}


//___________________________________________________
void AliESDCaloCells::CreateContainer(Short_t nCells)
{
  // function that creates container to store calorimeter cell data

  DeleteContainer();
  
  if (nCells <= 0) 
  {
    fNCells = 0;
    return;
  }

  fNCells = nCells;

  fCellNumber = new Short_t[fNCells];
  fAmplitude  = new Double32_t[fNCells];
  fTime       = new Double32_t[fNCells];
  fMCLabel    = new Int_t[fNCells];
  fEFraction  = new Double32_t[fNCells];

  // set to zero
  for(int i = 0;i<fNCells;++i)
  {
    fAmplitude [i] =  0.; 
    fCellNumber[i] = -1 ; 
    fEFraction [i] =  0.;
    fTime      [i] = -1.;
    fMCLabel   [i] = -1 ;
  }
}

//_____________________________________
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
    fAmplitude = NULL;
  }

  if (fTime)
  {
    delete[] fTime;
    fTime = NULL;
  }
  
  if (fMCLabel)
  {
    delete[] fMCLabel;
    fMCLabel = NULL;
  }
  
  if (fEFraction)
  {
    delete[] fEFraction;
    fEFraction = NULL;
  }
  
  fNCells = 0;
  fIsSorted = kFALSE;
  
}

//__________________________
void AliESDCaloCells::Sort() 
{
  // sort the cell array by cell number
  
  Int_t *idxArray = new Int_t[fNCells];
  TMath::Sort(fNCells,fCellNumber,idxArray,kFALSE);
  
  Short_t    *newIndex     = new Short_t[fNCells];
  Double32_t *newAmplitude = new Double32_t[fNCells];
  Double32_t *newTime      = new Double32_t[fNCells];
  
  Int_t    *newMCLabel   = 0 ;
  Double32_t *newEFraction = 0 ; 
  if(fMCLabel)   newMCLabel   = new Int_t[fNCells];
  if(fEFraction) newEFraction = new Double32_t[fNCells];
    
  for (Int_t i=0; i < fNCells; i++) 
  {
    newIndex[i]     = fCellNumber[idxArray[i]];
    newAmplitude[i] = fAmplitude [idxArray[i]];
    newTime[i]      = fTime      [idxArray[i]];
    if(fMCLabel)   newMCLabel[i]   = fMCLabel  [idxArray[i]];
    if(fEFraction) newEFraction[i] = fEFraction[idxArray[i]];
  }
  
  delete [] fCellNumber;
  delete [] fAmplitude;
  delete [] fTime;
  delete [] fMCLabel;
  delete [] fEFraction;

  fCellNumber = newIndex;
  fAmplitude  = newAmplitude;
  fTime       = newTime;
  if(fMCLabel)   fMCLabel    = newMCLabel;
  if(fEFraction) fEFraction  = newEFraction;

  delete [] idxArray;
  
  fIsSorted = kTRUE;
} 

//________________________________________________________________________________________
Bool_t AliESDCaloCells::SetCell(Short_t pos,     Short_t cellNumber, Double32_t amplitude,  
                                Double32_t time, Int_t mclabel,    Double32_t efrac)
{
  // Sets a cell at the given position

  if (pos>=0 && pos < fNCells) 
  {
    fCellNumber[pos] = cellNumber;
    fAmplitude[pos]  = amplitude;
    fTime[pos]       = time;
    
    // note: initialize (can't use memset for non-0 values)
    //       plus sizeof(Double32_t) is 0
    if(!fMCLabel){
      fMCLabel   = new Int_t[fNCells];

      for( Int_t i = 0; i < fNCells; i++ )
        fMCLabel[i] = -1;
    }
    if(!fEFraction){
      fEFraction = new Double32_t[fNCells];

      for( Int_t i = 0; i < fNCells; i++ )
        fEFraction[i] = 0;
    }

    fMCLabel[pos]    = mclabel;
    fEFraction[pos]  = efrac;
    
    fIsSorted = kFALSE;
    
    return kTRUE;
    
  } else {
    return kFALSE;
  }
}
