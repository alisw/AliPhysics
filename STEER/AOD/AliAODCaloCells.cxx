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

#include "AliAODCaloCells.h"

/// \cond CLASSIMP
ClassImp(AliAODCaloCells) ;
/// \endcond

///
/// Default constructor.
///
//_______________________________________________________
AliAODCaloCells::AliAODCaloCells() : 
    AliVCaloCells(), fNCells(0), fHGLG(0), fCellNumber(0), 
    fAmplitude(0), fTime(0), fEFraction(0), fMCLabel(0),
    fIsSorted(kTRUE), fType(kUndef)
{
}

///
/// Constructor.
///
//_____________________________________________________________________________________
AliAODCaloCells::AliAODCaloCells(const char* name, const char* title, VCells_t ttype) :
    AliVCaloCells(name, title), fNCells(0), fHGLG(0), fCellNumber(0), 
    fAmplitude(0), fTime(0), fEFraction(0), fMCLabel(0),
    fIsSorted(kTRUE), fType(ttype)
{
}

///
/// Copy constructor.
///
//________________________________________________________________
AliAODCaloCells::AliAODCaloCells(const AliAODCaloCells& cells) :
    AliVCaloCells(cells), fNCells(cells.fNCells), fHGLG(0), fCellNumber(0),
    fAmplitude(0), fTime(0), fEFraction(0), fMCLabel(0),
    fIsSorted(cells.fIsSorted), fType(cells.fType)
{
  fCellNumber = new Short_t[fNCells];
  fAmplitude  = new Double32_t[fNCells]; 
  fTime       = new Double32_t[fNCells]; 
  fMCLabel    = new Int_t[fNCells]; 
  fEFraction  = new Double32_t[fNCells]; 

  for (Int_t i = 0; i < fNCells; i++) {
    fCellNumber[i]    = cells.fCellNumber[i];
    fAmplitude[i]     = cells.fAmplitude[i];
    if(cells.fTime)  fTime[i]      = cells.fTime[i];
    if(cells.fMCLabel)  fMCLabel[i]   = cells.fMCLabel[i];
    if(cells.fEFraction)fEFraction[i] = cells.fEFraction[i];    
  }
  if(cells.fHGLG){
    fHGLG       = new Bool_t[fNCells] ;
    for (Int_t i = 0; i < fNCells; i++) {
      fHGLG[i]          = cells.fHGLG[i];     
    }
  }
}

///
/// Assignment operator.
///
//________________________________________________________________________
AliAODCaloCells& AliAODCaloCells::operator=(const AliAODCaloCells& source)
{
  if(this != &source) 
  {
    AliVCaloCells::operator=(source);
    
    if(fNCells != source.fNCells) 
    {
      if(fHGLG)
        delete [] fHGLG ;
      delete [] fCellNumber;
      delete [] fAmplitude;
      delete [] fTime;
      delete [] fMCLabel;
      delete [] fEFraction;

      fNCells     = source.fNCells;
      
      if(source.fHGLG)
        fHGLG       = new Bool_t[fNCells] ;
      fCellNumber = new Short_t[fNCells];
      fAmplitude  = new Double32_t[fNCells];
      fTime       = new Double32_t[fNCells];
      fMCLabel    = new Int_t[fNCells];
      fEFraction  = new Double32_t[fNCells];
    }
    
    if(source.fHGLG)
      memcpy(fCellNumber,source.fHGLG,fNCells*sizeof(Bool_t));
    memcpy(fCellNumber,source.fCellNumber, fNCells*sizeof(Short_t));
    memcpy(fAmplitude, source.fAmplitude,  fNCells*sizeof(Double32_t));
    if(source.fTime      && fTime)      memcpy(fTime,      source.fTime,      fNCells*sizeof(Double32_t));
    if(source.fMCLabel   && fMCLabel)   memcpy(fMCLabel,   source.fMCLabel,   fNCells*sizeof(Int_t));
    if(source.fEFraction && fEFraction) memcpy(fEFraction, source.fEFraction, fNCells*sizeof(Double32_t));

    fIsSorted = source.fIsSorted;
    fType     = source.fType;
  }

  return *this;
  
}

//_________________________________
AliAODCaloCells::~AliAODCaloCells()
{
  // destructor

  DeleteContainer();
}

///
/// Clear arrays.
///
//__________________________________________
void AliAODCaloCells::Clear(const Option_t*)
{  
  DeleteContainer();
}

///
/// This overwrites the virtual TObject::Copy()
/// to allow run time copying without casting
/// in AliAODEvent.
///
//____________________________________________
void AliAODCaloCells::Copy(TObject &obj) const 
{
  if(this==&obj)return;
  AliAODCaloCells *robj = dynamic_cast<AliAODCaloCells*>(&obj);
  if(!robj)return; // not an AliAODCaloCells
  *robj = *this;
}

///
/// Copy the calo cells into a new object. If option all=FALSE, just the object type, 
/// for mixing.
///
//_____________________________________________________________________
AliVCaloCells *AliAODCaloCells::CopyCaloCells(Bool_t all = kTRUE) const 
{  
  AliVCaloCells *obj =  new AliAODCaloCells();
  
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
      
      Float_t time = -1;
      if(fTime) time = fTime[i];
      
      obj->SetCell(i,fCellNumber[i],fAmplitude[i],time,mclabel,efrac);
    }
  }
  
  return obj;
  
}

///
/// Function that creates container to store calorimeter cell data.
///
//___________________________________________________
void AliAODCaloCells::CreateContainer(Short_t nCells)
{
  DeleteContainer();
  
  if (nCells <= 0) 
  {
    fNCells = 0;
    return;
  }

  fNCells = nCells;

  fHGLG       = new Bool_t[fNCells];
  fCellNumber = new Short_t[fNCells];
  fAmplitude  = new Double32_t[fNCells];
  fTime       = new Double32_t[fNCells];
  fMCLabel    = new Int_t[fNCells];
  fEFraction  = new Double32_t[fNCells];

  // set to zero
  for(int i = 0;i<fNCells;++i)
  {
    fHGLG[i]=kFALSE ;
    fAmplitude [i] =  0.; 
    fCellNumber[i] = -1 ; 
    fEFraction [i] =  0.;
    fTime      [i] = -1.;
    fMCLabel   [i] = -1 ;
  }
}

///
/// Deletes allocated memory.
///
//_____________________________________
void AliAODCaloCells::DeleteContainer()
{
  if(fHGLG){
    delete[] fHGLG;
    fHGLG = 0 ;
  }

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

///
/// Sort the cell array by cell number.
///
//__________________________
void AliAODCaloCells::Sort() 
{  
  Int_t *idxArray = new Int_t[fNCells];
  TMath::Sort(fNCells,fCellNumber,idxArray,kFALSE);
  
  Bool_t     *newHGLG =0x0 ;
  if(fHGLG) newHGLG = new Bool_t[fNCells];
  Short_t    *newIndex     = new Short_t[fNCells];
  Double32_t *newAmplitude = new Double32_t[fNCells];
  
  Double32_t *newTime      = 0 ; 
  Int_t      *newMCLabel   = 0 ;
  Double32_t *newEFraction = 0 ; 
  if(fTime)      newTime      = new Double32_t[fNCells];
  if(fMCLabel)   newMCLabel   = new Int_t[fNCells];
  if(fEFraction) newEFraction = new Double32_t[fNCells];
  
  for (Int_t i=0; i < fNCells; i++) 
  {
    newIndex[i]     = fCellNumber[idxArray[i]];
    newAmplitude[i] = fAmplitude [idxArray[i]];
    if(fTime)      newTime[i]      = fTime     [idxArray[i]];
    if(fMCLabel)   newMCLabel[i]   = fMCLabel  [idxArray[i]];
    if(fEFraction) newEFraction[i] = fEFraction[idxArray[i]];  
  }

  if(fHGLG)
  {
    for (Int_t i=0; i < fNCells; i++) 
    {
      newHGLG[i] = fHGLG[idxArray[i]];
    }
    delete [] fHGLG;
  }
  
  delete [] fCellNumber;
  delete [] fAmplitude;
  delete [] fTime;
  delete [] fMCLabel;
  delete [] fEFraction;

  fHGLG       = newHGLG;
  fCellNumber = newIndex;
  fAmplitude  = newAmplitude;
  fTime       = newTime;
  fMCLabel    = newMCLabel;
  fEFraction  = newEFraction;

  delete [] idxArray;
  
  fIsSorted = kTRUE;
} 

///
/// Sets a cell at the given position.
/// \param pos: cell position in array.
/// \param cellNumber: Cell absolute Id. number.
/// \param amplitude: Cell signal (GeV).
/// \param time: Cell time (s).
/// \param mclabel: MC particle index in kine array.
/// \param efrac: Fraction of energy from embedding.
/// \param isHG: bool true if cell is from high gain.
///
//________________________________________________________________________________________
Bool_t AliAODCaloCells::SetCell(Short_t pos,     Short_t cellNumber, Double32_t amplitude, 
                                Double32_t time, Int_t mclabel,    Double32_t efrac, Bool_t isHG)
{
  if (pos>=0 && pos < fNCells) 
  {
    if(fHGLG)
      fHGLG[pos]=isHG ;
    fCellNumber[pos] = cellNumber;
    fAmplitude[pos]  = amplitude;
    
    // note: initialize (can't use memset for non-0 values)
    //       plus sizeof(Double32_t) is 0
    if(!fTime){
      fTime  = new Double32_t[fNCells];
      
      for( Int_t i = 0; i < fNCells; i++ )
        fTime[i] = -1;
    }
    if(!fMCLabel){
      fMCLabel = new Int_t[fNCells];
      
      for( Int_t i = 0; i < fNCells; i++ )
        fMCLabel[i] = -1;
    }
    if(!fEFraction){
      fEFraction = new Double32_t[fNCells];
      
      for( Int_t i = 0; i < fNCells; i++ )
        fEFraction[i] = 0;
    }
    
    fTime[pos]       = time;
    fMCLabel[pos]    = mclabel;
    fEFraction[pos]  = efrac;
    
    fIsSorted        = kFALSE;
    
    return kTRUE;
    
  } else {
    return kFALSE;
  }
}




