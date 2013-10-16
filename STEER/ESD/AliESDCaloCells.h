#ifndef ALIESDCALOCELLS_H
#define ALIESDCALOCELLS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
/* $Log $ */

//-------------------------------------------------------------------------
//     ESD class to store calorimeter cell data
//     Clone of AliAODCaloCells made by Markus Oldenburg, CERN
//     Author: Gustavo Conesa Balbastre INFN-LNF
//
//-------------------------------------------------------------------------


#include <AliVCaloCells.h>
#include <TMath.h>

class AliESDCaloCells : public AliVCaloCells
{
 public:

  AliESDCaloCells();
  AliESDCaloCells(const char* name, const char* title, VCells_t ttype=kUndef);
  AliESDCaloCells(const AliESDCaloCells & cells);
  AliESDCaloCells & operator=(const AliESDCaloCells& source);
  virtual ~AliESDCaloCells();
  
  virtual AliVCaloCells * CopyCaloCells(Bool_t all) const;
  virtual void    Copy(TObject &obj) const;
  void            Clear(const Option_t*);
  void            CreateContainer(Short_t nCells);
  void            DeleteContainer();
  void            Sort();
  
  Bool_t          IsEMCAL()  const { return (fType == kEMCALCell); }
  Bool_t          IsPHOS()   const { return (fType == kPHOSCell) ; }
  Char_t          GetType()  const { return  fType               ; }
  void            SetType(Char_t t){ fType = t                   ; }
  
  inline Bool_t   GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude, Double_t &time, Int_t &mclabel,      Double_t &efrac) const;
  Bool_t          SetCell(Short_t pos, Short_t  cellNumber, Double_t  amplitude, Double_t  time, Int_t  mclabel = -1, Double_t  efrac = 0.);
  
  Short_t         GetNumberOfCells() const  { return fNCells ; }
  void            SetNumberOfCells(Int_t n) { fNCells = n    ; }
  
  
  inline Double_t GetCellAmplitude(Short_t cellNumber);
  inline Short_t  GetCellPosition(Short_t cellNumber);
  inline Double_t GetCellTime(Short_t cellNumber);
  
  inline Double_t GetAmplitude(Short_t pos) const;
  inline Double_t GetTime(Short_t pos) const;
  inline Short_t  GetCellNumber(Short_t pos) const;

  // MC & embedding
  inline Int_t    GetCellMCLabel(Short_t cellNumber) ;
  inline Int_t    GetMCLabel(Short_t pos) const ;
  
  inline Double_t GetCellEFraction(Short_t cellNumber) ;
  inline Double_t GetEFraction(Short_t pos) const ;
  
  inline void     SetEFraction    (Short_t pos,         Double32_t efrac) ;
  inline void     SetCellEFraction(Short_t cellNumber,  Double32_t efrac) ;
  
 protected:
  
  Int_t       fNCells;       // Number of cells
  Short_t    *fCellNumber;   //[fNCells] array of cell numbers
  Double32_t *fAmplitude;    //[fNCells][0.,0.,16] array with cell amplitudes (= energy!)
  Double32_t *fTime;         //[fNCells][0.,0.,16] array with cell times
  Double32_t *fEFraction;    //[fNCells][0.,0.,16] array with fraction of MC energy and data - for embedding
  Int_t      *fMCLabel;      //[fNCells] array of MC labels
  Bool_t      fIsSorted;     //! true if cell arrays are sorted by index
  Char_t      fType;         // Cell type

  ClassDef(AliESDCaloCells, 4);
};


Bool_t AliESDCaloCells::GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude, 
                                Double_t & time, Int_t & mclabel, Double_t & efrac) const 
{ 
  if (pos>=0 && pos<fNCells) 
  {
    cellNumber = fCellNumber[pos];
    amplitude  = fAmplitude[pos];
    time       = fTime[pos];
    
    if(fMCLabel)   mclabel = fMCLabel[pos];
    else           mclabel =-1 ; 
    if(fEFraction) efrac   = fEFraction[pos];
    else           efrac   = 0 ;
    
    return kTRUE;
    
  } else 
  {
    return kFALSE;
  }
}

Double_t AliESDCaloCells::GetCellAmplitude(Short_t cellNumber)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }

  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && pos < fNCells && fCellNumber[pos] == cellNumber ) {
    return fAmplitude[pos];
  } else {
    return 0.;
  }
}

Double_t AliESDCaloCells::GetCellTime(Short_t cellNumber)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }

  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && pos < fNCells && fCellNumber[pos] == cellNumber) {
    return fTime[pos];
  } else {
    return -1.;
  }
}

Double_t AliESDCaloCells::GetAmplitude(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    return fAmplitude[pos];
  } else {
    return 0.;
  }
}

Double_t AliESDCaloCells::GetTime(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    return fTime[pos];
  } else {
    return -1.;
  }
}

Short_t AliESDCaloCells::GetCellNumber(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    return fCellNumber[pos];
  } else {
    return fNCells;
  }
}

Short_t AliESDCaloCells::GetCellPosition(Short_t cellNumber)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }
  
  Int_t nabove, nbelow, middle;
  Short_t pos = -1;
  
  nabove = fNCells + 1;
  nbelow = 0;
  while (nabove - nbelow > 1) {
    middle = (nabove + nbelow) / 2;
    if (cellNumber == fCellNumber[middle-1]) {
      pos =   middle - 1;
      break;
    }
    if (cellNumber  < fCellNumber[middle-1]) nabove = middle;
    else                                     nbelow = middle;
  }
  
  return pos;
}

Int_t AliESDCaloCells::GetMCLabel(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells && fMCLabel) {
    return fMCLabel[pos];
  } else {
    return -1 ;
  }
}

Double_t AliESDCaloCells::GetEFraction(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells && fEFraction) {
    return fEFraction[pos];
  } else {
    return 0.;
  }
}

Int_t AliESDCaloCells::GetCellMCLabel(Short_t cellNumber)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }
  
  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && fCellNumber[pos] == cellNumber && fMCLabel) {
    return fMCLabel[pos];
  } else {
    return -1;
  }
}

Double_t AliESDCaloCells::GetCellEFraction(Short_t cellNumber)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }
  
  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && pos < fNCells && fCellNumber[pos] == cellNumber && fEFraction) {
    return fEFraction[pos];
  } else {
    return 0.;
  }
}

void AliESDCaloCells::SetEFraction(Short_t pos,  Double32_t efrac)
{
  // Sets the fraction of energy from MC with respect to data at the given position
  
  if (pos>=0 && pos < fNCells) 
  {
    if(!fEFraction) fEFraction = new Double32_t[fNCells];
    fEFraction[pos]  = efrac;
  } 
}

void AliESDCaloCells::SetCellEFraction(Short_t cellNumber, Double32_t efrac)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }
  
  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && pos < fNCells && fCellNumber[pos] == cellNumber) 
  {
    if(!fEFraction) fEFraction = new Double32_t[fNCells];
    fEFraction[pos] = efrac;
  } 
}


#endif
