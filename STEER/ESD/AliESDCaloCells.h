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
  virtual void Copy(TObject &obj) const;
  virtual AliVCaloCells * CopyCaloCells(Bool_t all) const;
  void Clear(const Option_t*);
  
  Bool_t IsEMCAL()  const { return (fType == kEMCALCell);}
  Bool_t IsPHOS()   const { return (fType == kPHOSCell) ;}
  Char_t  GetType() const { return  fType;}
  void    SetType(Char_t ttype) { fType = ttype; }

  void CreateContainer(Short_t nCells);
  void DeleteContainer();
  void Sort();
  
  Bool_t SetCell(Short_t pos, Short_t cellNumber, Double_t amplitude, Double_t time);
  
  Short_t GetNumberOfCells() const { return fNCells; }
  void SetNumberOfCells(Int_t n) { fNCells = n ; }
  inline Bool_t   GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude, Double_t &time) const;
  inline Double_t GetCellAmplitude(Short_t cellNumber);
  inline Double_t GetCellTime(Short_t cellNumber);
  inline Double_t GetAmplitude(Short_t pos) const;
  inline Double_t GetTime(Short_t pos) const;
  inline Short_t  GetCellNumber(Short_t pos) const;

 protected:
  Int_t       fNCells;       // Number of cells
  Short_t    *fCellNumber;   //[fNCells] array of cell numbers
  Double32_t *fAmplitude;    //[fNCells][0.,0.,16] array with cell amplitudes (= energy!)
  Double32_t *fTime;         //[fNCells][0.,0.,16] array with cell times
  Bool_t      fIsSorted;     //! true if cell arrays are sorted by index
  Char_t      fType;         // Cell type

  ClassDef(AliESDCaloCells, 2);
};


Bool_t AliESDCaloCells::GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude, Double_t & time) const 
{ 
  if (pos>=0 && pos<fNCells) {
    cellNumber = fCellNumber[pos];
    amplitude = fAmplitude[pos];
    time = fTime[pos];
    return kTRUE;
  } else {
    Error("GetCell","Invalid cell array index %d", pos);
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
    Warning("GetCellAmplitude","Wrong cell array index %d", pos);
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
    Warning("GetCellTime","Wrong cell array index %d", pos);
    return 0.;
  }
}

Double_t AliESDCaloCells::GetAmplitude(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    return fAmplitude[pos];
  } else {
    Error("GetAmplitude","Invalid cell array index %d", pos);
    return 0.;
  }
}

Double_t AliESDCaloCells::GetTime(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    return fTime[pos];
  } else {
    Error("GetTime","Invalid cell array index %d", pos);
    return 0.;
  }
}

Short_t AliESDCaloCells::GetCellNumber(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    return fCellNumber[pos];
  } else {
    Error("GetCellNumber","Invalid cell array index %d", pos);
    return fNCells;
  }
}

#endif
