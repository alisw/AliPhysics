/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD class to store calorimeter cell data
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#ifndef ALIAODCELLS_H
#define ALIAODCELLS_H

#include <TNamed.h>
#include <TMath.h>

class AliAODCaloCells : public TNamed 
{
 public:
  enum AODTwrs_t {kUndef = -1, 
		  kEMCAL, 
		  kPHOS};

  AliAODCaloCells();
  AliAODCaloCells(const char* name, const char* title, AODTwrs_t ttype=kUndef);
  
  virtual ~AliAODCaloCells();
  
  void CreateContainer(Short_t nCells);
  void DeleteContainer();
  void Sort();
  
  Bool_t SetCell(Short_t pos, Short_t cellNumber, Double_t amplitude);
  
  Short_t GetNumberOfCells() const { return fNCells; }
  inline Bool_t   GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude) const;
  inline Double_t GetCellAmplitude(Short_t cellNumber);
  inline Double_t GetAmplitude(Short_t pos) const;
  inline Short_t  GetCellNumber(Short_t pos) const;

  Char_t  GetType() const { return fType;}
  void    SetType(AODTwrs_t ttype) { fType=ttype; }

 protected:
  Int_t       fNCells;       // Number of cells
  Short_t    *fCellNumber;   //[fNCells] array of cell numbers
  Double32_t *fAmplitude;    //[fNCells][0.,600.,16] array with cell amplitudes (= energy!)
  Bool_t      fIsSorted;     //! true if cell arrays are sorted by index
  Char_t      fType;         // Cell type
  
 private:
  AliAODCaloCells(const AliAODCaloCells& tow); 
  AliAODCaloCells& operator=(const AliAODCaloCells& tow);
  
  ClassDef(AliAODCaloCells, 1);
};


Bool_t AliAODCaloCells::GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude) const 
{ 
  if (pos>=0 && pos<fNCells) {
    cellNumber = fCellNumber[pos];
    amplitude = fAmplitude[pos];
    return kTRUE;
  } else {
    Error("GetCell","Invalid cell array index %d", pos);
    return kFALSE;
  }
}


Double_t AliAODCaloCells::GetCellAmplitude(Short_t cellNumber)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }

  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && pos == cellNumber) {
    return fAmplitude[pos];
  } else {
    Error("GetCellAmplitude","Wrong cell array index %d", pos);
    return 0.;
  }
}


Double_t AliAODCaloCells::GetAmplitude(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    return fAmplitude[pos];
  } else {
    Error("GetAmplitude","Invalid cell array index %d", pos);
    return 0.;
  }
}


Short_t AliAODCaloCells::GetCellNumber(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    return fCellNumber[pos];
  } else {
    Error("GetCellNumber","Invalid cell array index %d", pos);
    return fNCells;
  }
}


#endif
