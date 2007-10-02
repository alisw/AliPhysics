/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD class to store tower data
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#ifndef ALIAODTOWERS_H
#define ALIAODTOWERS_H

#include <TNamed.h>
#include <TMath.h>

class AliAODTowers : public TNamed 
{
 public:
  enum AODTwrs_t {kUndef = -1, 
		  kEMCAL, 
		  kPHOS};

  AliAODTowers();
  AliAODTowers(const char* name, const char* title, AODTwrs_t ttype=kUndef);
  
  virtual ~AliAODTowers();
  
  void CreateContainer(Short_t nTowers);
  void DeleteContainer();
  
  Bool_t SetTower(Short_t pos, Short_t index, Double_t amplitude);
  
  Short_t GetNumberOfTowers() const { return fNTowers; }
  inline Bool_t   GetTower(Short_t pos, Short_t &index, Double_t &amplitude) const;
  inline Double_t GetTowerAmplitude(Short_t towerIndex);
  inline Double_t GetAmplitude(Short_t pos) const;
  inline Short_t  GetIndex(Short_t pos) const;

  Char_t  GetType() const { return fType;}
  void    SetType(AODTwrs_t ttype) { fType=ttype; }

 protected:
  Short_t     fNTowers;      // Number of towers
  Short_t    *fIndex;        //[fNTowers] array of tower indices
  Double32_t *fAmplitude;    //[fNTowers][0.,600.,16] array with tower amplitudes (= energy!)
  Bool_t      fIsSorted;     //! true if tower arrays are sorted by index
  Char_t      fType;         // Towers type
  
 private:
  AliAODTowers(const AliAODTowers& tow); 
  AliAODTowers& operator=(const AliAODTowers& tow);
  
  ClassDef(AliAODTowers, 1);
};


Bool_t AliAODTowers::GetTower(Short_t pos, Short_t &index, Double_t &amplitude) const 
{ 
  if (pos>=0 && pos<fNTowers) {
    index = fIndex[pos];
    amplitude = fAmplitude[pos];
    return kTRUE;
  } else {
    Error("GetTower","Invalid tower index %d", pos);
    return kFALSE;
  }
}


Double_t AliAODTowers::GetTowerAmplitude(Short_t towerIndex)
{ 
  if (!fIsSorted) {
    // sort it here!
    Warning("GetTowerAmplitude","Array sort has to be implemented");
    fIsSorted = kTRUE;
  }

Short_t pos = TMath::BinarySearch(fNTowers, fIndex, towerIndex);
  if (pos>=0 && pos == towerIndex) {
    return fAmplitude[pos];
  } else {
    Error("GetTowerAmplitude","Wrong tower index %d",pos);
    return 0.;
  }
}


Double_t AliAODTowers::GetAmplitude(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNTowers) {
    return fAmplitude[pos];
  } else {
    Error("GetAmplitude","Invalid tower index %d",pos);
    return 0.;
  }
}


Short_t AliAODTowers::GetIndex(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNTowers) {
    return fIndex[pos];
  } else {
    Error("GetIndex","Invalid tower index %d",pos);
    return fNTowers;
  }
}


#endif
