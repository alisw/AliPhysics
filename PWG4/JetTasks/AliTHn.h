#ifndef AliTHn_H
#define AliTHn_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// optimized data container reusing functionality of AliCFContainer / THnsparse

#include "TObject.h"
#include "TString.h"
#include "AliCFContainer.h"

class TArrayF;
class TCollection;

class AliTHn : public AliCFContainer
{
 public:
  AliTHn();
  AliTHn(const Char_t* name, const Char_t* title,const Int_t nSelStep, const Int_t nVarIn, const Int_t* nBinIn);
  
  virtual ~AliTHn();
  
  virtual void  Fill(const Double_t *var, Int_t istep, Double_t weight=1.) ;
  virtual void  FillParent();
  virtual void  FillContainer(AliCFContainer* cont);
  
  TArrayF* GetValues(Int_t step) { return fValues[step]; }
  TArrayF* GetSumw2(Int_t step)  { return fSumw2[step]; }
  
  void DeleteContainers();
  void ReduceAxis();
  
  AliTHn(const AliTHn &c);
  AliTHn& operator=(const AliTHn& corr);
  virtual void Copy(TObject& c) const;

  virtual Long64_t Merge(TCollection* list);
  
protected:
  void Init();
  Long64_t GetGlobalBinIndex(const Int_t* binIdx);
  
  Long64_t fNBins;   // number of total bins
  Int_t    fNVars;   // number of variables
  Int_t    fNSteps;  // number of selection steps
  TArrayF **fValues; //[fNSteps] data container
  TArrayF **fSumw2;  //[fNSteps] data container
  TAxis** axisCache; //! cache axis pointers (about 50% of the time in Fill is spent in GetAxis otherwise)
  
  ClassDef(AliTHn, 3) // THn like container
};

#endif
