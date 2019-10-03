#ifndef AliTHn_H
#define AliTHn_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// optimized data container reusing functionality of AliCFContainer / THnsparse
//
// Use AliTHn instead of AliCFContainer and your memory consumption will be drastically reduced
// As AliTHn derives from AliCFContainer, you can just replace your current AliCFContainer object by AliTHn
// Once you have the merged output, call FillParent() and you can use AliCFContainer as usual

#include "TObject.h"
#include "TString.h"
#include "AliCFContainer.h"

class TArray;
class TArrayF;
class TArrayD;
class TCollection;

class AliTHnBase : public AliCFContainer
{
public:
  AliTHnBase() : AliCFContainer() { }
  AliTHnBase(const Char_t* name, const Char_t* title,const Int_t nSelStep, const Int_t nVarIn, const Int_t* nBinIn) : AliCFContainer(name, title, nSelStep, nVarIn, nBinIn) { }
  
  virtual void Fill(const Double_t *var, Int_t istep, Double_t weight=1.) = 0;
  virtual void FillParent() = 0;
  virtual void FillContainer(AliCFContainer* cont) = 0;

  virtual TArray* GetValues(Int_t step) = 0;
  virtual TArray* GetSumw2(Int_t step) = 0;

  virtual void DeleteContainers() = 0;
  virtual void ReduceAxis() = 0;  
  
  ClassDef(AliTHnBase, 1) // AliTHn base class
};

template <class TemplateArray, typename TemplateType>
class AliTHnT : public AliTHnBase
{
 public:
  AliTHnT();
  AliTHnT(const Char_t* name, const Char_t* title,const Int_t nSelStep, const Int_t nVarIn, const Int_t* nBinIn);
  
  virtual ~AliTHnT();
  
  virtual void Fill(const Double_t *var, Int_t istep, Double_t weight=1.) ;
  virtual void FillParent();
  virtual void FillContainer(AliCFContainer* cont);
  
  virtual TArray* GetValues(Int_t step) { return fValues[step]; }
  virtual TArray* GetSumw2(Int_t step)  { return fSumw2[step]; }
  
  virtual void DeleteContainers();
  virtual void ReduceAxis();
  
  AliTHnT(const AliTHnT &c);
  AliTHnT& operator=(const AliTHnT& corr);
  virtual void Copy(TObject& c) const;

  virtual Long64_t Merge(TCollection* list);
  
protected:
  void Init();
  Long64_t GetGlobalBinIndex(const Int_t* binIdx);
  
  Long64_t fNBins;   // number of total bins
  Int_t    fNVars;   // number of variables
  Int_t    fNSteps;  // number of selection steps
  TemplateArray **fValues;  //[fNSteps] data container
  TemplateArray **fSumw2;   //[fNSteps] data container
  
  TAxis** axisCache; //! cache axis pointers (about 50% of the time in Fill is spent in GetAxis otherwise)
  Int_t* fNbinsCache; //! cache Nbins per axis
  Double_t* fLastVars; //! caching of last used bins (in many loops some vars are the same for a while)
  Int_t* fLastBins; //! caching of last used bins (in many loops some vars are the same for a while)
  
  ClassDef(AliTHnT, 5) // THn like container
};

typedef AliTHnT<TArrayF, Float_t> AliTHn;
typedef AliTHnT<TArrayD, Double_t> AliTHnD;

#endif
