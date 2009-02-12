#ifndef ALIMULTIDIMVECTOR_H
#define ALIMULTIDIMVECTOR_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to store number of signal and background candidates     //
// in bins of cut variables                                      //
// Origin:       Elena Bruna (bruna@to.infn.it)                  //
// Updated:      Sergey Senyukov (senyukov@to.infn.it)           //
// Last updated: Francesco Prino (prino@to.infn.it)              //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "TArrayF.h"
#include "TArrayI.h"
#include "TNamed.h"
#include "TH2.h"
#include "TString.h"

class AliMultiDimVector :  public TNamed{

 public:
  AliMultiDimVector();
  AliMultiDimVector(const AliMultiDimVector &mv);
  AliMultiDimVector(const char *name, const char *title, const Int_t npars, Int_t nptbins, Int_t *nofcells, Float_t *loosecuts, Float_t *tightcuts, TString *axisTitles);
  virtual ~AliMultiDimVector(){};

  ULong64_t GetNTotCells()            const {return fNTotCells;}
  Int_t     GetNVariables()           const {return fNVariables;}
  Int_t     GetNPtBins()              const {return fNPtBins;}
  Int_t     GetNCutSteps(Int_t iVar)  const {return fNCutSteps[iVar];}
  Float_t   GetMinLimit(Int_t iVar)   const {return fMinLimits[iVar];}
  Float_t   GetMaxLimit(Int_t iVar)   const {return fMaxLimits[iVar];}
  Float_t   GetCutStep(Int_t iVar)    const {return (fMaxLimits[iVar]-fMinLimits[iVar])/(Float_t)fNCutSteps[iVar];}
  TString   GetAxisTitle(Int_t iVar)  const {return fAxisTitles[iVar];}
  Bool_t    IsIntegrated()            const {return fIsIntegrated;}

  void CopyStructure(const AliMultiDimVector* mv);

  Float_t   GetCutValue(Int_t iVar, Int_t iCell) const{
    if(fGreaterThan[iVar]) return fMinLimits[iVar]+(Float_t)iCell*GetCutStep(iVar);
    else return fMaxLimits[iVar]-(Float_t)iCell*GetCutStep(iVar);
  }
  Float_t   GetElement(ULong64_t globadd) const {return fVett[globadd];}
  Float_t   GetElement(Int_t *ind, Int_t ptbin) const {
    ULong64_t elem=GetGlobalAddressFromIndices(ind,ptbin);
    return fVett[elem];
  }
  void      GetEntireMultiDimVector(Float_t *vett) const {
    for(ULong64_t i=0; i<fNTotCells; i++) vett[i]=fVett[i];
  }

  Bool_t    GetIndicesFromGlobalAddress(ULong64_t globadd, Int_t *ind, Int_t &ptbin) const;
  ULong64_t GetGlobalAddressFromIndices(Int_t *ind, Int_t ptbin) const;
  Bool_t    GetIndicesFromValues(Float_t *values, Int_t *ind) const;
  ULong64_t GetGlobalAddressFromValues(Float_t *values, Int_t ptbin) const;
  Bool_t    GetCutValuesFromGlobalAddress(ULong64_t globadd, Float_t *cuts, Int_t &ptbin) const;

  void SetElement(ULong64_t globadd,Float_t val) {fVett[globadd]=val;}
  void SetElement(Int_t *ind, Int_t ptbin, Float_t val){
    ULong64_t elem=GetGlobalAddressFromIndices(ind,ptbin);
    fVett[elem]=val;
  }
  void IncrementElement(Int_t *ind, Int_t ptbin){
    SetElement(ind,ptbin,GetElement(ind,ptbin)+1);
  }
  void IncrementElement(ULong64_t globadd){
    SetElement(globadd,GetElement(globadd)+1.);
  }

  void Fill(Float_t* values, Int_t ptbin);
  void FillAndIntegrate(Float_t* values, Int_t ptbin);
  void Integrate();

  void Reset(){
    for(ULong64_t i=0; i<fNTotCells; i++) fVett[i]=0.;
  }
  void MultiplyBy(Float_t factor);
  void Multiply(AliMultiDimVector* mv,Float_t factor);
  void Multiply(AliMultiDimVector* mv1, AliMultiDimVector* mv2);
  void Add(AliMultiDimVector* mv);
  void Sum(AliMultiDimVector* mv1, AliMultiDimVector* mv2);
  void LinearComb(AliMultiDimVector* mv1, Float_t norm1, AliMultiDimVector* mv2, Float_t norm2);
  void DivideBy(AliMultiDimVector* mv);
  void Divide(AliMultiDimVector* mv1, AliMultiDimVector* mv2);
  void Sqrt();
  void Sqrt(AliMultiDimVector* mv);
  
  void FindMaximum(Float_t& max_value, Int_t *ind, Int_t ptbin); 

  TH2F*  Project(Int_t firstVar, Int_t secondVar, Int_t* fixedVars, Int_t ptbin, Float_t norm=1.);

  void SuppressZeroBKGEffect(AliMultiDimVector* BKG);
   AliMultiDimVector* ShrinkPtBins(Int_t firstBin, Int_t lastBin);

 protected:
  void GetIntegrationLimits(Int_t iVar, Int_t iCell, Int_t& minbin, Int_t& maxbin) const;
  void GetFillRange(Int_t iVar, Int_t iCell, Int_t& minbin, Int_t& maxbin) const;
  Bool_t    GetGreaterThan(Int_t iVar) const {return fGreaterThan[iVar];}
  Float_t   CountsAboveCell(ULong64_t globadd) const;

 private:
  static const Int_t fgkMaxNVariables=10;  // max. n. of selection variables

  Int_t     fNVariables;                   // n. of selection variables
  Int_t     fNPtBins;                      // n. of pt bins
  Int_t     fNCutSteps[fgkMaxNVariables];  // n. of cut step for each variable
  Float_t   fMinLimits[fgkMaxNVariables];  // lower cut value for each variable
  Float_t   fMaxLimits[fgkMaxNVariables];  // higher cut value for each variable
  Bool_t    fGreaterThan[fgkMaxNVariables];// sign of the cut (> or <)
  TString   fAxisTitles[fgkMaxNVariables]; // titles for variables
  TArrayF   fVett;                   // array with n. of candidates vs. cuts
  ULong64_t fNTotCells;              // total number of matrix elements
  Bool_t    fIsIntegrated;           // flag for integrated matrix 

  ClassDef(AliMultiDimVector,1); // a multi-dimensional vector class

};

#endif
