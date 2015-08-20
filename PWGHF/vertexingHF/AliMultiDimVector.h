#ifndef ALIMULTIDIMVECTOR_H
#define ALIMULTIDIMVECTOR_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to store number of signal and background candidates     //
// in bins of cut variables                                      //
// Origin:       Elena Bruna (bruna@to.infn.it)                  //
// Updated:      Sergey Senyukov (senyukov@to.infn.it)           //
//               Francesco Prino (prino@to.infn.it)              //
// Last Updated: Giacomo Ortona (ortona@to.infn.it)              //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "TArrayF.h"
#include "TArrayI.h"
#include "TNamed.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"

class AliMultiDimVector :  public TNamed{

 public:
  AliMultiDimVector();
  AliMultiDimVector(const AliMultiDimVector &mv);
  AliMultiDimVector& operator=(const AliMultiDimVector &mv); 
  AliMultiDimVector(const char *name, const char *title, const Int_t nptbins, 
		    const Float_t* ptlimits, const Int_t npars, const Int_t *nofcells, const Float_t *loosecuts, const Float_t *tightcuts, const TString *axisTitles);
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
  Float_t   GetElement(const Int_t *ind, Int_t ptbin) const {
    ULong64_t elem=GetGlobalAddressFromIndices(ind,ptbin);
    return fVett[elem];
  }
  Float_t   GetPtLimit(Int_t i) const{return fPtLimits[i];}
  Int_t     GetPtBin(const Float_t pt) const{
    Int_t theBin=TMath::BinarySearch(fNPtBins+1,fPtLimits,pt);
    if(theBin>=fNPtBins) theBin=-1;
    return theBin;
  }
  void      GetEntireMultiDimVector(Float_t *vett) const {
    for(ULong64_t i=0; i<fNTotCells; i++) vett[i]=fVett[i];
  }

  Bool_t    GetIndicesFromGlobalAddress(ULong64_t globadd, Int_t *ind, Int_t &ptbin) const;
  ULong64_t GetGlobalAddressFromIndices(const Int_t *ind, Int_t ptbin) const;
  Bool_t    GetIndicesFromValues(const Float_t *values, Int_t *ind) const;
  ULong64_t GetGlobalAddressFromValues(const Float_t *values, Int_t ptbin) const;
  Bool_t    GetCutValuesFromGlobalAddress(ULong64_t globadd, Float_t *cuts, Int_t &ptbin) const;
  
  ULong64_t* GetGlobalAddressesAboveCuts(const Float_t *values, Float_t pt, Int_t& nVals) const{
    Int_t theBin=GetPtBin(pt);
    if(theBin>=0) return GetGlobalAddressesAboveCuts(values,theBin,nVals);
    else return 0x0;
  }
  ULong64_t* GetGlobalAddressesAboveCuts(const Float_t *values, Int_t ptbin, Int_t& nVals) const;
  Bool_t    GetGreaterThan(Int_t iVar) const {return fGreaterThan[iVar];}

  void SetElement(ULong64_t globadd,Float_t val) {fVett[globadd]=val;}
  void SetElement(Int_t *ind, Int_t ptbin, Float_t val){
    ULong64_t elem=GetGlobalAddressFromIndices(ind,ptbin);
    if(elem>fNTotCells){
      printf("SetElement: indices %d %d %d  ptbin %d elem %d\n",ind[0],ind[1],ind[2],ptbin,(Int_t)elem);
    }
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
  void Multiply(const AliMultiDimVector* mv,Float_t factor);
  void Multiply(const AliMultiDimVector* mv1, const AliMultiDimVector* mv2);
  void Add(const AliMultiDimVector* mv);
  void Sum(const AliMultiDimVector* mv1, const AliMultiDimVector* mv2);
  void LinearComb(const AliMultiDimVector* mv1, Float_t norm1, const AliMultiDimVector* mv2, Float_t norm2);
  void DivideBy(const AliMultiDimVector* mv);
  void Divide(const AliMultiDimVector* mv1, const AliMultiDimVector* mv2);
  void Sqrt();
  void Sqrt(const AliMultiDimVector* mv);
  
  void FindMaximum(Float_t& max_value, Int_t *ind, Int_t ptbin); 
  Int_t* FindLocalMaximum(Float_t& maxValue, Int_t *numFixed,Int_t* indFixed, Int_t nfixed,Int_t ptbin);

  TH2F*  Project(Int_t firstVar, Int_t secondVar, const Int_t* fixedVars, Int_t ptbin, Float_t norm=1.);

  void SuppressZeroBKGEffect(const AliMultiDimVector* BKG);
  AliMultiDimVector* ShrinkPtBins(Int_t firstBin, Int_t lastBin);

  void SetNewLimits(Float_t* loose,Float_t* tight);
  void SwapLimits(Int_t ilim);


  void PrintStatus();

 protected:
  void GetIntegrationLimits(Int_t iVar, Int_t iCell, Int_t& minbin, Int_t& maxbin) const;
  void GetFillRange(Int_t iVar, Int_t iCell, Int_t& minbin, Int_t& maxbin) const;
  Float_t   CountsAboveCell(ULong64_t globadd) const;

  //void SetMinLimits(Int_t nvar, Float_t* minlim);
  //void SetMaxLimits(Int_t nvar, Float_t* maxlim);
 private:
  static const Int_t fgkMaxNVariables=10;  // max. n. of selection variables
  static const Int_t fgkMaxNPtBins=10;     // max. n. of Pt bins

  Int_t     fNVariables;                   // n. of selection variables
  Int_t     fNPtBins;                      // n. of pt bins
  Float_t   fPtLimits[fgkMaxNPtBins+1];    // limits of pt bins
  Int_t     fNCutSteps[fgkMaxNVariables];  // n. of cut step for each variable
  Float_t   fMinLimits[fgkMaxNVariables];  // lower cut value for each variable
  Float_t   fMaxLimits[fgkMaxNVariables];  // higher cut value for each variable
  Bool_t    fGreaterThan[fgkMaxNVariables];// sign of the cut (> or <)
  TString   fAxisTitles[fgkMaxNVariables]; // titles for variables
  TArrayF   fVett;                   // array with n. of candidates vs. cuts
  ULong64_t fNTotCells;              // total number of matrix elements
  Bool_t    fIsIntegrated;           // flag for integrated matrix 

  ClassDef(AliMultiDimVector,2); // a multi-dimensional vector class

};

#endif
