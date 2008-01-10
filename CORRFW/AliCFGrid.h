#ifndef ALICFGRID_H
#define ALICFGRID_H

/* $Id$ */

//--------------------------------------------------------------------//
//                                                                    //
// AliCFGrid.cxx Class                                                //
// Class to handle N-dim maps for the correction Framework            // 
//                                                                    //
//--------------------------------------------------------------------//

#include "AliCFFrame.h"

class TH1F;
class TH2F;
class TH3F;

class AliCFGrid : public AliCFFrame
{
 public:
  AliCFGrid();
  AliCFGrid(const Char_t* name,const Char_t* title);
  AliCFGrid(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t* nBinIn, const Float_t  *binLimitsIn=0);
  AliCFGrid(const AliCFGrid & c);
  
  virtual ~AliCFGrid();
  AliCFGrid& operator=(const AliCFGrid& corr);
  virtual void  Fill(Float_t *var, Float_t weight=1.);
  virtual Float_t GetOverFlows(Int_t var) const;
  virtual Float_t GetUnderFlows(Int_t var)const ;
  virtual Float_t GetOverFlows()const ;
  virtual Float_t GetUnderFlows()const ;
  virtual Float_t GetEntries()const ;
  virtual Int_t GetEmptyBins()const ;
  virtual Int_t CheckEfficiencyStats(Float_t thr) const;
  virtual Int_t GetSumW2()const {return fSumW2;};
  virtual Float_t GetElement(Int_t iel)const; 
  virtual Float_t GetElement(Int_t *bin)const; 
  virtual Float_t GetElement(Float_t *var)const; 
  virtual Float_t GetElementError(Int_t iel)const; 
  virtual Float_t GetElementError(Int_t *bin)const; 
  virtual Float_t GetElementError(Float_t *var)const; 
  virtual void SetElement(Int_t iel, Float_t val); 
  virtual void SetElement(Int_t *bin, Float_t val); 
  virtual void SetElement(Float_t *var, Float_t val); 
  virtual void SetElementError(Int_t iel, Float_t val); 
  virtual void SetElementError(Int_t *bin, Float_t val); 
  virtual void SetElementError(Float_t *var, Float_t val); 
  virtual void Scale(Int_t iel, Float_t *fact); 
  virtual void Scale(Int_t* bin, Float_t *fact); 
  virtual void Scale(Float_t* var, Float_t *fact); 
  virtual void Scale(Float_t *fact); // To normalize MC to int lumi, for ex. 
  virtual Int_t GetEmptyBins(Float_t *varMin,Float_t *varMax) const ;
  virtual Float_t GetIntegral() const ;
  virtual Float_t GetIntegral(Int_t *binMin,Int_t *binMax) const ;
  virtual Float_t GetIntegral(Float_t *varMin,Float_t *varMax) const ;
  virtual TH1F* Project( Int_t ivar) const;
  virtual TH2F* Project( Int_t ivar1, Int_t ivar2) const;
  virtual TH3F* Project( Int_t ivar1, Int_t ivar2,Int_t ivar3) const;
  virtual TH1F* Slice( Int_t ivar, Float_t *varMin, Float_t *varMax) const;
  virtual Float_t GetSum(Int_t ivar, Int_t *binMin, Int_t* binMax) const; 

  //basic operations

  virtual void SumW2();
  virtual void Copy(TObject& c) const;
  virtual void Add(AliCFGrid* aGrid, Float_t c=1.);
  virtual void Add(AliCFGrid* aGrid1 ,AliCFGrid* aGrid2, Float_t c1=1.,Float_t c2=1.);
  virtual void Multiply(AliCFGrid* aGrid, Float_t c=1.);
  virtual void Multiply(AliCFGrid* aGrid1,AliCFGrid* aGrid2, Float_t c1=1.,Float_t c2=1.);
  virtual void Divide(AliCFGrid* aGrid, Float_t c=1.,Option_t *option=0);
  virtual void Divide(AliCFGrid* aGrid1, AliCFGrid* aGrid2, Float_t c1=1., Float_t c2=1.,Option_t *option=0);
  virtual Long64_t Merge(TCollection* list);

  
 protected:
  Bool_t   fSumW2;//flag to check if calculation of squared weights enabled
  Float_t  fNunflTot;//Total number of underflows
  Float_t  fNovflTot;//Total number of underflows
  Float_t  fNentriesTot;//Total number of entries 
  Float_t  *fNunfl;//[fNVar] underflows in each dimension
  Float_t  *fNovfl;//[fNVar] overflows in each dimension

  Float_t  *fData;//[fNDim] The data Container
  Float_t  *fErr2;//[fNDim] The squared weights Container (only with SumW2)

  
  ClassDef(AliCFGrid,1);
};
    
#endif

