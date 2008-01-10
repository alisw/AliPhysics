#ifndef ALICFCONTAINER_H
#define ALICFCONTAINER_H

/* $Id$ */

//--------------------------------------------------------------------//
//                                                                    //
// AliCFContainer Class                                               //
// Class to handle input data for correction Framework                // 
//                                                                    //
//--------------------------------------------------------------------//

#include "AliCFFrame.h"

class TH1F;
class TH2F;
class TH3F;
class AliCFGrid;

class AliCFContainer : public AliCFFrame
{
 public:
  AliCFContainer();
  AliCFContainer(const Char_t* name,const Char_t* title);
  AliCFContainer(const Char_t* name, const Char_t* title,const Int_t nSelStep, const Int_t nVarIn, const Int_t* nBinIn, const Float_t  *binLimitsIn=0);
  AliCFContainer(const AliCFContainer& c);
  
  virtual ~AliCFContainer();
  AliCFContainer& operator=(const AliCFContainer& corr);
  virtual Int_t GetNStep() const {return fNStep;};
  virtual void  SetBinLimits(Int_t varindex, Float_t * array);
  virtual void  Fill(Float_t *var, Int_t istep, Float_t weight=1.);

  virtual Float_t GetOverFlows(Int_t var,Int_t istep) const;
  virtual Float_t GetUnderFlows(Int_t var,Int_t istep)const ;
  virtual Float_t GetOverFlows(Int_t istep)const ;
  virtual Float_t GetUnderFlows(Int_t istep)const ;
  virtual Float_t GetEntries(Int_t istep)const ;
  virtual Int_t   GetEmptyBins(Int_t istep)const ;
  virtual Int_t   GetEmptyBins(Int_t istep, Float_t *varMin,Float_t *varMax) const ;
  virtual Float_t GetIntegral(Int_t istep) const ;
  virtual Float_t GetIntegral(Int_t istep, Float_t *varMin,Float_t *varMax) const ;
  //basic operations

  virtual void Copy(TObject& c) const;
  virtual void Add(AliCFContainer* aContainerToAdd, Float_t c=1.);
  virtual Long64_t Merge(TCollection* list);

  virtual TH1F* ShowProjection( Int_t ivar, Int_t istep) const;
  virtual TH2F* ShowProjection( Int_t ivar1, Int_t ivar2, Int_t istep) const;
  virtual TH3F* ShowProjection( Int_t ivar1, Int_t ivar2,Int_t ivar3, Int_t istep) const;
  virtual TH1F* ShowSlice( Int_t ivar, Float_t *varMin, Float_t *varMax, Int_t istep) const;
  virtual AliCFGrid * GetGrid(Int_t istep) const {return fGrid[istep];};
  
 private:
  Int_t    fNStep; //number of selection steps
  AliCFGrid **fGrid;//[fNStep]
  
  ClassDef(AliCFContainer,1);
};
    
#endif

