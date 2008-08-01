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

class TH1D;
class TH2D;
class TH3D;
class AliCFVGrid;
class TCollection;

class AliCFContainer : public AliCFFrame
{
 public:
  AliCFContainer();
  AliCFContainer(const Char_t* name,const Char_t* title);
  AliCFContainer(const Char_t* name, const Char_t* title,const Int_t nSelStep, const Int_t nVarIn, const Int_t* nBinIn, const Double_t  *binLimitsIn=0, const Bool_t useSparse = kTRUE);
  AliCFContainer(const AliCFContainer& c);
  
  virtual ~AliCFContainer();
  AliCFContainer& operator=(const AliCFContainer& corr);
  virtual Int_t GetNStep() const {return fNStep;};
  virtual void  SetBinLimits(Int_t varindex, Double_t * array);
  virtual void  Fill(Double_t *var, Int_t istep, Double_t weight=1.);

  virtual void   SetExcludeOffEntriesInProj(Bool_t in){fExclOffEntriesInProj=in;}; 
  virtual Bool_t GetExcludeOffEntriesInProj( ) const {return fExclOffEntriesInProj;}; 
  virtual Float_t GetOverFlows(Int_t var,Int_t istep) const;
  virtual Float_t GetUnderFlows(Int_t var,Int_t istep)const ;
  virtual Float_t GetEntries(Int_t istep)const ;
  virtual Int_t   GetEmptyBins(Int_t istep)const ;
  virtual Int_t   GetEmptyBins(Int_t istep, Double_t *varMin,Double_t *varMax) const ;
  virtual Double_t GetIntegral(Int_t istep) const ;
  virtual Double_t GetIntegral(Int_t istep, Double_t *varMin,Double_t *varMax) const ;
  //basic operations

  virtual void Copy(TObject& c) const;
  virtual void Add(AliCFContainer* aContainerToAdd, Double_t c=1.);
  virtual Long64_t Merge(TCollection* list);

  virtual TH1D* ShowProjection( Int_t ivar, Int_t istep) const;
  virtual TH2D* ShowProjection( Int_t ivar1, Int_t ivar2, Int_t istep) const;
  virtual TH3D* ShowProjection( Int_t ivar1, Int_t ivar2,Int_t ivar3, Int_t istep) const;
  virtual TH1D* ShowSlice( Int_t ivar, Double_t *varMin, Double_t *varMax, Int_t istep) const;
  virtual AliCFVGrid * GetGrid(Int_t istep) const {return (AliCFVGrid*)fGrid[istep];};
  
 private:
  Int_t    fNStep; //number of selection steps
  Bool_t fExclOffEntriesInProj; // exclude under/overflows in 
  AliCFVGrid **fGrid;//[fNStep]
  
  ClassDef(AliCFContainer,3);
};
    
#endif

