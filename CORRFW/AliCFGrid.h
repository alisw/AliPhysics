#ifndef ALICFGRID_H
#define ALICFGRID_H

//--------------------------------------------------------------------//
//                                                                    //
// AliCFGrid.cxx Class                                                //
// Class to handle N-dim maps for the correction Framework            // 
// The class uses a une-dimensional array of floats to store the grid //     
//                                                                    //
//--------------------------------------------------------------------//

#include "AliCFVGrid.h"
#include "AliLog.h"

class TH1D;
class TH2D;
class TH3D;

class AliCFGrid : public AliCFVGrid
{
 public:
  AliCFGrid();
  AliCFGrid(const Char_t* name,const Char_t* title);
  AliCFGrid(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t* nBinIn, const Double_t  *binLimitsIn=0);
  AliCFGrid(const AliCFGrid & c);
  
  virtual ~AliCFGrid();
  AliCFGrid& operator=(const AliCFGrid& corr);

  virtual void  Fill(Double_t *var, Double_t weight=1.);

  virtual void   SetExcludeOffEntriesInProj(Bool_t in); 
  virtual Bool_t GetExcludeOffEntriesInProj( ) const; 
  virtual Float_t GetOverFlows(Int_t var) const;
  virtual Float_t GetUnderFlows(Int_t var)const ;
  virtual Float_t GetEntries()const ;

  virtual Float_t GetElement(Int_t iel)const; 
  virtual Float_t GetElement(Int_t *bin)const; 
  virtual Float_t GetElement(Double_t *var)const; 
  virtual Float_t GetElementError(Int_t iel)const; 
  virtual Float_t GetElementError(Int_t *bin)const; 
  virtual Float_t GetElementError(Double_t *var)const; 
  virtual void SetElement(Int_t iel, Float_t val); 
  virtual void SetElement(Int_t *bin, Float_t val); 
  virtual void SetElement(Double_t *var, Float_t val); 
  virtual void SetElementError(Int_t iel, Float_t val); 
  virtual void SetElementError(Int_t *bin, Float_t val); 
  virtual void SetElementError(Double_t *var, Float_t val); 

  virtual TH1D* Project( Int_t ivar) const;
  virtual TH2D* Project( Int_t ivar1, Int_t ivar2) const;
  virtual TH3D* Project( Int_t ivar1, Int_t ivar2,Int_t ivar3) const;
  virtual TH1D* Slice( Int_t ivar, Double_t *varMin, Double_t *varMax) const;
  virtual TH2D* Slice( Int_t /*ivar1*/, Int_t /*ivar2*/, Double_t */*varMin*/, Double_t */*varMax*/) const 
  {AliWarning("not implemented"); return 0x0;}
  virtual TH3D* Slice( Int_t /*ivar1*/, Int_t /*ivar2*/, Int_t /*ivar3*/, Double_t */*varMin*/, Double_t */*varMax*/) const 
  {AliWarning("not implemented"); return 0x0;}

  //basic operations

  virtual void SumW2();
  virtual void Copy(TObject& c) const;
  virtual void Add(AliCFVGrid* aGrid, Double_t c=1.);
  virtual void Add(AliCFVGrid* aGrid1 ,AliCFVGrid* aGrid2, Double_t c1=1.,Double_t c2=1.);
  virtual void Multiply(AliCFVGrid* aGrid, Double_t c=1.);
  virtual void Multiply(AliCFVGrid* aGrid1,AliCFVGrid* aGrid2, Double_t c1=1.,Double_t c2=1.);
  virtual void Divide(AliCFVGrid* aGrid, Double_t c=1.);
  virtual void Divide(AliCFVGrid* aGrid1, AliCFVGrid* aGrid2, Double_t c1=1., Double_t c2=1.,Option_t *option=0);

  void Rebin(const Int_t* group);

  
 protected:

  Float_t  fNentriesTot;//Total number of entries 
  Float_t  *fNunfl;//[fNVar] underflows in each dimension
  Float_t  *fNovfl;//[fNVar] overflows in each dimension
  Float_t  *fData;//[fNDim] The data Container
  Float_t  *fErr2;//[fNDim] The squared weights Container (only with SumW2)

  
  ClassDef(AliCFGrid,4);
};
    
#endif

