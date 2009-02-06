#ifndef ALICFVGRID_H
#define ALICFVGRID_H
//--------------------------------------------------------------------//
//                                                                    //
// AliCFGrid.cxx Class                                                //
// Just an interface to handle both AliCFGrid and AliCFGridSparse     //
// implementations                                                    //
//                                                                    //
//--------------------------------------------------------------------//

#include "AliCFFrame.h"

class TCollection;
class TH1D;
class TH2D;
class TH3D;

class AliCFVGrid : public AliCFFrame
{
 public:
  AliCFVGrid();
  AliCFVGrid(const AliCFVGrid & c);
  AliCFVGrid(const Char_t* name,const Char_t* title);
  AliCFVGrid(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t* nBinIn, const Double_t  *binLimitsIn=0);
  
  virtual ~AliCFVGrid();
  AliCFVGrid& operator=(const AliCFVGrid& corr);
 
  //abstract stuff

  virtual void  Fill(Double_t *var, Double_t weight=1.) = 0;

  virtual Float_t GetOverFlows(Int_t var) const = 0;
  virtual Float_t GetUnderFlows(Int_t var)const = 0;
  virtual Float_t GetEntries()const = 0 ;

  virtual Float_t GetElement(Int_t iel)const = 0; 
  virtual Float_t GetElement(Int_t *bin)const = 0; 
  virtual Float_t GetElement(Double_t *var)const = 0; 
  virtual Float_t GetElementError(Int_t iel)const = 0; 
  virtual Float_t GetElementError(Int_t *bin)const = 0; 
  virtual Float_t GetElementError(Double_t *var)const = 0; 

  virtual void SetElement(Int_t iel, Float_t val) = 0; 
  virtual void SetElement(Int_t *bin, Float_t val) = 0; 
  virtual void SetElement(Double_t *var, Float_t val) = 0; 
  virtual void SetElementError(Int_t iel, Float_t val) = 0; 
  virtual void SetElementError(Int_t *bin, Float_t val) = 0; 
  virtual void SetElementError(Double_t *var, Float_t val) = 0; 

  virtual TH1D*       Project( Int_t ivar) const = 0;
  virtual TH2D*       Project( Int_t ivar1, Int_t ivar2) const = 0;
  virtual TH3D*       Project( Int_t ivar1, Int_t ivar2,Int_t ivar3) const = 0;
  virtual AliCFVGrid* Project( Int_t nVars, Int_t* vars, Double_t* varMin, Double_t* varMax) const = 0;
  virtual TH1D* Slice(Int_t ivar, Double_t *varMin, Double_t *varMax) const = 0;
  virtual TH2D* Slice(Int_t ivar1, Int_t ivar2, Double_t *varMin, Double_t *varMax) const = 0;
  virtual TH3D* Slice(Int_t ivar1, Int_t ivar2, Int_t ivar3, Double_t *varMin, Double_t *varMax) const = 0;

  virtual void  UseAxisRange(Bool_t b) const = 0 ;

  //basic operations
  virtual void SumW2()=0;
  virtual void Add(AliCFVGrid* aGrid, Double_t c=1.) = 0;
  virtual void Add(AliCFVGrid* aGrid1 ,AliCFVGrid* aGrid2, Double_t c1=1.,Double_t c2=1.) = 0;
  virtual void Multiply(AliCFVGrid* aGrid, Double_t c=1.) = 0;
  virtual void Multiply(AliCFVGrid* aGrid1,AliCFVGrid* aGrid2, Double_t c1=1.,Double_t c2=1.) = 0;
  virtual void Divide(AliCFVGrid* aGrid, Double_t c=1.) = 0;
  virtual void Divide(AliCFVGrid* aGrid1, AliCFVGrid* aGrid2, Double_t c1=1., Double_t c2=1.,Option_t *option=0) = 0;
  virtual void Rebin(const Int_t* group) = 0;



  //implemented in AliCFVGrid

  virtual void Scale(Int_t iel, Double_t *fact); 
  virtual void Scale(Int_t* bin, Double_t *fact); 
  virtual void Scale(Double_t* var, Double_t *fact); 
  virtual void Scale(Double_t *fact); // To normalize MC to int lumi, for ex. 
  virtual Int_t    GetEmptyBins()const;
  virtual Int_t    CheckStats(Double_t thr) const;
  virtual Int_t    GetSumW2()const {return fSumW2;};
  virtual Int_t    GetEmptyBins(Double_t *varMin,Double_t *varMax) const;
  virtual Double_t GetIntegral() const;
  virtual Double_t GetIntegral(Int_t *binMin,Int_t *binMax) const;
  virtual Double_t GetIntegral(Double_t *varMin,Double_t *varMax) const;
  virtual Long64_t Merge(TCollection* list);
  virtual void Copy(TObject& c) const;

 protected:

  Double_t GetSum(Int_t ivar, Int_t *binMin, Int_t* binMax) const; 
  void GetScaledValues(Double_t *fact, Double_t *in, Double_t *out) const;
  //'hidden dimensions' when performing projections, 
  // default is kTRUE. please notice that 
  // if you you use AliCFGrid instead of AliCFGridSparse, 
  // only option kTRUE is actually available  
  Bool_t   fSumW2;//flag to check if calculation of squared weights enabled
   
  ClassDef(AliCFVGrid,3);
};
    
#endif

