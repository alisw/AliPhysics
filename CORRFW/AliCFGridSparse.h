#ifndef ALICFGRIDSPARSE_H
#define ALICFGRIDSPARSE_H
//--------------------------------------------------------------------//
//                                                                    //
// AliCFGridSparse.cxx Class                                          //
// Class to handle N-dim maps for the correction Framework            // 
// uses a THnSparse to store the grid                                 //
// Author:S.Arcelli, silvia.arcelli@cern.ch
//--------------------------------------------------------------------//

#include "AliCFVGrid.h"
#include "THnSparse.h"
#include "AliLog.h"
class TH1D;
class TH2D;
class TH3D;

class AliCFGridSparse : public AliCFVGrid
{
 public:
  AliCFGridSparse();
  AliCFGridSparse(const Char_t* name,const Char_t* title);
  AliCFGridSparse(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t* nBinIn, const Double_t  *binLimitsIn=0);
  AliCFGridSparse(const AliCFGridSparse & c);
  
  virtual ~AliCFGridSparse();
  AliCFGridSparse& operator=(const AliCFGridSparse& corr) ;
  virtual void  SetBinLimits(Int_t ivar, Double_t * array);
  
  virtual void  Fill(Double_t *var, Double_t weight=1.);

  virtual Float_t GetOverFlows(Int_t var) const;
  virtual Float_t GetUnderFlows(Int_t var)const;
  virtual Float_t GetEntries()const;

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
  virtual void SetElementError(Int_t *bin, Float_t val) ; 
  virtual void SetElementError(Double_t *var, Float_t val); 

  virtual TH1D* Project( Int_t ivar) const;
  virtual TH2D* Project( Int_t ivar1, Int_t ivar2) const;
  virtual TH3D* Project( Int_t ivar1, Int_t ivar2,Int_t ivar3) const;
  virtual TH1D* Slice(Int_t ivar, Double_t* varMin, Double_t* varMax) const ; 
  virtual TH2D* Slice(Int_t ivar1, Int_t ivar2, Double_t *varMin, Double_t *varMax) const ;
  virtual TH3D* Slice(Int_t ivar1, Int_t ivar2, Int_t ivar3, Double_t *varMin, Double_t *varMax) const ;
  virtual void  SetRangeUser(Double_t* varMin, Double_t* varMax) ;


  //basic operations

  virtual void SumW2();
  virtual void Add(AliCFVGrid* aGrid, Double_t c=1.);
  virtual void Add(AliCFVGrid* aGrid1 ,AliCFVGrid* aGrid2, Double_t c1=1.,Double_t c2=1.);
  virtual void Multiply(AliCFVGrid* aGrid, Double_t c=1.);
  virtual void Multiply(AliCFVGrid* aGrid1,AliCFVGrid* aGrid2, Double_t c1=1.,Double_t c2=1.);
  virtual void Divide(AliCFVGrid* aGrid, Double_t c=1.);
  virtual void Divide(AliCFVGrid* aGrid1, AliCFVGrid* aGrid2, Double_t c1=1., Double_t c2=1.,Option_t *option=0);

  virtual void Rebin(const Int_t* group);

  THnSparse  *GetGrid() const {return fData;};//  Getter for the data Container: a THnSparse

  virtual void Copy(TObject& c) const;

  
 protected:

  THnSparse  *fData;//  The data Container: a THnSparse  
  ClassDef(AliCFGridSparse,2);
};
    
#endif

