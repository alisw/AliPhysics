#ifndef ALICFGRIDSPARSE_H
#define ALICFGRIDSPARSE_H
//--------------------------------------------------------------------//
//                                                                    //
// AliCFGridSparse.cxx Class                                          //
// Class to handle N-dim maps for the correction Framework            // 
// uses a THnSparse to store the grid                                 //
// Author:S.Arcelli, silvia.arcelli@cern.ch
//--------------------------------------------------------------------//

#include "AliCFFrame.h"
#include "THnSparse.h"
#include "AliLog.h"
#include "TAxis.h"

class TH1D;
class TH2D;
class TH3D;

class AliCFGridSparse : public AliCFFrame
{
 public:
  AliCFGridSparse();
  AliCFGridSparse(const Char_t* name, const Char_t* title);
  AliCFGridSparse(const Char_t* name, const Char_t* title, Int_t nVarIn, const Int_t* nBinIn);
  AliCFGridSparse(const AliCFGridSparse& c);
  virtual ~AliCFGridSparse();
  AliCFGridSparse& operator=(const AliCFGridSparse& c);
  virtual void Copy(TObject& c) const;

  // AliCFFrame functions
  virtual Int_t      GetNVar() const {return fData->GetNdimensions();}
  virtual void       PrintBinLimits() const ;
  virtual void       PrintNBins() const ; 
  virtual void       SetBinLimits(Int_t ivar, Double_t min, Double_t max); // for uniform bin width only
  virtual void       SetBinLimits(Int_t ivar, const Double_t * array);     // for variable or uniform bin width
  virtual void       GetBinLimits(Int_t ivar, Double_t * array) const ;
  virtual Double_t * GetBinLimits(Int_t ivar) const ;
  virtual Long_t     GetNBinsTotal() const ;
  virtual Long_t     GetNFilledBins() const {return fData->GetNbins();}
  virtual Int_t      GetNBins(Int_t ivar) const {return fData->GetAxis(ivar)->GetNbins();}
  virtual Int_t *    GetNBins() const ;
  virtual Float_t    GetBinCenter(Int_t ivar,Int_t ibin) const ;
  virtual Float_t    GetBinSize  (Int_t ivar,Int_t ibin) const ;
  //virtual void       GetBinCenters(const Int_t *ibin, Float_t *binCenter) const ;
  //virtual void       GetBinSizes  (const Int_t *ibin, Float_t *binSizes)  const ;
  virtual TAxis    * GetAxis(Int_t ivar) const {return fData->GetAxis(ivar);}

  virtual void          SetVarTitle(Int_t ivar, const Char_t* lab) {fData->GetAxis(ivar)->SetTitle(lab);}
  virtual const Char_t* GetVarTitle(Int_t ivar) const {return GetAxis(ivar)->GetTitle();}
  virtual Int_t         GetVar(const Char_t* title) const ; // returns the variable corresponding to the given title

  // probably not needed anymore
  //virtual Int_t      GetBinIndex(const Int_t *ibin) const ;
  //virtual void       GetBinIndex(Int_t iel, const Int_t *ibin) const ;
  //virtual Int_t      GetBinIndex(Int_t ivar, Int_t ind) const ;

  virtual void    Fill(const Double_t *var, Double_t weight=1.);
  virtual Float_t GetEntries()const;
  virtual Float_t GetElement(Long_t iel)               const; 
  virtual Float_t GetElement(const Int_t *bin)         const; 
  virtual Float_t GetElement(const Double_t *var)      const; 
  virtual Float_t GetElementError(Long_t iel)          const;  
  virtual Float_t GetElementError(const Int_t *bin)    const; 
  virtual Float_t GetElementError(const Double_t *var) const; 
  virtual void    SetElement(Long_t iel, Float_t val); 
  virtual void    SetElement(const Int_t *bin, Float_t val); 
  virtual void    SetElement(const Double_t *var, Float_t val); 
  virtual void    SetElementError(Long_t iel, Float_t val); 
  virtual void    SetElementError(const Int_t *bin, Float_t val) ; 
  virtual void    SetElementError(const Double_t *var, Float_t val); 

  virtual TH1*             Project(Int_t ivar1, Int_t ivar2=-1, Int_t ivar3=-1) const {return Slice(ivar1,ivar2,ivar3,0x0,0x0,kFALSE);}
  virtual TH1*             Slice(Int_t ivar1, Int_t ivar2=-1, Int_t ivar3=-1, 
				 const Double_t *varMin=0x0, const Double_t *varMax=0x0, Bool_t useBins=0) const ; 
  virtual AliCFGridSparse* MakeSlice(Int_t nVars, const Int_t* vars,
				   const Double_t* varMin, const Double_t* varMax, Bool_t useBins=0) const ;

  virtual void             SetRangeUser(Int_t iVar, Double_t varMin, Double_t varMax, Bool_t useBins=kFALSE) const ;
  virtual void             SetRangeUser(const Double_t* varMin, const Double_t* varMax, Bool_t useBins=kFALSE) const ;
  virtual void             Smooth() ;

  //basic operations
  virtual void     SumW2();
  virtual void     Add(const AliCFGridSparse* aGrid, Double_t c=1.);
  virtual void     Add(const AliCFGridSparse* aGrid1 ,const AliCFGridSparse* aGrid2, Double_t c1=1.,Double_t c2=1.);
  virtual void     Multiply(const AliCFGridSparse* aGrid, Double_t c=1.);
  virtual void     Multiply(const AliCFGridSparse* aGrid1,const AliCFGridSparse* aGrid2, Double_t c1=1.,Double_t c2=1.);
  virtual void     Divide(const AliCFGridSparse* aGrid, Double_t c=1.);
  virtual void     Divide(const AliCFGridSparse* aGrid1, const AliCFGridSparse* aGrid2, Double_t c1=1., Double_t c2=1.,Option_t *option=0);
  virtual void     Rebin(const Int_t* group);
  virtual void     Scale(Long_t iel, const Double_t *fact); 
  virtual void     Scale(const Int_t* bin, const Double_t *fact); 
  virtual void     Scale(const Double_t* var, const Double_t *fact); 
  virtual void     Scale(const Double_t *fact); // To normalize MC to int lumi, for ex. 
  virtual Int_t    CheckStats(Double_t thr) const;
  virtual Int_t    GetSumW2() const {return fSumW2;};
  virtual Double_t GetIntegral() const;
  //virtual Double_t GetIntegral(const Double_t *varMin, const Double_t *varMax) const;
  virtual Long64_t Merge(TCollection* list);

  virtual void     SetGrid(THnSparse* grid) {if (fData) delete fData ; fData=grid;}
  THnSparse   *    GetGrid() const {return fData;}

  virtual Float_t GetOverFlows (Int_t var, Bool_t excl=kFALSE) const;
  virtual Float_t GetUnderFlows(Int_t var, Bool_t excl=kFALSE) const;
  virtual Long_t  GetEmptyBins() const;

  /*  FUNCTIONS TO REMOVE   */
  virtual AliCFGridSparse* Project(Int_t nVars, const Int_t* vars, const Double_t* varMin, const Double_t* varMax, Bool_t useBins=0) const 
  {return MakeSlice(nVars,vars,varMin,varMax,useBins);}


 protected:

  //protected functions
  void     GetScaledValues(const Double_t *fact, const Double_t *in, Double_t *out) const;
  void     SetAxisRange(TAxis* axis, Double_t min, Double_t max, Bool_t useBins) const;
  void     GetProjectionName (TString& s,Int_t var0, Int_t var1=-1, Int_t var2=-1) const;
  void     GetProjectionTitle(TString& s,Int_t var0, Int_t var1=-1, Int_t var2=-1) const;

  // data members:
  Bool_t      fSumW2    ; // Flag to check if calculation of squared weights enabled
  THnSparse  *fData     ; // The data Container: a THnSparse  

  ClassDef(AliCFGridSparse,3);
};


//inline functions :

inline Long_t AliCFGridSparse::GetNBinsTotal() const {
  Long_t n=1;
  for (Int_t iVar=0; iVar<GetNVar(); iVar++) {
    n *= fData->GetAxis(iVar)->GetNbins();
  }
  return n ;
}

inline void AliCFGridSparse::PrintNBins() const
{
  //
  // printing the array containing the # of bins  
  //
  for (Int_t i=0;i<GetNVar();i++) {
    AliInfo(Form("bins in axis %i are: %i",i,fData->GetAxis(i)->GetNbins()));
  }
} 

inline void AliCFGridSparse::PrintBinLimits() const
{
  //
  // printing the bin limits for each variable  
  //
  for (Int_t iVar=0; iVar<GetNVar(); iVar++) {
    AliInfo(Form("variable %d :",iVar));
    const Int_t nBins = GetNBins(iVar) ;
    Double_t *array = new Double_t[nBins+1];
    GetBinLimits(iVar,array);
    for (Int_t iBin=0; iBin<nBins; iBin++) {
      AliInfo(Form("    bin limit index %i is: %e",iBin,array[iBin]));
    }
    delete [] array ;
  }
}

inline void AliCFGridSparse::GetBinLimits(Int_t ivar, Double_t * array) const {
  TAxis * axis = fData->GetAxis(ivar) ;
  Int_t nBins = axis->GetNbins();
  for (Int_t iBin=0; iBin<nBins; iBin++) array[iBin] = axis->GetBinLowEdge(iBin+1);
  array[nBins] = axis->GetBinUpEdge(nBins);
}

inline Int_t* AliCFGridSparse::GetNBins() const {
  Int_t *bins = new Int_t[GetNVar()];
  for (Int_t iVar=0; iVar<GetNVar(); iVar++) {
    bins[iVar] = GetNBins(iVar) ;
  }
  return bins;
}

inline Double_t* AliCFGridSparse::GetBinLimits(Int_t ivar) const {
  Double_t * binLimits = new Double_t[GetNBins(ivar)+1] ;
  GetBinLimits(ivar,binLimits);
  return binLimits;
}

inline Int_t AliCFGridSparse::GetVar(const Char_t* title) const {
  TString str(title);
  for (Int_t iVar=0; iVar<GetNVar(); iVar++) {
    if (!str.CompareTo(GetVarTitle(iVar))) return iVar;
  }
  AliError("Variable not found");
  return -1;
}

inline void AliCFGridSparse::GetProjectionName (TString& s, Int_t var0, Int_t var1, Int_t var2) const {
  s.Form("%s_proj-%s",GetName(),GetVarTitle(var0));
  if (var1>=0) {
    s.Append(Form("-%s",GetVarTitle(var1)));
    if (var2>=0) s.Append(Form("-%s",GetVarTitle(var2)));
  }
}

inline void AliCFGridSparse::GetProjectionTitle(TString& s, Int_t var0, Int_t var1, Int_t var2) const {
  s.Form("%s: projection on %s",GetTitle(),GetVarTitle(var0));
  if (var1>=0) {
    s.Append(Form("-%s",GetVarTitle(var1)));
    if (var2>=0) s.Append(Form("-%s",GetVarTitle(var2)));
  }
}

#endif

