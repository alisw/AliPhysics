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
#include "AliCFGridSparse.h"

class TH1D;
class TH2D;
class TH3D;
class TCollection;

class AliCFContainer : public AliCFFrame
{
 public:
  AliCFContainer();
  AliCFContainer(const Char_t* name, const Char_t* title,const Int_t nSelStep, const Int_t nVarIn, const Int_t* nBinIn);
  AliCFContainer(const AliCFContainer& c);
  AliCFContainer& operator=(const AliCFContainer& corr);
  virtual void Copy(TObject& c) const;

  virtual ~AliCFContainer();

  // AliCFFrame functions
  virtual Int_t      GetNVar()                                       const {return fGrid[0]->GetNVar();}
  virtual void       PrintBinLimits()                                const {fGrid[0]->PrintBinLimits();}
  virtual void       PrintNBins()                                    const {fGrid[0]->PrintNBins();}
  virtual void       SetBinLimits(Int_t ivar, Double_t min, Double_t max) ; // for uniform bin width only
  virtual void       SetBinLimits(Int_t ivar, const Double_t * array) ;     // for variable or uniform bin width
  virtual void       SetBinLabel (Int_t ivar,Int_t ibin, const Char_t* label);// set a label to bin ibin on axis ivar
  virtual void       SetBinContent(Int_t* bin, Int_t step, Double_t value);
  virtual void       SetBinError  (Int_t* bin, Int_t step, Double_t value);
  virtual void       GetBinLimits(Int_t ivar, Double_t * array)      const {return fGrid[0]->GetBinLimits(ivar,array);}
  virtual Double_t * GetBinLimits(Int_t ivar)                        const {return fGrid[0]->GetBinLimits(ivar);}
  virtual Long_t     GetNBinsTotal()                                 const {return fGrid[0]->GetNBinsTotal()*fNStep;}
  virtual Int_t      GetNBins(Int_t ivar)                            const {return fGrid[0]->GetNBins(ivar);}
  virtual Int_t    * GetNBins()                                      const {return fGrid[0]->GetNBins();}
  virtual Float_t    GetBinCenter(Int_t ivar,Int_t ibin)             const {return fGrid[0]->GetBinCenter(ivar,ibin);}
  virtual Float_t    GetBinSize  (Int_t ivar,Int_t ibin)             const {return fGrid[0]->GetBinSize  (ivar,ibin);}
  virtual Float_t    GetBinContent(const Int_t* coordinates, Int_t step) const {return fGrid[step]->GetGrid()->GetBinContent(coordinates);}
  virtual Float_t    GetBinError  (const Int_t* coordinates, Int_t step) const {return fGrid[step]->GetGrid()->GetBinError  (coordinates);}
  virtual const Char_t* GetBinLabel (Int_t ivar,Int_t ibin)          const {return GetAxis(ivar,0)->GetBinLabel(ibin);}

  virtual void       Print(const Option_t*) const ;

  virtual TAxis       * GetAxis(Int_t ivar, Int_t istep) const {return fGrid[istep]->GetAxis(ivar);}
  virtual void          SetVarTitle (Int_t ivar,  const Char_t* title) ;
  virtual void          SetStepTitle(Int_t istep, const Char_t* title) ;
  virtual const Char_t* GetVarTitle (Int_t ivar)  const {return GetAxis(ivar,0)->GetTitle();}
  virtual const Char_t* GetStepTitle(Int_t istep) const {return fGrid[istep]->GetTitle();}
  virtual Int_t         GetStep(const Char_t* title) const ; // returns the step     corresponding to the given title
  virtual Int_t         GetVar (const Char_t* title) const ; // returns the variable corresponding to the given title

  virtual Int_t GetNStep() const {return fNStep;};
  virtual void  SetNStep(Int_t nStep) {fNStep=nStep;}
  virtual void  Fill(const Double_t *var, Int_t istep, Double_t weight=1.) ;

  virtual Float_t  GetOverFlows (Int_t var,Int_t istep,Bool_t excl=kFALSE) const;
  virtual Float_t  GetUnderFlows(Int_t var,Int_t istep,Bool_t excl=kFALSE) const ;
  virtual Float_t  GetEntries  (Int_t istep) const ;
  virtual Long_t   GetEmptyBins(Int_t istep) const {return fGrid[istep]->GetEmptyBins();}
  virtual Double_t GetIntegral (Int_t istep) const ;

  //basic operations
  virtual void     Add(const AliCFContainer* aContainerToAdd, Double_t c=1.);
  virtual Long64_t Merge(TCollection* list);

  virtual TH1* Project (Int_t istep, Int_t ivar1, Int_t ivar2=-1 ,Int_t ivar3=-1) const;
  virtual AliCFContainer* MakeSlice(Int_t nVars, const Int_t* vars, const Double_t* varMin=0x0, const Double_t* varMax=0x0, Bool_t useBins=0) const ;
  virtual AliCFContainer* MakeSlice(Int_t nStep, const Int_t* steps, 
				    Int_t nVars, const Int_t* vars, const Double_t* varMin=0x0, const Double_t* varMax=0x0, 
				    Bool_t useBins=0) const ;
  virtual void  Smooth(Int_t istep) {GetGrid(istep)->Smooth();}

  virtual void  SetRangeUser(Int_t ivar, Double_t varMin, Double_t varMax, Bool_t useBins=kFALSE) const ;
  virtual void  SetRangeUser(const Double_t* varMin, const Double_t* varMax, Bool_t useBins=kFALSE) const ;

  virtual void  SetGrid(Int_t step, AliCFGridSparse* grid) {if (fGrid[step]) delete fGrid[step]; fGrid[step]=grid;}
  virtual AliCFGridSparse * GetGrid(Int_t istep) const {return fGrid[istep];};

  virtual void  Scale(Double_t factor) const;

  /****   TO BE REMOVED SOON ******/
  virtual TH1D* ShowProjection( Int_t ivar,  Int_t istep)                          const {return (TH1D*)Project(istep,ivar);}
  virtual TH2D* ShowProjection( Int_t ivar1, Int_t ivar2, Int_t istep)             const {return (TH2D*)Project(istep,ivar1,ivar2);}
  virtual TH3D* ShowProjection( Int_t ivar1, Int_t ivar2,Int_t ivar3, Int_t istep) const {return (TH3D*)Project(istep,ivar1,ivar2,ivar3);}
  
 private:
  Int_t    fNStep; //number of selection steps
  AliCFGridSparse **fGrid;//[fNStep]
  
  ClassDef(AliCFContainer,5);
};

inline void AliCFContainer::SetBinLimits(Int_t ivar, const Double_t* array) {
  for (Int_t iStep=0; iStep<GetNStep(); iStep++) {
    fGrid[iStep]->SetBinLimits(ivar,array);
  }
}

inline void AliCFContainer::SetBinLimits(Int_t ivar, Double_t min, Double_t max) {
  for (Int_t iStep=0; iStep<GetNStep(); iStep++) {
    fGrid[iStep]->SetBinLimits(ivar,min,max);
  }
}

inline void AliCFContainer::SetVarTitle(Int_t ivar, const Char_t* title) {
  for (Int_t iStep=0; iStep<fNStep; iStep++) {
    GetAxis(ivar,iStep)->SetTitle(title);
  }
}

inline void AliCFContainer::SetStepTitle(Int_t istep, const Char_t* title) {
  fGrid[istep]->SetTitle(title);
}

inline Int_t AliCFContainer::GetStep(const Char_t* title) const {
  TString str(title);
  for (Int_t iStep=0; iStep<fNStep; iStep++) {
    if (!str.CompareTo(GetStepTitle(iStep))) return iStep;
  }
  AliError("Step not found");
  return -1;
}

inline Int_t AliCFContainer::GetVar(const Char_t* title) const {
  return fGrid[0]->GetVar(title);
}

inline void AliCFContainer::SetBinLabel(Int_t iVar, Int_t iBin, const Char_t* label) {
  for (Int_t iStep=0; iStep<GetNStep(); iStep++) GetAxis(iVar,iStep)->SetBinLabel(iBin,label);
}

inline void  AliCFContainer::Scale(Double_t factor) const {
  Double_t fact[2] = {factor,0} ;
  for (Int_t iStep=0; iStep<fNStep; iStep++) fGrid[iStep]->Scale(fact);
}

inline void AliCFContainer::SetBinContent(Int_t* bin, Int_t step, Double_t value) {
  // sets the content 'value' to the current container, at step 'step'
  // 'bin' is the array of the bin coordinates
  GetGrid(step)->GetGrid()->SetBinContent(bin,value);
}

inline void AliCFContainer::SetBinError(Int_t* bin, Int_t step, Double_t value) {
  // sets the error 'value' to the current container, at step 'step'
  // 'bin' is the array of the bin coordinates
  GetGrid(step)->GetGrid()->SetBinError(bin,value);
}

#endif

