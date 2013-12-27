#ifndef ALIPIDPERFCONTAINER_H
#define ALIPIDPERFCONTAINER_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliPIDperfContainer.h 49869 2013-07-10 04:49:51Z fnoferin $ */

/////////////////////////////////////////////////
//                                             //
//         output v2-VZERO Class               //
//           noferini@bo.infn.it               //
/////////////////////////////////////////////////

#include "TH2F.h"
#include "TClonesArray.h"
#include "TArrayI.h"
#include "TArrayF.h"

class AliPIDperfContainer : public TNamed
{
 public:
  AliPIDperfContainer(const char *name, Int_t nvar,const Int_t* binVar);
  AliPIDperfContainer();
  ~AliPIDperfContainer();
  AliPIDperfContainer(const AliPIDperfContainer &old);
  AliPIDperfContainer& operator=(const AliPIDperfContainer &source);

  Int_t GetNhistos() const {return fQA->GetEntries();};
  Int_t GetNspecies() const;
  TH2F *GetQA(Int_t histo) const {return ((TH2F *) fQA->At(histo));};
  TH2F *GetQA(Int_t species,Float_t x[]) const;
  TH2F *GetQA(Int_t species,Float_t xMin[],Float_t xMax[]) const;
  void DirectFill(Int_t histo,Float_t var1,Float_t var2){GetQA(histo)->Fill(var1,var2);};
  void Fill(Int_t species,Float_t var1,Float_t var2,Float_t x[]);

  void AddSpecies(const char *name,Int_t nXbin,const Double_t *xbin,Int_t nYbin,const Double_t *ybin);

  const char *GetSpeciesName(Int_t species){if(species<GetNspecies()) return GetQA(species*GetNhistos()/GetNspecies())->GetName();else return "";};

  Int_t Add(const AliPIDperfContainer *oth);

  Int_t GetNvar() const {return fNbinVar->GetSize();};
  Int_t GetNbinVar(Int_t ivar) const {return (*fNbinVar)[ivar];};

  void SetVarRange(Int_t ivar,Float_t xMin,Float_t xMax);
  void SetVarName(Int_t ivar,const char *name){TNamed *atemp = (TNamed *) fNameVar->At(ivar); atemp->SetName(name);};

  Float_t GetXmin(Int_t ivar) const {return (*fXmin)[ivar];};
  Float_t GetXmax(Int_t ivar) const {return (*fXmax)[ivar];};
  const char *GetVarName(Int_t ivar) const {TNamed *atemp = (TNamed *) fNameVar->At(ivar); return atemp->GetName();};

  Int_t GetBin(Int_t ivar,Float_t x) const {return Int_t((x-(*fXmin)[ivar])/((*fXmax)[ivar]-(*fXmin)[ivar])*(*fNbinVar)[ivar]);};

  Long64_t Merge(TCollection* list);

  void Reset();

  void SetTitleX(const char *title){snprintf(fTitleX,100,"%s",title);};
  void SetTitleY(const char *title){snprintf(fTitleY,100,"%s",title);};
  const char *GetTitleX() const {return fTitleX;};
  const char *GetTitleY() const {return fTitleY;};

 private:
  TArrayI *fNbinVar;
  TArrayF *fXmin,*fXmax;
  TClonesArray *fNameVar;

  TClonesArray *fQA;

  char fTitleX[100]; // title for X-axis
  char fTitleY[100]; // title for Y-axis

  ClassDef(AliPIDperfContainer,1)  // qa vzero outuput object
};
#endif


