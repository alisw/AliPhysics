#ifndef ALIFLOWVZEROQA_H
#define ALIFLOWVZEROQA_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliFlowVZEROQA.h 49869 2011-05-18 04:49:51Z fnoferin $ */

/////////////////////////////////////////////////
//                                             //
//         output v2-VZERO Class               //
//           noferini@bo.infn.it               //
/////////////////////////////////////////////////

#include "TH2F.h"
#include "TClonesArray.h"
#include "TArrayI.h"
#include "TArrayF.h"

class AliFlowVZEROQA : public TNamed
{
 public:
  AliFlowVZEROQA(const char *name,const Int_t nvar,const Int_t* binVar);
  AliFlowVZEROQA();
  ~AliFlowVZEROQA();
  AliFlowVZEROQA(const AliFlowVZEROQA &old);
  AliFlowVZEROQA& operator=(const AliFlowVZEROQA &source);

  Int_t GetNhistos() const {return fQA->GetEntries();};
  Int_t GetNspecies() const;
  TH2F *GetQA(Int_t histo) const {return ((TH2F *) fQA->At(histo));};
  TH2F *GetQA(Int_t species,Float_t x[]) const;
  TH2F *GetQA(Int_t species,Float_t xMin[],Float_t xMax[]) const;
  void DirectFill(Int_t histo,Float_t var1,Float_t var2){GetQA(histo)->Fill(var1,var2);};
  void Fill(Int_t species,Float_t var1,Float_t var2,Float_t x[]);

  void AddSpecies(const char *name,Int_t nXbin,const Double_t *xbin,Int_t nYbin,const Double_t *ybin);

  const char *GetSpeciesName(Int_t species){return GetQA(species)->GetName();};

  Int_t Add(const AliFlowVZEROQA *oth);

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

 private:
  TArrayI *fNbinVar;
  TArrayF *fXmin,*fXmax;
  TClonesArray *fNameVar;

  TClonesArray *fQA;


  ClassDef(AliFlowVZEROQA,1)  // qa vzero outuput object
};
#endif


