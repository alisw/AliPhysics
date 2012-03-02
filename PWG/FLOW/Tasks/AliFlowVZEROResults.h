#ifndef ALIFLOWVZERORESULTS_H
#define ALIFLOWVZERORESULTS_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliFlowVZEROResults.h 49869 2011-05-18 04:49:51Z fnoferin $ */

/////////////////////////////////////////////////
//                                             //
//         output v2-VZERO Class               //
//           noferini@bo.infn.it               //
/////////////////////////////////////////////////

#include "TProfile.h"
#include "TClonesArray.h"
#include "TArrayI.h"
#include "TArrayF.h"

class AliFlowVZEROResults : public TNamed
{
 public:
  AliFlowVZEROResults(const char *name,const Int_t nvar,const Int_t* binVar);
  AliFlowVZEROResults();
  ~AliFlowVZEROResults();
  AliFlowVZEROResults(const AliFlowVZEROResults &old);
  AliFlowVZEROResults& operator=(const AliFlowVZEROResults &source);

  Int_t GetNhistos() const {return fV2->GetEntries();};
  Int_t GetNspecies() const;
  TProfile *GetV2(Int_t histo) const {return ((TProfile *) fV2->At(histo));};
  TProfile *GetV2(Int_t species,Float_t x[]) const;
  TProfile *GetV2(Int_t species,Float_t xMin[],Float_t xMax[]) const;
  void DirectFill(Int_t histo,Float_t pt,Float_t v2){GetV2(histo)->Fill(pt,v2);};
  void Fill(Int_t species,Float_t pt,Float_t v2,Float_t x[]);

  void AddSpecies(const char *name,Int_t nXbin,const Double_t *bin);

  const char *GetSpeciesName(Int_t species){return GetV2(species)->GetName();};

  Int_t Add(const AliFlowVZEROResults *oth);

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

  TClonesArray *fV2;


  ClassDef(AliFlowVZEROResults,1)  // v2 vzero outuput object
};
#endif


