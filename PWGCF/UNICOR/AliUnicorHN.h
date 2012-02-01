#ifndef ALIUNICORHN_H
#define ALIUNICORHN_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// multidimensional histogram 
//=============================================================================

#include <TH1.h>
class TH2D;
class TAxis;

//=============================================================================
class AliUnicorHN : public TH1D {

 public:
  AliUnicorHN(const char *nam="muhi", Int_t ndim=0, TAxis **ax=0);     // constructor
  AliUnicorHN(TRootIOCtor *) : TH1D(), fNdim(0) {for (int i=0; i<fgkMaxNdim; i++) fNbins[i]=fMbins[i]=0;}  // default constructor
  virtual ~AliUnicorHN() {}                                            // destructor
  static AliUnicorHN* Retrieve(const char *filnam, const char *nam);   // read from file

  Int_t GetNdim() const                     {return fNdim;}
  TAxis *GetAxis(Int_t i) const             {return (TAxis*) &fAxis[i];}

  Int_t Fill(Double_t *xx, Double_t y=1);   // fill histo
  Int_t Fill(Double_t)                      {return -1;} // insufficient number of arguments
  Int_t Fill(Double_t x0, Double_t w)       {Double_t x[1]={x0}; return Fill(x,w);} // 1-dim histo fill
  Int_t Fill(Double_t x0, Double_t x1, ...);// 2 or more dim histo fill
  Int_t Fill(const char*, Double_t)         {return -1;} // overload TH1

  Int_t Save() const;                      // save histo and axis on file 

  // project along (integrate over) one axis
  AliUnicorHN  *ProjectAlong(const char *nam, Int_t dim, Int_t first=0, Int_t last=0);
  // project on 1-dim histogram
  TH1D *ProjectOn(const char *nam, Int_t dim, const Int_t * const first=0, const Int_t * const last=0) const;
  TH1D *ProjectOn(const char *nam, Int_t dim, const Double_t * const first, const Double_t * const last) const;
  // project on 2-dim histogram
  TH2D *ProjectOn(const char *nam, Int_t dim0, Int_t dim1, const Int_t * const first=0, const Int_t * const last=0) const;
  TH2D *ProjectOn(const char *nam, Int_t dim0, Int_t dim1, const Double_t * const first, const Double_t * const last) const;
  void OneToMul(Int_t n, Int_t *k) const;      // calc n-dim indices from 1-dim index

 protected:

  static const Int_t fgkMaxNdim=10;             // maximum number of dimensions
  Int_t              fNdim;                     // number of dimensions
  TAxis              fAxis[fgkMaxNdim];         // axes
  Int_t              fNbins[fgkMaxNdim];        // {fAxis[0]->GetNbins(),fAxis[1]->...
  Int_t              fMbins[fgkMaxNdim];        // {...[fNdim-2]*fNbins[fNdim-1],fNbins[fNdim-1],1}
  static Int_t Albins(Int_t n, TAxis **ax);     // product of nbins of ax[0]...ax[n-1]
  Int_t MulToOne(const Int_t * const k) const;  // calc 1-dim index from n-dim indices
  Int_t MulToOne(Double_t *x);                  // calc 1-dim index from n-dim vector

  ClassDef(AliUnicorHN,1)
};
//=============================================================================
#endif
