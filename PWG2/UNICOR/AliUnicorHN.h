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
  AliUnicorHN() : TH1D(), fNdim(0) {}               // default contructor
  AliUnicorHN(Char_t *nam, Int_t ndim, TAxis **ax); // constructor from scratch
  AliUnicorHN(Char_t *filename, Char_t *name);      // constructor from file
  virtual ~AliUnicorHN() {}                         // destructor
  Int_t GetNdim() const                     {return fNdim;}
  TAxis *GetAxis(Int_t i) const             {return (TAxis*) &fAxis[i];}

  Int_t Fill(Double_t *xx, Double_t y=1);   // fill histo
  Int_t Fill(Double_t)                      {return -1;} // insufficient number of arguments
  Int_t Fill(Double_t x0, Double_t w)       {return Fill(&x0,w);} // 1-dim histo fill
  Int_t Fill(Double_t x0, Double_t x1, ...);// 2 or more dim histo fill
  Int_t Fill(const char*, Double_t)         {return -1;} // overload TH1

  Int_t Write() const;                      // save histo and axis on file 
  Int_t Write()                             {return ((const AliUnicorHN*)this)->Write();}
  Int_t Write(const char *, Int_t, Int_t)   {return Write();} // overload TObject
  Int_t Write(const char *, Int_t, Int_t) const {return Write();} 

  // project along (integrate over) one axis
  AliUnicorHN  *ProjectAlong(char *nam, Int_t dim, Int_t first=-1, Int_t last=-1);
  // project on 1-dim histogram
  TH1D *ProjectOn(char *nam, Int_t dim, const Int_t * const first=0, const Int_t * const last=0) const;
  // project on 1-dim histogram
  TH1D *ProjectOn(char *nam, Int_t dim, const Double_t * const first, const Double_t * const last);
  // project on 2-dim histogram
  TH2D *ProjectOn(char *nam, Int_t dim0, Int_t dim1, const Int_t * const first=0, const Int_t * const last=0) const;

 protected:

  static const Int_t fgkMaxNdim=10;             // maximum number of dimensions
  Int_t              fNdim;                     // number of dimensions
  TAxis              fAxis[fgkMaxNdim];         // axes
  Int_t              fNbins[fgkMaxNdim];        // {fAxis[0]->GetNbins(),fAxis[1]->...
  Int_t              fMbins[fgkMaxNdim];        // {...[fNdim-2]*fNbins[fNdim-1],fNbins[fNdim-1],1}
  static Int_t Albins(Int_t n, TAxis **ax);     // product of nbins of ax[0]...ax[n-1]
  Int_t MulToOne(const Int_t * const k) const;  // calc 1-dim index from n-dim indices
  Int_t MulToOne(Double_t *x);                  // calc 1-dim index from n-dim vector
  void  OneToMul(Int_t n, Int_t *k) const;      // calc n-dim indices from 1-dim index

  ClassDef(AliUnicorHN,1)
};
//=============================================================================
#endif
