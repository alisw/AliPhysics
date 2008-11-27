// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// multidimensional histogram 
//=============================================================================

#ifndef AliDHN_H
#define AliDHN_H

#include <TH1.h>
class TH2D;
class TAxis;

const Int_t fMaxNdim=10;                     // maximum number of dimensions

//=============================================================================
class AliDHN : public TH1D {

 public:
  AliDHN() : TH1D(), fNdim(0)                  {printf("AliDHN object created\n");}
  AliDHN(Char_t *nam, Int_t ndim, TAxis **ax); // constructor from scratch
  AliDHN(Char_t *filename, Char_t *name);      // constructor from file
  virtual ~AliDHN()                            {printf("AliDHN object %s deleted\n",GetName());}
  Int_t GetNdim() const                     {return fNdim;}
  TAxis *GetAxis(Int_t i) const             {return (TAxis*) &fAxis[i];}

  Int_t Fill(Double_t *xx, Double_t y=1);   // fill histo
  Int_t Fill(Double_t)                      {return -1;} // insufficient number of arguments
  Int_t Fill(Double_t x0, Double_t w)       {return Fill(&x0,w);} // 1-dim histo fill
  Int_t Fill(Double_t x0, Double_t x1, ...);// 2 or more dim histo fill
  Int_t Fill(const char*, Double_t)         {return -1;} // overload TH1

  Int_t Write() const;                      // save histo and axis on file 
  Int_t Write()                             {return ((const AliDHN*)this)->Write();}
  Int_t Write(const char *, Int_t, Int_t)   {return Write();} // overload TObject
  Int_t Write(const char *, Int_t, Int_t) const {return Write();} 

  // project along (integrate over) one axis
  AliDHN  *ProjectAlong(char *nam, Int_t dim, Int_t first=-1, Int_t last=-1);
  // project on 1-dim histogram
  TH1D *ProjectOn(char *nam, Int_t dim, Int_t *first=0, Int_t *last=0) const;
  // project on 1-dim histogram
  TH1D *ProjectOn(char *nam, Int_t dim, Double_t *first, Double_t *last);
  // project on 2-dim histogram
  TH2D *ProjectOn(char *nam, Int_t dim0, Int_t dim1, Int_t *first=0, Int_t *last=0) const;

 protected:

  Int_t       fNdim;                        // number of dimensions
  TAxis       fAxis[fMaxNdim];              // axes
  Int_t       fNbins[fMaxNdim];             // {fAxis[0]->GetNbins(),fAxis[1]->...
  Int_t       fMbins[fMaxNdim];             // {...[fNdim-2]*fNbins[fNdim-1],fNbins[fNdim-1],1}

  static Int_t Albins(Int_t n, TAxis **ax); // product of nbins of ax[0]...ax[n-1]
  Int_t MulToOne(Int_t *k) const;           // calc 1-dim index from n-dim indices
  Int_t MulToOne(Double_t *x);              // calc 1-dim index from n-dim vector
  void  OneToMul(Int_t n, Int_t *k) const;  // calc n-dim indices from 1-dim index

  ClassDef(AliDHN,1)
};
//=============================================================================
#endif
