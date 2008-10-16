// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

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
  TAxis *GetAxis(Int_t i)                   {return &fAxis[i];}
  void Fill(Double_t *xx, Double_t y=1);    // fill histo
  void Fill(Double_t x0=0, Double_t x1=0, 
	    Double_t x2=0, Double_t x3=0, 
	    Double_t x4=0, Double_t x5=0, 
	    Double_t x6=0, Double_t x7=0, 
	    Double_t x8=0, Double_t x9=0,
	    Double_t x10=0) {               // fill histo; fNdim-th arg is weight
    Double_t xx[fMaxNdim+1] = {x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10};
    Fill(xx,xx[fNdim]);}
  Int_t Write();                            // save histo and axis on file 
  Int_t Write(const char *, Int_t, Int_t)   {return Write();} 
  // project along (integrate over) one axis
  AliDHN  *ProjectAlong(char *nam, Int_t dim, Int_t first=-1, Int_t last=-1);  
  // project on 1-dim histogram
  TH1D *ProjectOn(char *nam, Int_t dim, Int_t *first=0, Int_t *last=0);
  // project on 1-dim histogram
  TH1D *ProjectOn(char *nam, Int_t dim, Double_t *first, Double_t *last);
  // project on 2-dim histogram
  TH2D *ProjectOn(char *nam, Int_t dim0, Int_t dim1, Int_t *first=0, Int_t *last=0);

 protected:

  Int_t       fNdim;                        // number of dimensions
  TAxis       fAxis[fMaxNdim];              // axes
  Int_t       fNbins[fMaxNdim];             // {fAxis[0]->GetNbins(),fAxis[1]->...
  Int_t       fMbins[fMaxNdim];             // {...[fNdim-2]*fNbins[fNdim-1],fNbins[fNdim-1],1}

  static Int_t Albins(Int_t n, TAxis **ax); // product of nbins of ax[0]...ax[n-1]
  Int_t MulToOne(Int_t *k) const;           // calc 1-dim index from n-dim indices
  Int_t MulToOne(Double_t *x);              // calc 1-dim index from n-dim vector
  void  OneToMul(Int_t n, Int_t *k);        // calc n-dim indices from 1-dim index

  ClassDef(AliDHN,1)
};
//=============================================================================
#endif
