#ifndef ALITRDARRAYSIGNAL_H
#define ALITRDARRAYSIGNAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice */ 

/* $Id: AliTRDarraySignal.h 23387 2008-01-17 17:25:16Z cblume $ */

/////////////////////////////////////////////
//                                         //
// Container Class for Signals             //
//                                         //
/////////////////////////////////////////////

#include <TObject.h>

class AliTRDarraySignal: public TObject
{

 public:

  AliTRDarraySignal();
  AliTRDarraySignal(Int_t nrow, Int_t ncol,Int_t ntime);
  AliTRDarraySignal(const AliTRDarraySignal &d); 
  ~AliTRDarraySignal();
  AliTRDarraySignal &operator=(const AliTRDarraySignal &d); 

  void    Allocate(Int_t nrow, Int_t ncol, Int_t ntime);
  void    SetNdet(Int_t ndet) {fNdet=ndet;};  
  Int_t   GetNdet()  const {return fNdet;};
  Int_t   GetNrow()  const {return fNrow;};
  Int_t   GetNcol()  const {return fNcol;};
  Int_t   GetNtime() const {return fNtime;};
  Float_t GetData(Int_t row, Int_t col, Int_t time) const
               {return fSignal[(row*fNcol+col)*fNtime+time];};
  void    SetData(Int_t row, Int_t col, Int_t time, Float_t value)
              {fSignal[(row*fNcol+col)*fNtime+time]=value;};
  Bool_t  HasData() const {return fNtime ? 1 : 0;};
  Int_t   GetDim() const {return fNdim;};
  Int_t   GetOverThreshold(Float_t threshold);
  void    Compress(Float_t minval);
  void    Expand();
  void    Reset();

 protected:

  Int_t    fNdet;      //ID number of the chamber
  Int_t    fNrow;      //Number of rows of the chamber
  Int_t    fNcol;      //Number of columns of the chamber
  Int_t    fNtime;     //Number of time bins
  Int_t    fNdim;      //Dimension of the array
  Float_t *fSignal;    //[fNdim]  //Pointer to signals 

  ClassDef(AliTRDarraySignal,1)  //Signal container class
    
};
#endif
