#ifndef ALITRDARRAYADC_H
#define ALITRDARRAYADC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice */ 

/* $Id: AliTRDarrayADC.h 23387 2008-01-17 17:25:16Z cblume $ */

///////////////////////////////////////////////
//                                           //
// Container class for ADC values            //
//                                           // 
///////////////////////////////////////////////

#include <TObject.h>

class AliTRDarrayADC: public TObject
{
 public:

  AliTRDarrayADC();
  AliTRDarrayADC(Int_t nrow, Int_t ncol, Int_t ntime);
  AliTRDarrayADC(const AliTRDarrayADC &b);
  ~AliTRDarrayADC();
  AliTRDarrayADC &operator=(const AliTRDarrayADC &b);

  void    Allocate(Int_t nrow, Int_t ncol, Int_t ntime);
  void    SetNdet(Int_t ndet) {fNdet=ndet;};  
  Int_t   GetNdet()  const {return fNdet;};
  void    SetData(Int_t nrow, Int_t ncol, Int_t ntime, Short_t value)
                        {fADC[(nrow*fNcol+ncol)*fNtime+ntime]=value;};
  Bool_t  HasData() const {return fNtime ? 1 : 0;};
  Short_t GetData(Int_t nrow, Int_t ncol, Int_t ntime) const
                       {return fADC[(nrow*fNcol+ncol)*fNtime+ntime];};
  Short_t GetDataB(Int_t nrow, Int_t ncol, Int_t ntime) const;
  UChar_t GetPadStatus(Int_t nrow, Int_t ncol, Int_t ntime) const;
  void    SetPadStatus(Int_t nrow, Int_t ncol, Int_t ntime, UChar_t status);
  Bool_t  IsPadCorrupted(Int_t nrow, Int_t ncol, Int_t ntime);
  void    Compress();
  void    Expand();
  Int_t   GetNtime() const {return fNtime;};
  Int_t   GetNrow() const {return fNrow;};
  Int_t   GetNcol() const {return fNcol;};
  Int_t   GetDim() const {return fNAdim;};
  void    DeleteNegatives();

 protected:

  Int_t fNdet;    //ID number of the chamber
  Int_t fNrow;    //Number of rows
  Int_t fNcol;    //Number of columns
  Int_t fNtime;   //Number of time bins
  Int_t fNAdim;   //Dimension of the ADC array
  Short_t* fADC;  //[fNAdim]   //Pointer to adc values

  ClassDef(AliTRDarrayADC,1) //ADC container class
    
};
#endif 
