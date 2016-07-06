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
  Float_t GetDataByAdcCol(Int_t row, Int_t col, Int_t time) const
               {return fSignal[(row*fNumberOfChannels+col)*fNtime+time];};
  void    SetDataByAdcCol(Int_t row, Int_t col, Int_t time, Float_t value)
              {fSignal[(row*fNumberOfChannels+col)*fNtime+time]=value;};
  Bool_t  HasData() const {return fNtime ? 1 : 0;};
  Int_t   GetDim() const {return fNdim;};
  Int_t   GetOverThreshold(Float_t threshold) const;
  void    Compress(Float_t minval);
  void    Expand();
  void    Reset();
  Float_t GetData(Int_t nrow, Int_t ncol, Int_t ntime) const;
  void    SetData(Int_t nrow, Int_t ncol, Int_t ntime, Float_t value);
  static  void    CreateLut(); 

 protected:

  Int_t    fNdet;      //ID number of the chamber
  Int_t    fNrow;      //Number of rows of the chamber
  Int_t    fNcol;      //Number of columns of the chamber
  Int_t    fNumberOfChannels;  //  Number of MCM channels per row
  Int_t    fNtime;     //Number of time bins
  Int_t    fNdim;      //Dimension of the array
  Float_t *fSignal;    //[fNdim]  //Pointer to signals
  static Short_t *fgLutPadNumbering;   //  [fNcol] Look Up Table        

  ClassDef(AliTRDarraySignal,3)  //Signal container class
    
};

// inline definitions

//________________________________________________________________________________
inline Float_t AliTRDarraySignal::GetData(Int_t nrow, Int_t ncol, Int_t ntime) const
{
  //
  // Get the data using the pad numbering.
  // To access data using the mcm scheme use instead
  // the method GetDataByAdcCol
  //

  Int_t corrcolumn = fgLutPadNumbering[ncol];

  return fSignal[(nrow*fNumberOfChannels+corrcolumn)*fNtime+ntime];

}
//________________________________________________________________________________
inline void AliTRDarraySignal::SetData(Int_t nrow, Int_t ncol, Int_t ntime, Float_t value)
{
  //
  // Set the data using the pad numbering.
  // To write data using the mcm scheme use instead
  // the method SetDataByAdcCol
  //

  Int_t colnumb = fgLutPadNumbering[ncol];

  fSignal[(nrow*fNumberOfChannels+colnumb)*fNtime+ntime]=value;

}

#endif
