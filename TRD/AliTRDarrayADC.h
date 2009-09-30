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

class AliTRDSignalIndex;
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
  void    SetDataByAdcCol(Int_t nrow, Int_t ncol, Int_t ntime, Short_t value)
                         {fADC[(nrow*fNumberOfChannels+ncol)*fNtime+ntime]=value;}
  Bool_t  HasData() const {return fNtime ? 1 : 0;};
  Short_t GetDataByAdcCol(Int_t nrow, Int_t ncol, Int_t ntime) const
                         {return fADC[(nrow*fNumberOfChannels+ncol)*fNtime+ntime];};
  inline  void GetData(Int_t r, Int_t c, Int_t t, Int_t n, Short_t *vals) const;
  Short_t GetDataBits(Int_t nrow, Int_t ncol, Int_t ntime) const;
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
  void    Reset();
  void    ConditionalReset(AliTRDSignalIndex* idx);
  inline  Short_t GetData(Int_t nrow, Int_t ncol, Int_t ntime) const;
  inline  void    SetData(Int_t nrow, Int_t ncol, Int_t ntime, Short_t value);
  static  void    CreateLut(); 

 protected:

  Int_t fNdet;    //ID number of the chamber
  Int_t fNrow;    //Number of rows
  Int_t fNcol;    //Number of columns(pads)
  Int_t fNumberOfChannels;  //  Number of MCM channels per row
  Int_t fNtime;   //Number of time bins
  Int_t fNAdim;   //Dimension of the ADC array
  Short_t* fADC;  //[fNAdim]   //Pointer to adc values
  static Short_t *fgLutPadNumbering;   //  [fNcol] Look Up Table

  ClassDef(AliTRDarrayADC,4) //ADC container class
    
};

//________________________________________________________________________________
Short_t AliTRDarrayADC::GetData(Int_t nrow, Int_t ncol, Int_t ntime) const
{
  //
  // Get the data using the pad numbering.
  // To access data using the mcm scheme use instead
  // the method GetDataByAdcCol
  //

  Int_t corrcolumn = fgLutPadNumbering[ncol];

  return fADC[(nrow*fNumberOfChannels+corrcolumn)*fNtime+ntime];

}
//________________________________________________________________________________
void AliTRDarrayADC::SetData(Int_t nrow, Int_t ncol, Int_t ntime, Short_t value)
{
  //
  // Set the data using the pad numbering.
  // To write data using the mcm scheme use instead
  // the method SetDataByAdcCol
  //

  Int_t colnumb = fgLutPadNumbering[ncol];

  fADC[(nrow*fNumberOfChannels+colnumb)*fNtime+ntime] = value;

}

void AliTRDarrayADC::GetData(Int_t r, Int_t c, Int_t t, Int_t n, Short_t *vals) const
{
  Int_t colNum = fgLutPadNumbering[c];
  for(Int_t ic=n, idx = (r*fNumberOfChannels+colNum)*fNtime+t; ic--; idx+=fNtime) vals[ic] = fADC[idx];
 }

#endif 

