#ifndef ALITRDDATAARRAYF_H
#define ALITRDDATAARRAYF_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDdataArrayF.h,v */

#include <TMath.h>
 
#include "AliTRDdataArray.h"

/////////////////////////////////////////////////////////////
//                                                         //
//  General container for float data from TRD detector     //
//  segments.                                              //
//  Adapted from AliDigits, origin M.Ivanov                //
//                                                         //
/////////////////////////////////////////////////////////////

class AliTRDarrayF;

class AliTRDdataArrayF : public AliTRDdataArray {

 public:

  AliTRDdataArrayF();
  AliTRDdataArrayF(Int_t nrow, Int_t ncol,Int_t ntime);
  AliTRDdataArrayF(const AliTRDdataArrayF &a);
  virtual ~AliTRDdataArrayF();
  AliTRDdataArrayF &operator=(const AliTRDdataArrayF &a);

  virtual void    Allocate(Int_t nrow, Int_t ncol,Int_t ntime);
  virtual void    Copy(TObject &a);
  virtual void    Compress(Int_t bufferType, Float_t threshold);
  virtual void    Compress(Int_t bufferType); 
  virtual void    Expand();
  virtual Bool_t  First();
  virtual Bool_t  Next(); 
  virtual void    Reset();

          void    SetData(Int_t row, Int_t col, Int_t time, Float_t value);
  virtual void    SetThreshold(Float_t threshold) { fThreshold = threshold; };

  virtual Float_t GetData(Int_t row, Int_t col, Int_t time) const;
  virtual Float_t GetThreshold() const            { return fThreshold;  };

  virtual Int_t   GetSize();
  virtual Int_t   GetDataSize(); 
  virtual Int_t   GetOverThreshold(Float_t threshold);  

 protected:

  inline  void    SetDataFast(Int_t idx1, Int_t idx2, Float_t value); 
  inline  Float_t GetDataFast(Int_t idx1, Int_t idx2) const; 

  Float_t         GetData1(Int_t idx1, Int_t idx2) const; 
  void            Expand1(); 
  void            Compress1(); 
  void            Expand2();
  void            Compress2();
  Bool_t          First0();
  Bool_t          Next0(); 
  Bool_t          First1();
  Bool_t          Next1();

  AliTRDarrayF  *fElements;        // Buffer of 4 bytes floats for the array content
  Float_t        fThreshold;       // Threshold for zero suppression
 
  ClassDef(AliTRDdataArrayF,1)     // Container for float data of one TRD detector segment

};
 
#endif

