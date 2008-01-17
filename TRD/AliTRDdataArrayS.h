#ifndef ALITRDDATAARRAYS_H
#define ALITRDDATAARRAYS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDdataArrayS.h,v Exp $ */
 
#include   "AliTRDdataArray.h"

/////////////////////////////////////////////////////////////
//                                                         //
//  General container for integer data from TRD detector   //
//  segments.                                              //
//  Adapted from AliDigits, origin M.Ivanov                //
//                                                         //
/////////////////////////////////////////////////////////////

class AliTRDarrayS;

class AliTRDdataArrayS : public AliTRDdataArray {

 public:

  AliTRDdataArrayS();
  AliTRDdataArrayS(Int_t nrow, Int_t ncol, Int_t ntime);
  AliTRDdataArrayS(const AliTRDdataArrayS &a);
  virtual ~AliTRDdataArrayS();
  AliTRDdataArrayS &operator=(const AliTRDdataArrayS &a);

  virtual void    Allocate(Int_t nrow, Int_t ncol, Int_t ntime);
  virtual void    Copy(TObject &a) const;
  virtual void    Compress(Int_t bufferType, Short_t threshold);
  virtual void    Compress(Int_t bufferType); 
  virtual void    Expand();
  virtual Bool_t  First();
  virtual Bool_t  Next(); 
  virtual void    Reset();

          void    SetData(Int_t row, Int_t col, Int_t time, Short_t value);
          void    SetDataUnchecked(Int_t row, Int_t col, Int_t time, Short_t value)
                                  { SetDataFast(GetIdx1Unchecked(row,col),time,value); };

  virtual void    SetThreshold(Short_t threshold) { fThreshold = threshold; };

  virtual Short_t GetData(Int_t row, Int_t col, Int_t time) const;
          Short_t GetDataUnchecked(Int_t row, Int_t col, Int_t time) const
                                 { return GetDataFast(GetIdx1Unchecked(row,col),time); };

  virtual Short_t GetThreshold() const            { return fThreshold;      };

  virtual Int_t   GetSize() const;
  virtual Int_t   GetDataSize() const; 
  virtual Int_t   GetOverThreshold(Short_t threshold);  

 protected:

          void    SetDataFast(Int_t idx1, Int_t idx2, Short_t value);
          Short_t GetDataFast(Int_t idx1, Int_t idx2) const;

  Short_t         GetData1(Int_t idx1, Int_t idx2) const; 
  void            Expand1(); 
  void            Compress1(); 
  void            Expand2();
  void            Compress2();
  Bool_t          First0();
  Bool_t          Next0(); 
  Bool_t          First1();
  Bool_t          Next1();
 
  AliTRDarrayS   *fElements;        // Buffer of 2 bytes integers for the array content
  Short_t         fThreshold;       // Threshold for zero suppression
 
  ClassDef(AliTRDdataArrayS,1)      // Container for short data of one TRD detector segment

};
 
#endif


