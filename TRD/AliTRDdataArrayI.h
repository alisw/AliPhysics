#ifndef ALITRDDATAARRAYI_H
#define ALITRDDATAARRAYI_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDdataArrayI.h,v */
 
#include   "AliTRDdataArray.h"

/////////////////////////////////////////////////////////////
//                                                         //
//  General container for integer data from TRD detector   //
//  segments.                                              //
//  Adapted from AliDigits, origin M.Ivanov                //
//                                                         //
/////////////////////////////////////////////////////////////

class AliTRDdataArrayI : public AliTRDdataArray {

 public:

  AliTRDdataArrayI();
  AliTRDdataArrayI(Int_t nrow, Int_t ncol, Int_t ntime);
  AliTRDdataArrayI(const AliTRDdataArrayI &a);
  virtual ~AliTRDdataArrayI();
  AliTRDdataArrayI &operator=(const AliTRDdataArrayI &a);

  virtual void   Allocate(Int_t nrow, Int_t ncol, Int_t ntime);
  virtual void   Copy(TObject &a);
  virtual void   Compress(Int_t bufferType, Int_t threshold);
  virtual void   Compress(Int_t bufferType); 
  virtual void   Expand();
  virtual Bool_t First();
  virtual Bool_t Next(); 
  virtual void   Reset();

          void   SetData(Int_t row, Int_t col, Int_t time, Int_t value);
  virtual void   SetThreshold(Int_t threshold) { fThreshold = threshold; };

  virtual Int_t  GetData(Int_t row, Int_t col, Int_t time) const;
  virtual Int_t  GetThreshold() const          { return fThreshold;  };

  virtual Int_t  GetSize();
  virtual Int_t  GetDataSize(); 
  virtual Int_t  GetOverThreshold(Int_t threshold);  

 protected:

  inline  void   SetDataFast(Int_t idx1, Int_t idx2, Int_t value); 
  inline  Int_t  GetDataFast(Int_t idx1, Int_t idx2) const; 

  Int_t          GetData1(Int_t idx1, Int_t idx2) const; 
  void           Expand1(); 
  void           Compress1(); 
  void           Expand2();
  void           Compress2();
  Bool_t         First0();
  Bool_t         Next0(); 
  Bool_t         First1();
  Bool_t         Next1();
 
  AliTRDarrayI  *fElements;        // Buffer of 4 bytes integers for the array content
  Int_t          fThreshold;       // Threshold for zero suppression
 
  ClassDef(AliTRDdataArrayI,1)     // Container for integer data of one TRD detector segment

};
 
#endif


