#ifndef TRDdataArrayI_H
#define TRDdataArrayI_H

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
  ~AliTRDdataArrayI();

  virtual void   Allocate(Int_t nrow, Int_t ncol, Int_t ntime);
  virtual void   Compress(Int_t bufferType, Int_t threshold);
  virtual void   Compress(Int_t bufferType); 
  virtual void   Expand();
  virtual Bool_t First();
  virtual Bool_t Next(); 
  virtual void   Reset();

  inline  void   SetData(Int_t row, Int_t col, Int_t time, Int_t value);
  virtual void   SetThreshold(Int_t threshold) { fThreshold = threshold; };

  virtual Int_t  GetData(Int_t row, Int_t col, Int_t time);
  virtual Int_t  GetThreshold()                { return fThreshold;  };

  virtual Int_t  GetSize();
  virtual Int_t  GetDataSize(); 
  virtual Int_t  GetOverThreshold(Int_t threshold);  

 protected:

  inline  void   SetDataFast(Int_t idx1, Int_t idx2, Int_t value); 
  inline  Int_t  GetDataFast(Int_t idx1, Int_t idx2); 

  Int_t          GetData1(Int_t idx1, Int_t idx2); 
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
 

//_____________________________________________________________________________
inline Int_t AliTRDdataArrayI::GetDataFast(Int_t idx1, Int_t idx2)
{
  //
  // Returns the value at a given position in the array
  //
  
  return fElements->At(fIndex->At(idx2) + idx1); 

}

//_____________________________________________________________________________
inline void AliTRDdataArrayI::SetData(Int_t row, Int_t col, Int_t time
                                                          , Int_t value)
{
  //
  // Sets the data value at a given position of the array
  // Includes boundary checking
  //

  if ((row >= 0) && (col >= 0) && (time >= 0)) {
    Int_t idx1 = GetIdx1(row,col);
    if ((idx1 >= 0) && (time < fNdim2)) {
       SetDataFast(idx1,time,value);
    }
    else {
      if (idx1 >= 0) {
        TObject::Error("SetData"
                      ,"time %d out of bounds (size: %d, this: 0x%08x)"
                      ,time,fNdim2,this);
      }
    }
  }

}

//_____________________________________________________________________________
inline void  AliTRDdataArrayI::SetDataFast(Int_t idx1, Int_t idx2, Int_t value)
{
  //
  // Set the value at a given position in the array
  //

  if ((idx1 < 0) || (idx1 >= fNdim1) || 
      (idx2 < 0) || (idx2 >= fNdim2)) { 
    TObject::Error("SetDataFast"
                  ,"idx1 %d  idx2 %d out of bounds (size: %d x %d, this: 0x%08x)"
                  ,idx1,idx2,fNdim1,fNdim2,this);
  }

  (*fElements)[fIndex->fArray[idx2] + idx1] = value; 

}

#endif

