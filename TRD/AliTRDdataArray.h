#ifndef TRDdataArray_H
#define TRDdataArray_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
 
#include   "AliTRDarrayI.h"
#include   "AliTRDarrayF.h"
#include   "AliTRDsegmentID.h"

/////////////////////////////////////////////////////////////
//  General container for data from TRD detector segments  //
//  Adapted from AliDigits, origin M.Ivanov                //
/////////////////////////////////////////////////////////////

class AliTRDdataArray : public AliTRDsegmentID {

 public:

  AliTRDdataArray();
  AliTRDdataArray(Int_t nrow, Int_t ncol,Int_t ntime);
  ~AliTRDdataArray();

  virtual void   Allocate(Int_t nrow, Int_t ncol,Int_t ntime);
  virtual void   Reset();

  virtual Int_t  GetNRow()                     { return fNrow;       };
  virtual Int_t  GetNCol()                     { return fNcol;       };
  virtual Int_t  GetNtime()                    { return fNtime;      };

  inline  Int_t  GetIndex(Int_t row, Int_t col, Int_t time);

 protected:

  inline  Int_t  GetIdx1(Int_t row, Int_t col);
  inline  Bool_t CheckBounds(const char *where, Int_t idx1, Int_t idx2);
  inline  Bool_t OutOfBoundsError(const char *where, Int_t idx1, Int_t idx2);
 
  Int_t          fNrow;            // Number of rows of the detector segement
  Int_t          fNcol;            // Number of columns of the detector segment
  Int_t          fNtime;           // Number of timebins of the detector segment

  Int_t          fNdim1;           // First dimension of the array (row * column)
  Int_t          fNdim2;           // Second dimension of the array (time, not compressed) 

  AliTRDarrayI  *fIndex;           // Index position of column
  Int_t          fBufType;         // Type of the buffer - defines the compression algorithm  
  Int_t          fNelems;          // Total number of elements 
  Int_t          fCurrentIdx1;     // !Current index 1
  Int_t          fCurrentIdx2;     // !Current index 2
  Int_t          fCurrentIndex;    // !Current index in field
 
  ClassDef(AliTRDdataArray,1)      // Data container for one TRD detector segment

};
 
//_____________________________________________________________________________
inline Int_t AliTRDdataArray::GetIndex(Int_t row, Int_t col, Int_t time)
{
  //
  // Maps the row/column/time into one number
  // 

  if (time > fNtime) {
    TObject::Error("GetIdx1"
                  ,"time %d out of bounds (size: %d, this: 0x%08x)"
                  ,time,fNtime,this);
    return -1;
  }  
  
  return time * fNrow*fNcol + GetIdx1(row,col);

}

//_____________________________________________________________________________
inline Int_t AliTRDdataArray::GetIdx1(Int_t row, Int_t col)
{
  //
  // Maps the two-dimensional row/column plane into an one-dimensional array.
  //

  if (row >= fNrow) {
    TObject::Error("GetIdx1"
                  ,"row %d out of bounds (size: %d, this: 0x%08x)"
                  ,row,fNrow,this);
    return -1;
  }  

  if (col >= fNcol) {
    TObject::Error("GetIdx1"
                  ,"col %d out of bounds (size: %d, this: 0x%08x)"
                  ,col,fNcol,this);
    return -1;
  }  

  return row + col * fNrow;

}

//_____________________________________________________________________________
inline Bool_t AliTRDdataArray::CheckBounds(const char *where
                                          , Int_t idx1, Int_t idx2) 
{
  //
  // Does the boundary checking
  //

  if ((idx2 >= fNdim2) || (idx2 < 0)) 
    return OutOfBoundsError(where,idx1,idx2);

  Int_t index = (*fIndex).At(idx2) + idx1;
  if ((index < 0) || (index > fNelems)) 
    return OutOfBoundsError(where,idx1,idx2);

  return kTRUE;  

}

//_____________________________________________________________________________
inline Bool_t AliTRDdataArray::OutOfBoundsError(const char *where
                                               , Int_t idx1, Int_t idx2) 
{
  //
  // Generate an out-of-bounds error. Always returns false.
  //

  TObject::Error(where, "idx1 %d  idx2 %d out of bounds (size: %d x %d, this: 0x%08x)"
	   ,idx1,idx2,fNdim1,fNdim2,this);

  return kFALSE;

}

#endif

