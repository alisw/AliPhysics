#ifndef ALITRDDATAARRAY_H
#define ALITRDDATAARRAY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
 
#include   "AliTRDsegmentID.h"

/////////////////////////////////////////////////////////////
//  General container for data from TRD detector segments  //
//  Adapted from AliDigits, origin M.Ivanov                //
/////////////////////////////////////////////////////////////

class AliTRDarrayI;

class AliTRDdataArray : public AliTRDsegmentID {

 public:

  AliTRDdataArray();
  AliTRDdataArray(Int_t nrow, Int_t ncol,Int_t ntime);
  AliTRDdataArray(const AliTRDdataArray &d);
  virtual ~AliTRDdataArray();
  AliTRDdataArray &operator=(const AliTRDdataArray &d);

  virtual void   Copy(TObject &d);
  virtual void   Allocate(Int_t nrow, Int_t ncol,Int_t ntime);
  virtual void   Reset();

  virtual Int_t  GetNRow() const               { return fNrow;       };
  virtual Int_t  GetNCol() const               { return fNcol;       };
  virtual Int_t  GetNtime() const              { return fNtime;      };
          Int_t  GetIndex(Int_t row, Int_t col, Int_t time) const;
          Int_t  GetBufType() const            { return fBufType;    };
  virtual Int_t  GetNelems() const             { return fNelems;     };

 protected:

          Int_t  GetIdx1(Int_t row, Int_t col) const;
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
 
#endif

