/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Base class of a general container for data of a TRD detector segment.    //
//  Adapted from AliDigits (origin: M.Ivanov).                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h> 

#include "TClass.h"
#include "TError.h"

#include "AliLog.h"

#include "AliTRDsegmentID.h"
#include "AliTRDarrayI.h"
#include "AliTRDdataArray.h"

ClassImp(AliTRDdataArray)

//_____________________________________________________________________________
AliTRDdataArray::AliTRDdataArray()
   :fNrow(0)
   ,fNcol(0)
   ,fNtime(0)
   ,fNdim1(-1)
   ,fNdim2(-1)
   ,fIndex(0)
   ,fBufType(-1)
   ,fNelems(-1)
   ,fCurrentIdx1(0)
   ,fCurrentIdx2(0)
   ,fCurrentIndex(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDdataArray::AliTRDdataArray(Int_t nrow, Int_t ncol, Int_t ntime)
   :fNrow(0)
   ,fNcol(0)
   ,fNtime(0)
   ,fNdim1(-1)
   ,fNdim2(-1)
   ,fIndex(0)
   ,fBufType(-1)
   ,fNelems(-1)
   ,fCurrentIdx1(0)
   ,fCurrentIdx2(0)
   ,fCurrentIndex(0)
{
  //
  // Creates a AliTRDdataArray with the dimensions <nrow>, <ncol>, and <ntime>.
  // The row- and column dimensions are compressible.
  //

  Allocate(nrow,ncol,ntime);

}

//_____________________________________________________________________________
AliTRDdataArray::AliTRDdataArray(const AliTRDdataArray &d)
   :AliTRDsegmentID(d)
   ,fNrow(d.fNrow)
   ,fNcol(d.fNcol)
   ,fNtime(d.fNtime)
   ,fNdim1(d.fNdim1)
   ,fNdim2(d.fNdim2)
   ,fIndex(d.fIndex)
   ,fBufType(d.fBufType)
   ,fNelems(d.fNelems)
   ,fCurrentIdx1(0)
   ,fCurrentIdx2(0)
   ,fCurrentIndex(0)
{
  //
  // AliTRDdataArray copy constructor
  //

}

//_____________________________________________________________________________
AliTRDdataArray::~AliTRDdataArray()
{
  //
  // AliTRDdataArray destructor
  //

  if (fIndex) delete fIndex;
  
}

//_____________________________________________________________________________
AliTRDdataArray &AliTRDdataArray::operator=(const AliTRDdataArray &d)
{
  //
  // Assignment operator
  //

  if (this != &d) ((AliTRDdataArray &) d).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDdataArray::Copy(TObject &d) const
{
  //
  // Copy function
  //

  ((AliTRDdataArray &) d).fNrow         = fNrow;
  ((AliTRDdataArray &) d).fNcol         = fNcol;
  ((AliTRDdataArray &) d).fNtime        = fNtime;

  ((AliTRDdataArray &) d).fNdim1        = fNdim1;
  ((AliTRDdataArray &) d).fNdim2        = fNdim2;

  ((AliTRDdataArray &) d).fBufType      = fBufType;
  ((AliTRDdataArray &) d).fNelems       = fNelems;

  ((AliTRDdataArray &) d).fCurrentIdx1  = 0;
  ((AliTRDdataArray &) d).fCurrentIdx2  = 0;
  ((AliTRDdataArray &) d).fCurrentIndex = 0;

  fIndex->Copy(*((AliTRDdataArray &) d).fIndex);

}

//_____________________________________________________________________________
void AliTRDdataArray::Allocate(Int_t nrow, Int_t ncol, Int_t ntime)
{
  //
  // Allocates memory for a AliTRDdataArray with the dimensions 
  // <nrow>, <ncol>, and <ntime>.
  // The row- and column dimensions are compressible.
  //

  if (nrow  <= 0) {
    AliError("The number of rows has to be positive");
    exit(1);
  }
  if (ncol  <= 0) {
    AliError("The number of columns has to be positive");
    exit(1);
  }
  if (ntime <= 0) {
    AliError("The number of timebins has to be positive");
    exit(1);
  }

  // The two-dimensional array row/column gets mapped into the first 
  // dimension of the array. The second array dimension, which is not compressible,
  // corresponds to the time direction
  fNdim1  = nrow * ncol;
  fNdim2  = ntime;
  fNelems = fNdim1 * fNdim2;

  fNrow   = nrow;
  fNcol   = ncol;
  fNtime  = ntime;

  if (fIndex) delete fIndex;
  fIndex = new AliTRDarrayI();
  fIndex->Set(fNdim2);
  for (Int_t i = 0, k = 0; i < fNdim2; i++, k += fNdim1) { 
    (*fIndex)[i] = k;
  }

  fBufType = 0;

}

//_____________________________________________________________________________
Bool_t AliTRDdataArray::CheckBounds(Int_t idx1, Int_t idx2) 
{
  //
  // Does the boundary checking
  //

  if ((idx2 >= fNdim2) || (idx2 < 0)) 
    return OutOfBoundsError(idx1,idx2);

  Int_t index = (*fIndex).At(idx2) + idx1;
  if ((index < 0) || (index > fNelems)) 
    return OutOfBoundsError(idx1,idx2);

  return kTRUE;  

}

//_____________________________________________________________________________
Bool_t AliTRDdataArray::OutOfBoundsError(Int_t idx1, Int_t idx2) 
{
  //
  // Generate an out-of-bounds error. Always returns false.
  //

  AliError(Form("idx1 %d  idx2 %d out of bounds (size: %d x %d, this: 0x%08x)"
	       ,idx1,idx2,fNdim1,fNdim2,this));

  return kFALSE;

}

//_____________________________________________________________________________
void AliTRDdataArray::Reset() 
{ 
  //
  // Reset the array (old content gets deleted)
  //

  if (fIndex) delete fIndex;
  fIndex = new AliTRDarrayI();
  fIndex->Set(0); 

  fNdim1   = -1;
  fNdim2   = -1;
  fNelems  = -1; 

  fBufType = -1;

  fNrow    =  0;
  fNcol    =  0;
  fNtime   =  0;

}

//_____________________________________________________________________________
Int_t AliTRDdataArray::GetIdx1(Int_t row, Int_t col) const
{
  //
  // Maps the two-dimensional row/column plane into an one-dimensional array.
  //

  if (row >= fNrow) {
    AliError(Form("row %d out of bounds (size: %d, this: 0x%08x)",row,fNrow,this));
    return -1;
  }  

  if (col >= fNcol) {
    AliError(Form("col %d out of bounds (size: %d, this: 0x%08x)",col,fNcol,this));
    return -1;
  }  

  return row + col * fNrow;

}

//_____________________________________________________________________________
Int_t AliTRDdataArray::GetIndex(Int_t row, Int_t col, Int_t time) const
{
  //
  // Maps the row/column/time into one number
  // 

  if (time > fNtime) {
    AliError(Form("time %d out of bounds (size: %d, this: 0x%08x)",time,fNtime,this));
    return -1;
  }  
  
  return time * fNrow*fNcol + GetIdx1(row,col);

}
 
