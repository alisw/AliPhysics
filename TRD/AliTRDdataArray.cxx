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

/*
$Log$
Revision 1.4  2000/06/07 16:27:01  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.3  2000/05/18 07:56:44  cblume
Added #include <stdlib.h>

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 15:13:59  cblume
Introduce boundary checking

Revision 1.1  2000/02/28 18:59:19  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Base class of a general container for data of a TRD detector segment.    //
//  Adapted from AliDigits (origin: M.Ivanov).                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h> 

#include "TClass.h"
#include "TError.h"
#include "AliTRDsegmentID.h"
#include "AliTRDarrayI.h"
#include "AliTRDdataArray.h"

ClassImp(AliTRDdataArray)

//_____________________________________________________________________________
AliTRDdataArray::AliTRDdataArray()
{
  //
  // Default constructor
  //

  fIndex   =  0;

  fNdim1   = -1;
  fNdim2   = -1;
  fNelems  = -1; 

  fBufType = -1;

  fNrow    =  0;
  fNcol    =  0;
  fNtime   =  0;

}

//_____________________________________________________________________________
AliTRDdataArray::AliTRDdataArray(Int_t nrow, Int_t ncol, Int_t ntime)
{
  //
  // Creates a AliTRDdataArray with the dimensions <nrow>, <ncol>, and <ntime>.
  // The row- and column dimensions are compressible.
  //

  Allocate(nrow,ncol,ntime);

}

//_____________________________________________________________________________
AliTRDdataArray::AliTRDdataArray(AliTRDdataArray &d)
{
  //
  // AliTRDdataArray copy constructor
  //

  d.Copy(*this);

}

//_____________________________________________________________________________
AliTRDdataArray::~AliTRDdataArray()
{
  //
  // AliTRDdataArray destructor
  //

  if (fIndex) fIndex->Delete();
  
}

//_____________________________________________________________________________
void AliTRDdataArray::Copy(AliTRDdataArray &d)
{
  //
  // Copy function
  //

  d.fNrow         = fNrow;
  d.fNcol         = fNcol;
  d.fNtime        = fNtime;

  d.fNdim1        = fNdim1;
  d.fNdim2        = fNdim2;

  d.fBufType      = fBufType;
  d.fNelems       = fNelems;

  d.fCurrentIdx1  = 0;
  d.fCurrentIdx2  = 0;
  d.fCurrentIndex = 0;

  fIndex->Copy(*d.fIndex);

}

//_____________________________________________________________________________
void AliTRDdataArray::Allocate(Int_t nrow, Int_t ncol,Int_t ntime)
{
  //
  // Allocates memory for a AliTRDdataArray with the dimensions 
  // <nrow>, <ncol>, and <ntime>.
  // The row- and column dimensions are compressible.
  //

  if (nrow  <= 0) {
    Error("AliTRDdataArray::Allocate","The number of rows has to be positive");
    exit(1);
  }
  if (ncol  <= 0) {
    Error("AliTRDdataArray::Allocate","The number of columns has to be positive");
    exit(1);
  }
  if (ntime <= 0) {
    Error("AliTRDdataArray::Allocate","The number of timebins has to be positive");
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
  fIndex = new AliTRDarrayI;
  fIndex->Set(fNdim2);
  for (Int_t i = 0, k = 0; i < fNdim2; i++, k += fNdim1) { 
    (*fIndex)[i] = k;
  }

  fBufType = 0;

}

//_____________________________________________________________________________
void AliTRDdataArray::Reset() 
{ 
  //
  // Reset the array (old content gets deleted)
  //

  if (fIndex) delete fIndex;
  fIndex = new AliTRDarrayI;
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
Int_t AliTRDdataArray::GetIdx1(Int_t row, Int_t col)
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
Int_t AliTRDdataArray::GetIndex(Int_t row, Int_t col, Int_t time)
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
 
