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
AliTRDdataArray::~AliTRDdataArray()
{
  //
  // Destructor
  //

  if (fIndex) fIndex->Delete();
  
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

 
