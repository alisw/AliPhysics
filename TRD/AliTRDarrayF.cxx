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

Revision 1.1.2.1  2000/05/08 14:35:38  cblume
Add float array

*/

///////////////////////////////////////////////////////////////////////
//                                                                   //  
// Added additional functionality to the original TArrayF.           //
//  - Multiple inheritance from TObject                              //
//  - Function Expand() allows to expand the array without           //
//    deleting the array contents                                    //
//                                                                   //
// Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk  // 
//                                                                   //  
///////////////////////////////////////////////////////////////////////

#include "AliTRDarrayF.h"

ClassImp(AliTRDarrayF)

//_____________________________________________________________________________
AliTRDarrayF::~AliTRDarrayF()
{
  //
  // Default destructor
  //

}

//_____________________________________________________________________________
void AliTRDarrayF::Copy(AliTRDarrayF &a)
{
  //
  // Copy function
  //

  TObject::Copy(a);
  TArrayF::Copy(a);

}

//_____________________________________________________________________________
void AliTRDarrayF::Expand(Int_t n)
{
  //
  // Sets the  array size of the TArrayF object to <n> integers and copies
  // the old array.
  // If n < 0 leave the array unchanged.
  // The user is responsible for the appropriate size of the array.
  //

  if (n < 0) return;  
  fArray = (Float_t*) TStorage::ReAlloc(fArray
                                       ,n  * sizeof(Float_t)
                                       ,fN * sizeof(Float_t));
  if (fArray != 0) fN = n;
 
}

