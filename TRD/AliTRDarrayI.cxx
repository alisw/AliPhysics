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
Revision 1.5  2000/11/01 14:53:20  cblume
Merge with TRD-develop


Revision 1.1.4.3  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.4.2  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.4  2000/06/09 11:10:07  cblume
Compiler warnings and coding conventions, next round

Revision 1.3  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 14:35:54  cblume
Update

Revision 1.4  2000/06/09 11:10:07  cblume
Compiler warnings and coding conventions, next round

Revision 1.3  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.4.1  2000/05/08 14:35:54  cblume
Update

Revision 1.1  2000/02/28 18:57:18  cblume
Add new TRD classes

*/

///////////////////////////////////////////////////////////////////////
//                                                                   //  
// Added additional functionality to the original TArrayI.           //
//  - Multiple inheritance from TObject                              //
//  - Function Expand() allows to expand the array without           //
//    deleting the array contents                                    //
//                                                                   //
// Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk  // 
//                                                                   //  
///////////////////////////////////////////////////////////////////////

#include "AliTRDarrayI.h"

ClassImp(AliTRDarrayI)

//_____________________________________________________________________________
AliTRDarrayI::AliTRDarrayI():TArrayI()
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDarrayI::~AliTRDarrayI()
{
  //
  // Default destructor
  //

}

//_____________________________________________________________________________
void AliTRDarrayI::Copy(TObject &a)
{
  //
  // Copy function
  //

  TObject::Copy(a);
  TArrayI::Copy(((TArrayI &) a));

}

//_____________________________________________________________________________
void AliTRDarrayI::Expand(Int_t n)
{
  //
  // Sets the  array size of the TArrayI object to <n> integers and copies
  // the old array.
  // If n < 0 leave the array unchanged.
  // The user is responsible for the appropriate size of the array.
  //

  if (n < 0) return;  
  fArray = (Int_t*) TStorage::ReAlloc(fArray
                                     ,n  * sizeof(Int_t)
                                     ,fN * sizeof(Int_t));
  if (fArray != 0) fN = n;
 
}

