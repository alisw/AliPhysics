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

// $Id$

#include "AliMUON1DArray.h"

#include "AliLog.h"
#include "TObjArray.h"

///
/// This class is simply a wrapper to a TObjArray, offering in addition a
/// control over the replacement policy when you add
/// something to it.
///

ClassImp(AliMUON1DArray)

//_____________________________________________________________________________
AliMUON1DArray::AliMUON1DArray(Int_t theSize)
: AliMUONV1DStore(),
  fArray(0x0)
{
  // 
  // Default ctor
  //
  if ( theSize ) 
  {
    fArray = new TObjArray(theSize);
  }
}

//_____________________________________________________________________________
AliMUON1DArray::AliMUON1DArray(const AliMUON1DArray& other)
: AliMUONV1DStore(),
  fArray(0x0)
{
  other.CopyTo(*this);
}

//_____________________________________________________________________________
AliMUON1DArray&
AliMUON1DArray::operator=(const AliMUON1DArray& other)
{
  other.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUON1DArray::~AliMUON1DArray()
{
  //
  // dtor, we're the owner of our internal array.
  //
  delete fArray;
}

//_____________________________________________________________________________
void
AliMUON1DArray::CopyTo(AliMUON1DArray& dest) const
{
  //
  // Make a deep copy
  //
  delete dest.fArray;
  dest.fArray = new TObjArray;
  dest.fArray->SetOwner(kTRUE);
  for ( Int_t i = 0; i < fArray->GetLast(); ++i )
  {
    dest.fArray->AddAt(fArray->At(i)->Clone(),i);
  }
}

//_____________________________________________________________________________
TObject* 
AliMUON1DArray::Get(Int_t i) const
{
  //
  // Get the object located at index i, if it exists, and if i is correct.
  //
  if ( i >= 0 && i < fArray->GetSize() )
  {
    return fArray->At(i);
  }
  AliError(Form("Index %d out of bounds (max %d)",i,fArray->GetSize()));
  return 0x0;
}

//_____________________________________________________________________________
Bool_t 
AliMUON1DArray::Set(Int_t i, TObject* object, Bool_t replace)
{
  //
  // Set the object located at i
  // If replace=kFALSE and there's already an object at location i,
  // this method fails and returns kFALSE, otherwise it returns kTRUE
  //
  if ( i >= 0 && i < fArray->GetSize() )
  {
    TObject* o = Get(i);
    if ( o && !replace )
    {
      AliError(Form("Object %p is already there for i=%d",o,i));
      return kFALSE;
    }
    if ( replace ) 
    {
      delete o;
    }
    fArray->AddAt(object,i);
    return kTRUE;
  }
  AliError(Form("Index %d out of bounds (max %d)",i,fArray->GetSize()));
  return kFALSE;
}



