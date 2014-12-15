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
#include <TClass.h>
#include <TObjArray.h>
#include <Riostream.h>

//-----------------------------------------------------------------------------
/// \class AliMUON1DArray
/// This class is simply a wrapper to a TObjArray, offering in addition a
/// control over the replacement policy when you add
/// something to it.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUON1DArray)
/// \endcond

//_____________________________________________________________________________
AliMUON1DArray::AliMUON1DArray(Int_t theSize)
: AliMUONVStore(),
  fArray(0x0)
{
    /// Default ctor

  if (theSize<=0) theSize=16;
        
  fArray = new TObjArray(theSize);
  fArray->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUON1DArray::AliMUON1DArray(const AliMUON1DArray& other)
: AliMUONVStore(),
  fArray(0x0)
{
/// Copy constructor

    AliDebug(1,Form("this=%p copy ctor",this));
  other.CopyTo(*this);
}

//_____________________________________________________________________________
AliMUON1DArray&
AliMUON1DArray::operator=(const AliMUON1DArray& other)
{
/// Assignment operator

  other.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUON1DArray::~AliMUON1DArray()
{
  /// dtor, we're the owner of our internal array.

  AliDebug(1,Form("this=%p",this));
  delete fArray;
}

//_____________________________________________________________________________
Bool_t
AliMUON1DArray::Add(TObject* object)
{
  /// Add an object to this, if its uniqueID is below maxsize
  if (!object) return kFALSE;
  
  Int_t i = (Int_t)object->GetUniqueID();
  if ( i >= fArray->GetSize() ) 
  {
    AliError(Form("Index out of bounds %u (max is %u)",i,fArray->GetSize()));
    return kFALSE;
  }

  Set(object->GetUniqueID(),object,kFALSE);
  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUON1DArray::Clear(Option_t* opt)
{
  /// Reset
  fArray->Clear(opt);
}

//_____________________________________________________________________________
void
AliMUON1DArray::CopyTo(AliMUON1DArray& dest) const
{
/// Make a deep copy

  delete dest.fArray;
  dest.fArray = 0;
  dest.fArray = new TObjArray(fArray->GetSize());
  dest.fArray->SetOwner(kTRUE);
  for ( Int_t i = 0; i < fArray->GetLast(); ++i )
  {
    dest.fArray->AddAt(fArray->At(i)->Clone(),i);
  }
}

//_____________________________________________________________________________
AliMUON1DArray* 
AliMUON1DArray::Create() const 
{
  /// Create an empty clone of this
  return new AliMUON1DArray(fArray->GetSize());
}

//_____________________________________________________________________________
TObject* 
AliMUON1DArray::FindObject(UInt_t i) const
{
  /// Get the object located at index i, if it exists, and if i is correct.

  if ( (Int_t)(i) < fArray->GetSize() )
  {
    return fArray->At(i);
  }
  AliError(Form("Index %d out of bounds (max %d)",i,fArray->GetSize()));
  return 0x0;
}

//_____________________________________________________________________________
TIterator* 
AliMUON1DArray::CreateIterator() const
{
  /// Return an iterator on this
  return fArray->MakeIterator();
}

//_____________________________________________________________________________
Bool_t 
AliMUON1DArray::Set(Int_t i, TObject* object, Bool_t replace)
{
/// Set the object located at i
/// If replace=kFALSE and there's already an object at location i,
/// this method fails and returns kFALSE, otherwise it returns kTRUE

  if ( i >= 0 && i < fArray->GetSize() )
  {
    if (((Int_t)(object->GetUniqueID()))!=i)
    {
      AliError(Form("object's UniqueID is %d, which is different from the expected %d",
                    object->GetUniqueID(),i));
      return kFALSE;
    }
    
    TObject* o = FindObject(i);
    if ( o && !replace )
    {
      AliError(Form("Object %p is already there for i=%d",o,i));
      return kFALSE;
    }
    if ( o && replace ) 
    {
      delete o;
    }
    fArray->AddAt(object,i);
    return kTRUE;
  }
  AliError(Form("Index %d out of bounds (max %d)",i,fArray->GetSize()));
  return kFALSE;
}

//_____________________________________________________________________________
Int_t 
AliMUON1DArray::GetSize() const
{
  /// Return the number of object we hold
  return fArray->GetEntries();
}
