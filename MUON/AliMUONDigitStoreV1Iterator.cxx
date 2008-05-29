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

//-----------------------------------------------------------------------------
/// \class AliMUONDigitStoreV1Iterator
///
/// Implementation of TIteraor for AliMUONVDigitStoreV1
/// Reuses the AliMUONTOTCAStoreIterator iterator
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONDigitStoreV1Iterator.h"

#include "AliLog.h"
#include "AliMpDEManager.h"
#include "AliMUONVDigit.h"
#include "TObjArray.h"

/// \cond CLASSIMP
ClassImp(AliMUONDigitStoreV1Iterator)
/// \endcond

//_____________________________________________________________________________
AliMUONDigitStoreV1Iterator::AliMUONDigitStoreV1Iterator(TObjArray* a,
                                                         Int_t firstDetElemId,
                                                         Int_t lastDetElemId,
                                                         Int_t cathode)
: AliMUONTOTCAStoreIterator(a,AliMpDEManager::GetChamberId(firstDetElemId),
                            AliMpDEManager::GetChamberId(lastDetElemId)),
fArray(a),
fFirstDetElemId(firstDetElemId),
fLastDetElemId(lastDetElemId),
fCathode(cathode)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONDigitStoreV1Iterator::AliMUONDigitStoreV1Iterator(const AliMUONDigitStoreV1Iterator& rhs)
: AliMUONTOTCAStoreIterator(rhs),
  fArray(rhs.fArray),
  fFirstDetElemId(rhs.fFirstDetElemId),
  fLastDetElemId(rhs.fLastDetElemId),
  fCathode(rhs.fCathode)
{
    /// copy ctor
}

//_____________________________________________________________________________
AliMUONDigitStoreV1Iterator& 
AliMUONDigitStoreV1Iterator::operator=(const TIterator& rhs)
{
  /// overriden assignment operator (imposed by Root's definition of TIterator ?)
  
  if ( this != &rhs )
  {
    if ( rhs.IsA() != AliMUONDigitStoreV1Iterator::Class() )
    {
      AliErrorGeneral("AliMUONDigitStoreV1Iterator::operator=","Wrong type");
    }
    else
    {
      const AliMUONDigitStoreV1Iterator& rhs1 = 
      static_cast<const AliMUONDigitStoreV1Iterator&>(rhs);
      
      AliMUONDigitStoreV1Iterator::operator=(rhs1);
    }
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONDigitStoreV1Iterator& 
AliMUONDigitStoreV1Iterator::operator=(const AliMUONDigitStoreV1Iterator& rhs)
{
  /// assignement operator
  if ( this != &rhs ) 
  {
    TIterator::operator=(rhs);
    fArray = rhs.fArray;
    fFirstDetElemId = rhs.fFirstDetElemId;
    fLastDetElemId = rhs.fLastDetElemId;
    fCathode = rhs.fCathode;
  }
  return *this;
}

//_____________________________________________________________________________
AliMUONDigitStoreV1Iterator::~AliMUONDigitStoreV1Iterator()
{
  /// dtor
}

//_____________________________________________________________________________
const TCollection* 
AliMUONDigitStoreV1Iterator::GetCollection() const
{
  /// Return the TObjArray we're iterating upon
  return fArray;
}

//_____________________________________________________________________________
TObject*
AliMUONDigitStoreV1Iterator::Next()
{
  /// Return the next digit (with its DE in [fFirstDetElemId,fLastDetElemId],
  /// and its cathode == fCathode (or any cathode if fCathode==2)
  /// in the store.
  
  TObject* object = 0x0;
  
  while ( (object = static_cast<AliMUONVDigit*>(AliMUONTOTCAStoreIterator::Next()) ) )
  {  
    AliMUONVDigit* digit = static_cast<AliMUONVDigit*>(object);

    if ( digit->DetElemId() >= fFirstDetElemId &&
         digit->DetElemId() <= fLastDetElemId ) 
    {
      if ( fCathode == 2 || digit->Cathode() == fCathode ) 
      {
        return digit;
      }
    }
  }
  
  return 0x0;
}
