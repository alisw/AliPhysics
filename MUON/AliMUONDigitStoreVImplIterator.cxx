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
/// \class AliMUONDigitStoreVImplIterator
///
/// Implementation of AliMUONVDataIterator for AliMUONDigitStoreVImpl
///
/// \author Laurent Aphecetche, Subatech
///
//-----------------------------------------------------------------------------

#include "AliMUONDigitStoreVImplIterator.h"

#include "AliMUONVDigit.h"
#include "AliMUONDigitStoreVImpl.h"
#include "AliMUON2DMap.h"
#include "AliMUONVCalibParam.h"
#include <TClonesArray.h>
#include <TError.h>

/// \cond CLASSIMP
ClassImp(AliMUONDigitStoreVImplIterator)
/// \endcond

//_____________________________________________________________________________
AliMUONDigitStoreVImplIterator::AliMUONDigitStoreVImplIterator(const AliMUONDigitStoreVImpl* store)
: TIterator(),
fkStore(store),
fFirstDetElemId(100),
fLastDetElemId(1417),
fCathode(2),
fStoreIterator(store->fMap->CreateIterator()),
fCurrentCalibParam(0x0),
fCurrentCalibParamIndex(-1)
{
  /// ctor for full iteration
}

//_____________________________________________________________________________
AliMUONDigitStoreVImplIterator::AliMUONDigitStoreVImplIterator(const AliMUONDigitStoreVImpl* store,
                                                               Int_t firstDE,
                                                               Int_t lastDE,
                                                               Int_t cathode)
: TIterator(),
fkStore(store),
fFirstDetElemId(firstDE),
fLastDetElemId(lastDE),
fCathode(cathode),
fStoreIterator(store->fMap->CreateIterator(firstDE,lastDE)),
fCurrentCalibParam(0x0),
fCurrentCalibParamIndex(-1)
{
  /// ctor for partial iteration
}

//_____________________________________________________________________________
AliMUONDigitStoreVImplIterator&
AliMUONDigitStoreVImplIterator::operator=(const TIterator&)
{
  // overriden assignment operator (imposed by Root's declaration of Titerator ?)
  Fatal("TIterator::operator=","Not implementeable"); // because there's no clone in TIterator :-(
  return *this;
}

//_____________________________________________________________________________
AliMUONDigitStoreVImplIterator::~AliMUONDigitStoreVImplIterator()
{
  /// dtor
  delete fStoreIterator;
}

//_____________________________________________________________________________
TObject*
AliMUONDigitStoreVImplIterator::Next()
{
  /// Return next digit in store
  if ( !fCurrentCalibParam ) 
  {
    fCurrentCalibParam = static_cast<AliMUONVCalibParam*>(fStoreIterator->Next());
    fCurrentCalibParamIndex = 0;
    if ( !fCurrentCalibParam ) return 0x0;
  }
  
  Int_t ix(-1);
  AliMUONVDigit* d(0x0);
  
  if ( fCathode == 2 ) 
  {
    while ( fCurrentCalibParamIndex < 64 && ix < 0 )
    {
      ix = fCurrentCalibParam->ValueAsInt(fCurrentCalibParamIndex++);
    };
    
    if (ix>=0)
    {
      d = static_cast<AliMUONVDigit*>(fkStore->fDigits->UncheckedAt(ix));
    }  
  }
  else
  {
    while ( d == 0x0 ) 
    {
      while ( fCurrentCalibParamIndex < 64 && ix < 0 )
      {
        ix = fCurrentCalibParam->ValueAsInt(fCurrentCalibParamIndex++);
      };
    
      if (ix>=0)
      {
        d = static_cast<AliMUONVDigit*>(fkStore->fDigits->UncheckedAt(ix));
        
        if (  fCathode == 2 || d->Cathode() == fCathode ) 
        {
          break;
        }
        d = 0x0;
        ix = -1;
      }
      else
      {
        break;
      }
    }
  }
  
  if (ix<0) 
  {
    fCurrentCalibParam = 0x0;
    return Next();
  }
  
  return d;
}

//_____________________________________________________________________________
void
AliMUONDigitStoreVImplIterator::Reset()
{
  /// Reset the iterator
  fCurrentCalibParam = 0x0;
  fCurrentCalibParamIndex = 0;  
  fStoreIterator->Reset();
}
