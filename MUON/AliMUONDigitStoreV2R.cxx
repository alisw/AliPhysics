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
/// \class AliMUONDigitStoreV2R
///
/// Concrete implementation of AliMUONVDigitStore for real digits, using
/// the AliMUONDigitStoreVImpl base implementation
/// 
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONDigitStoreV2R.h"

#include "AliLog.h"
#include "AliMUONRealDigit.h"
#include <TClonesArray.h>

/// \cond CLASSIMP
ClassImp(AliMUONDigitStoreV2R)
/// \endcond

//_____________________________________________________________________________
AliMUONDigitStoreV2R::AliMUONDigitStoreV2R()
: AliMUONDigitStoreVImpl("AliMUONRealDigit")
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONDigitStoreV2R::~AliMUONDigitStoreV2R()
{
  /// dtor
}

//_____________________________________________________________________________
AliMUONVDigit*
AliMUONDigitStoreV2R::AddConcreteDigit(TClonesArray& a, 
                                       const AliMUONVDigit& digit,
                                       Int_t index)
{
  /// Add a digit to this store
  
  const AliMUONRealDigit* d = dynamic_cast<const AliMUONRealDigit*>(&digit);
  
  if ( !d ) 
  {
    AliError(Form("Digit (of class %s) is not of the expected type AliMUONRealDigit",
                  digit.ClassName()));
    return 0x0;
  }
  
  return new(a[index]) AliMUONRealDigit(*d);
}

//_____________________________________________________________________________
AliMUONVDigitStore* 
AliMUONDigitStoreV2R::Create() const
{
  /// Create an empty store
  return new AliMUONDigitStoreV2R;
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreV2R::CreateDigit(Int_t detElemId, Int_t manuId,                                 
                                  Int_t manuChannel, Int_t cathode) const
{
  /// Create a digit that is compatible with this store
  return new AliMUONRealDigit(detElemId,manuId,manuChannel,cathode);
}


