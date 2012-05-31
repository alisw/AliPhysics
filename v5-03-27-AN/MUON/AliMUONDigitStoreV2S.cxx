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

#include "AliMUONDigitStoreV2S.h"

//-----------------------------------------------------------------------------
/// \class AliMUONDigitStoreV2S
///
/// Concrete implementation of AliMUONVDigitStore for simulated digits, using
/// the AliMUONDigitStoreVImpl base implementation
/// 
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------


#include "AliMUONDigit.h"
#include <TClonesArray.h>

/// \cond CLASSIMP
ClassImp(AliMUONDigitStoreV2S)
/// \endcond

//_____________________________________________________________________________
AliMUONDigitStoreV2S::AliMUONDigitStoreV2S()
: AliMUONDigitStoreVImpl("AliMUONDigit")
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONDigitStoreV2S::~AliMUONDigitStoreV2S()
{
  /// dtor
}

//_____________________________________________________________________________
AliMUONVDigit*
AliMUONDigitStoreV2S::AddConcreteDigit(TClonesArray& a, 
                                       const AliMUONVDigit& digit,
                                       Int_t index)
{
  /// add a digit to this store  
  
  if ( digit.IsA() != AliMUONDigit::Class() ) 
  {
    AliMUONDigit d(digit.DetElemId(),digit.ManuId(),digit.ManuChannel(),digit.Cathode());
    d.SetCharge(digit.Charge());
    d.SetADC(digit.ADC());
    d.SetPadXY(digit.PadX(),digit.PadY());
    d.ChargeInFC();
    d.Converted();
    return new(a[index]) AliMUONDigit(d);
  }

  return new(a[index]) AliMUONDigit(static_cast<const AliMUONDigit&>(digit));
}

//_____________________________________________________________________________
AliMUONVDigitStore* 
AliMUONDigitStoreV2S::Create() const
{
  /// create an empty store
  return new AliMUONDigitStoreV2S;
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitStoreV2S::CreateDigit(Int_t detElemId, Int_t manuId,                                 
                                  Int_t manuChannel, Int_t cathode) const
{
  /// create a digit, compatible with this store
  return new AliMUONDigit(detElemId,manuId,manuChannel,cathode);
}


