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

#include "AliMUONRealDigit.h"

//-----------------------------------------------------------------------------
/// \class AliMUONRealDigit
///
/// Implementation of AliMUONVDigit for real digit.
/// 
/// This class should store the bare minimum in order to save disk space
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONRealDigit)
/// \endcond

//_____________________________________________________________________________
AliMUONRealDigit::AliMUONRealDigit() 
  : AliMUONVDigit(),
  fCharge(0),
  fPadXY(0),
  fADC(0),
  fStatusMap(0)
{
      /// default ctor
}

//_____________________________________________________________________________
AliMUONRealDigit::AliMUONRealDigit(Int_t detElemId, Int_t manuId, 
                                   Int_t manuChannel, Int_t cathode)
: AliMUONVDigit(detElemId,manuId,manuChannel,cathode),
fCharge(0),
fPadXY(0),
fADC(0),
fStatusMap(0)
{
  /// normal ctor
}

//_____________________________________________________________________________
AliMUONRealDigit::~AliMUONRealDigit()
{
  /// empty ctor
}

//_____________________________________________________________________________
Bool_t
AliMUONRealDigit::MergeWith(const AliMUONVDigit& src)
{
  /// Merge with src.
  
  Bool_t check = ( src.DetElemId() == DetElemId() &&
                   src.PadX() == PadX() &&
                   src.PadY() == PadY() &&
                   src.Cathode() == Cathode() );
  if (!check)
  {
    return kFALSE;
  }
  
  AddCharge(src.Charge());
  return kTRUE;
}

//_____________________________________________________________________________
Int_t
AliMUONRealDigit::PadX() const
{
  /// Return (integer) position in x (within the detection element)
  return fPadXY & 0xFFFF;
}

//_____________________________________________________________________________
Int_t
AliMUONRealDigit::PadY() const
{
  /// Return (integer) position in y (within the detection element)
  return ( fPadXY & 0xFFFF0000 ) >> 16;
}

//_____________________________________________________________________________
void
AliMUONRealDigit::SetPadXY(Int_t padx, Int_t pady)
{
  /// Set the pad (integer) positions
  fPadXY = ( padx | (pady << 16) );
}

