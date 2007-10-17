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

#include "AliMpPadUID.h"

/// \class AliMpPadUID
///
/// Unique ID for pads
/// 
/// \author Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMpPadUID)
/// \endcond

//_____________________________________________________________________________
AliMpPadUID::AliMpPadUID(UInt_t uid)
: TObject()
{
  /// ctor
  SetUniqueID(uid);
}

//_____________________________________________________________________________
AliMpPadUID::AliMpPadUID(Int_t detElemId, Int_t manuId, Int_t manuChannel)
: TObject()
{
  /// ctor
  SetUniqueID(BuildUniqueID(detElemId,manuId,manuChannel));
}

//_____________________________________________________________________________
AliMpPadUID::~AliMpPadUID()
{
  /// dtor
}

//_____________________________________________________________________________
UInt_t 
AliMpPadUID::BuildUniqueID(Int_t detElemId, Int_t manuId, 
                          Int_t manuChannel)
{
  /// Build a single integer with id information
  return ( ( detElemId ) | ( manuId << 12 ) | ( manuChannel << 24 ) );
}

//_____________________________________________________________________________
Int_t
AliMpPadUID::DetElemId(UInt_t uniqueID)
{
  /// Return detection element id part of the uniqueID
  return uniqueID & 0xFFF;
}

//_____________________________________________________________________________
Int_t
AliMpPadUID::ManuChannel(UInt_t uniqueID)
{
  /// Return manuChannel part of the uniqueID
  return ( uniqueID & 0x3F000000 ) >> 24;
}

//_____________________________________________________________________________
Int_t
AliMpPadUID::ManuId(UInt_t uniqueID)
{
  /// Return manuId part of the uniqueID
  return ( uniqueID & 0xFFF000 ) >> 12;
}
