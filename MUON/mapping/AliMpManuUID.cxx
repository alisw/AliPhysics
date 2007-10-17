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

#include "AliMpManuUID.h"

/// \class AliMpManuUID
///
/// Unique ID for manus
/// 
/// \author Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMpManuUID)
/// \endcond

//_____________________________________________________________________________
AliMpManuUID::AliMpManuUID()
  : TObject(), fNofChannels(0)
{
    /// default ctor
}

//_____________________________________________________________________________
AliMpManuUID::AliMpManuUID(Int_t detElemId, Int_t manuId, Int_t nofChannels)
: TObject(), fNofChannels(nofChannels)
{
  /// normal ctor
  SetUniqueID(BuildUniqueID(detElemId,manuId));                            
}

//_____________________________________________________________________________
AliMpManuUID::~AliMpManuUID()
{
  /// dtor
}

//_____________________________________________________________________________
UInt_t 
AliMpManuUID::BuildUniqueID(Int_t detElemId, Int_t manuId)
{
  /// Build a unique id from (de,manu) pair
  
  return ( ( detElemId ) | ( manuId << 12 ) );
}

//_____________________________________________________________________________
Int_t AliMpManuUID::DetElemId(UInt_t uniqueID)
{
  /// Return detection element id part of the uniqueID
  return uniqueID & 0xFFF;
}

//_____________________________________________________________________________
Int_t AliMpManuUID::ManuId(UInt_t uniqueID)
{  
  /// Return manuId part of the uniqueID
  return ( uniqueID & 0xFFF000 ) >> 12;
}

