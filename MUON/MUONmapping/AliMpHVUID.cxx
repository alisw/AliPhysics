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

#include "AliMpHVUID.h"

/// \class AliMpHVUID
/// 
/// A utility class to assign a unique ID to a given HV channel
///
/// \author: Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMpHVUID)
/// \endcond

//_____________________________________________________________________________
AliMpHVUID::AliMpHVUID() : TObject()
{
  /// ctor
}

//_____________________________________________________________________________
AliMpHVUID::~AliMpHVUID()
{
  /// dtor
}

//_____________________________________________________________________________
UInt_t
AliMpHVUID::BuildUniqueID(Int_t detElemId, Int_t index)
{
  /// Build a single index from the pair (de,index)
  return ( index | ( detElemId << 16 ) );
}

//_____________________________________________________________________________
Int_t
AliMpHVUID::Index(UInt_t uniqueID)
{
  /// Extract index from uniqueID
  return uniqueID & 0xFFFF;
}

//_____________________________________________________________________________
Int_t
AliMpHVUID::DetElemId(UInt_t uniqueID)
{
  /// Extract detElemId from uniqueID
  return ( uniqueID & 0xFFFF0000 ) >> 16;
}

