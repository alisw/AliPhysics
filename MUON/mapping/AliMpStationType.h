/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpStationType.h,v 1.6 2006/01/11 10:07:18 ivana Exp $

/// \ingroup basic
/// \enum AliMpStationType
/// Enumeration for refering to a MUON station
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
 
#ifndef ALI_MP_STATION_TYPE_H
#define ALI_MP_STATION_TYPE_H

#include <TString.h>

#include "AliLog.h"
 
enum AliMpStationType
{
  kStation1,           ///< station 1 (quadrants)
  kStation2,           ///< station 2 (quadrants)
  kStation345,         ///< station 3,4,5 (slats)
  kStationTrigger      ///< trigger stations (slats)
};

inline 
TString StationTypeName(AliMpStationType stationType)
{
  switch ( stationType ) {
    case kStation1:       return "st1";     break;
    case kStation2:       return "st2";     break;
    case kStation345:     return "slat";    break;
    case kStationTrigger: return "trigger"; break;
  }
  
  // Cannot reach this line
  AliFatalGeneral("AliMpStationType.h", "Unknown station type"); 
  return "invalidStation";
}

#endif //ALI_MP_STATION_TYPE_H
