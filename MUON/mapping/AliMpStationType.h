/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpStationType.h,v 1.5 2005/10/28 15:05:04 ivana Exp $

/// \ingroup basic
/// \enum AliMpStationType
/// Enumeration for refering to a MUON station
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
 
#ifndef ALI_MP_STATION_TYPE_H
#define ALI_MP_STATION_TYPE_H
 
enum AliMpStationType
{
  kStationInvalid = -1,///< invalid station
  kStation1 = 0,       ///< station 1 (quadrants)
  kStation2,           ///< station 2 (quadrants)
  kStation345,         ///< station 3,4,5 (slats)
  kStationTrigger      ///< trigger stations (slats)

};

inline 
const char* StationTypeName(AliMpStationType stationType)
{
  switch ( stationType )
  {
    case kStation1:
      return "st1";
      break;
    case kStation2:
      return "st2";
      break;
    case kStation345:
      return "slat";
      break;
    case kStationTrigger:
      return "trigger";
      break;
    case kStationInvalid:
    default:
      return "invalid";
      break;
  }
  return "unknown";
}

#endif //ALI_MP_STATION_TYPE_H
