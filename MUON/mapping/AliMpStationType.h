/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpStationType.h,v 1.9 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \enum AliMpStationType
/// Enumeration for refering to a MUON station
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay
 
#ifndef ALI_MP_STATION_TYPE_H
#define ALI_MP_STATION_TYPE_H

#include <TString.h>

namespace AliMp {

  enum StationType
  {
    kStation1,           ///< station 1 (quadrants)
    kStation2,           ///< station 2 (quadrants)
    kStation345,         ///< station 3,4,5 (slats)
    kStationTrigger      ///< trigger stations (slats)
  };

  /// Return name for given stationType
  TString StationTypeName(AliMp::StationType stationType);
}

#endif //ALI_MP_STATION_TYPE_H
