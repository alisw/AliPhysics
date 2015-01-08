/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup basic
/// \enum AliMq::Station12Type
/// Enumeration for refering to a MUON station12.
/// It is defined in a different namespace than StationType
/// in order to prevent from unwanted mixing of both types
/// in switch command.
///
/// \author Ivana Hrivnacova; IPN Orsay
 
#ifndef ALI_MP_STATION_12_TYPE_H
#define ALI_MP_STATION_12_TYPE_H

#include <TString.h>

namespace AliMq {

  enum Station12Type
  {
    kStation1,  ///< station 1
    kStation2,  ///< station 2
    kNotSt12    ///< value for all non sector stations
  };

  /// Return name for given stationType
  TString Station12TypeName(AliMq::Station12Type station12Type);
}

#endif //ALI_MP_SECTOR_STATION_TYPE_H
