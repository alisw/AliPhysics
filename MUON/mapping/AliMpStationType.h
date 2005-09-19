/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpStationType.h,v 1.4 2005/08/26 15:43:36 ivana Exp $

/// \ingroup basic
/// \enum AliMpStationType
/// Enumeration for refering to a MUON station
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
 
#ifndef ALI_MP_STATION_TYPE_H
#define ALI_MP_STATION_TYPE_H
 
enum AliMpStationType
{
  kStation1,  ///< station 1 (quadrants)
  kStation2,  ///< station 2 (quadrants)
  kStation345 ///< station 3,4,5 (slats)
};

#endif //ALI_MP_STATION_TYPE_H
