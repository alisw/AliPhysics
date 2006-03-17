/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPlaneType.h,v 1.6 2006/03/17 11:34:29 ivana Exp $

/// \ingroup basic
/// \enum AliMpPlaneType
/// Enumeration for refering to bending and non-bending planes.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
 
#ifndef ALI_MP_PLANE_TYPE_H
#define ALI_MP_PLANE_TYPE_H

#include "AliLog.h"
 
#include <TString.h>

enum AliMpPlaneType
{
  kBendingPlane,    ///< bending plane
  kNonBendingPlane  ///< non-bending plane
};

inline 
TString PlaneTypeName(AliMpPlaneType planeType)
{
  switch ( planeType ) {
    case kBendingPlane:    return "bp";  break;
    case kNonBendingPlane: return "nbp"; break;
  }
  
  // Cannot reach this line
  AliFatalGeneral("AliMpPlaneType.h", "Unknown plane type"); 
  return "invalidPlane";
}       

#endif //ALI_MP_PLANE_TYPE_H
