/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPlaneType.h,v 1.4 2005/10/28 15:03:46 ivana Exp $

/// \ingroup basic
/// \enum AliMpPlaneType
/// Enumeration for refering to bending and non-bending planes.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
 
#ifndef ALI_MP_PLANE_TYPE_H
#define ALI_MP_PLANE_TYPE_H
 
enum AliMpPlaneType
{
  kBendingPlane,    ///< bending plane
  kNonBendingPlane  ///< non-bending plane
};

inline 
const char* PlaneTypeName(AliMpPlaneType planeType)
{
  switch ( planeType ) {
    case kBendingPlane:
      return "BendingPlane";
      break;
    case kNonBendingPlane:
      return "NonBendingPlane";
      break;
  }
  return "invalidPlane";
}       

#endif //ALI_MP_PLANE_TYPE_H
