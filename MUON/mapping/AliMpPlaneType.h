/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPlaneType.h,v 1.3 2005/08/26 15:43:36 ivana Exp $

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

#endif //ALI_MP_PLANE_TYPE_H
