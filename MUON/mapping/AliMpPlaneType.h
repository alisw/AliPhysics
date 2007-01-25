/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPlaneType.h,v 1.8 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \enum AliMpPlaneType
/// Enumeration for refering to bending and non-bending planes.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay
 
#ifndef ALI_MP_PLANE_TYPE_H
#define ALI_MP_PLANE_TYPE_H

#include <TString.h>

namespace AliMp {

  enum PlaneType
  {
    kBendingPlane,    ///< bending plane
    kNonBendingPlane  ///< non-bending plane
  };

  /// Return name for given planeType
  TString PlaneTypeName(AliMp::PlaneType planeType);


  /// Return the other plane type
  AliMp::PlaneType OtherPlaneType(AliMp::PlaneType planeType);
}  

#endif //ALI_MP_PLANE_TYPE_H
