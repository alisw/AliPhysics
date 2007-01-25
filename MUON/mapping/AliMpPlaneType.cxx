/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPlaneType.h,v 1.8 2006/05/24 13:58:07 ivana Exp $

// \enum AliMpPlaneType
// Enumeration for refering to bending and non-bending planes.
//
// Author: David Guez, Ivana Hrivnacova; IPN Orsay
 
#include "AliMpPlaneType.h"

#include "AliLog.h" 
 
//_____________________________________________________________________________
TString AliMp::PlaneTypeName(PlaneType planeType)
{
  switch ( planeType ) {
    case kBendingPlane:    return "bp";  break;
    case kNonBendingPlane: return "nbp"; break;
  }
  
  // Cannot reach this line
  AliFatalGeneral("AliMpPlaneType.h", "Unknown plane type"); 
  return "invalidPlane";
}       

//_____________________________________________________________________________
AliMp::PlaneType AliMp::OtherPlaneType(PlaneType planeType)
{
  switch ( planeType ) {
    case kBendingPlane:    return kNonBendingPlane;  break;
    case kNonBendingPlane: return kBendingPlane;     break;
  }
  
  // Cannot reach this line
  AliFatalGeneral("AliMpPlaneType.h", "Unknown plane type"); 
  return kBendingPlane;
}       
