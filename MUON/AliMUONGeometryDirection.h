/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Enum AliMUONGeometryDirection
// -----------------------------
// Enumeration for the directions in a plane.
// This method gives the direction where the pad size of the chamber is constant.
// In other words this corresponds with the coordinate where the spatial resolution is the best.
// Normally kDirY will correspond with cathode segmentation for the bending plane and 
// kDirX  with cathode segmentation for the non bending plane
//
// Author:Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_DIRECTION_H
#define ALI_MUON_GEOMETRY_DIRECTION_H

enum AliMUONGeometryDirection 
{
  kDirX,          // direction in x
  kDirY,          // direction in y
  kDirUndefined   // undefined direction
};

#endif //ALI_MUON_GEOMETRY_DIRECTION_H
