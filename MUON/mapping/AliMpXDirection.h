/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpXDirection.h,v 1.4 2005/08/26 15:43:36 ivana Exp $

/// \ingroup basic
/// \enum  AliMpDirection
/// Enumeration for the directions in a plane.
///
/// Authors: David Guez, IPN Orsay

#ifndef ALI_MP_X_DIRECTION_H
#define ALI_MP_X_DIRECTION_H

enum AliMpXDirection
{
  kLeft,  ///< special row segments built from right to left
  kRight  ///< special row segments built from left to right
};

#endif //ALI_MP_X_DIRECTION_H
