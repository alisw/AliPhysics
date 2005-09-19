/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpBasicTypes.h,v 1.4 2005/08/26 15:43:36 ivana Exp $

/// \ingroup basic
///
/// AliMpBasicTypes
/// System dependent types definitions for basic category.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_BASIC_TYPES_H
#define ALI_MP_BASIC_TYPES_H

#include <utility>
#include "AliMpPad.h"

#ifdef __HP_aCC
  typedef pair<AliMpPad, AliMpPad> PadPair;
#else
  typedef std::pair<AliMpPad, AliMpPad> PadPair;
#endif

#endif //ALI_MP_BASIC_TYPES_H
