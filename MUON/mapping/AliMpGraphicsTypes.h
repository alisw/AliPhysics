/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpGraphicsTypes.h,v 1.5 2005/08/26 15:43:36 ivana Exp $

/// \ingroup graphics
/// AliMpGraphicsTypes
/// System dependent types definitions for graphics category.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_GRAPHICS_TYPES_H
#define ALI_MP_GRAPHICS_TYPES_H

#include "AliMpContainers.h"

#ifdef WITH_STL
  #include <vector>
#endif

#ifdef WITH_ROOT
  #include <TObjArray.h>
#endif

class AliMpGraphContext;

#ifdef WITH_STL
#ifdef __HP_aCC
  typedef vector<AliMpGraphContext*> GraphContextVector;
#else
  typedef std::vector<AliMpGraphContext*> GraphContextVector;
#endif
#endif

#ifdef WITH_ROOT
  typedef TObjArray GraphContextVector;
#endif

#endif //ALI_MP_GRAPHICS_TYPES_H
