// $Id$
// Category: graphics
//
// AliMpGraphicsTypes
// ------------------
// Sytem dependent types definitions for graphics category.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

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
