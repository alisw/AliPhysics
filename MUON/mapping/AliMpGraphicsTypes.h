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

#include <vector>

class AliMpGraphContext;

#ifdef __HP_aCC
  typedef vector<AliMpGraphContext*> GraphContextVector;
#else
  typedef std::vector<AliMpGraphContext*> GraphContextVector;
#endif

#endif //ALI_MP_GRAPHICS_TYPES_H
