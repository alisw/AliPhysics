// $Id$
// Category: basic
//
// AliMpBasicTypes
// ---------------
// Sytem dependent types definitions for basic category.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

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
