// $Id$
// Category: plane
//
// AliMpPlaneTypes
// ---------------
// Sytem dependent types definitions for plane category.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PLANE_TYPES_H
#define ALI_MP_PLANE_TYPES_H

#include "AliMpContainers.h"

#ifdef WITH_STL
  #include <vector>
#endif

#ifdef WITH_ROOT
  #include <TObjArray.h>
#endif

class AliMpSectorPosition;
class AliMpTransformPadIterator;
class AliMpTransformer;

#ifdef WITH_STL
#ifdef __HP_aCC
  typedef vector<AliMpSectorPosition*>  SectorPositionVector;
  typedef vector<AliMpTransformPadIterator*>  PadIteratorVector;
  typedef PadIteratorVector::iterator PadIteratorVectorIterator;
  typedef vector<AliMpTransformer*>  TransformerVector;
#else
  typedef std::vector<AliMpSectorPosition*>  SectorPositionVector;
  typedef std::vector<AliMpTransformPadIterator*>  PadIteratorVector;
  typedef PadIteratorVector::iterator PadIteratorVectorIterator;
  typedef std::vector<AliMpTransformer*>  TransformerVector;
#endif
#endif

#ifdef WITH_ROOT
  typedef  TObjArray  SectorPositionVector;
  typedef  TObjArray  PadIteratorVector;
  typedef  Int_t      PadIteratorVectorIterator;
  typedef  TObjArray  TransformerVector;
#endif


#endif //ALI_MP_PLANE_TYPES_H
