// $Id$
// Category: plane
//
// Class AliMpSectorPosition
// -------------------------
// Class that represents a placed sector.
// Only translation + reflection transformations can
// be applied.
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpSectorPosition.h"

ClassImp(AliMpSectorPosition)

//______________________________________________________________________________
AliMpSectorPosition::AliMpSectorPosition(const AliMpSector* sector,
                                         const TVector2& offset, 
				         const AliMpIntPair& scale) 
  : TObject(),
    fkSector(sector),
    fOffset(offset),
    fScale(scale)
{
//
}

//_____________________________________________________________________________
AliMpSectorPosition::AliMpSectorPosition(const AliMpSectorPosition& right) 
  : TObject(right) {
// 
  Fatal("AliMpSectorPosition", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpSectorPosition::AliMpSectorPosition() 
  : TObject(),
    fkSector(),
    fOffset(),
    fScale()
{
//
}

//_____________________________________________________________________________
AliMpSectorPosition::~AliMpSectorPosition() {
// 
}

// operators

//_____________________________________________________________________________
AliMpSectorPosition& 
AliMpSectorPosition::operator=(const AliMpSectorPosition& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

