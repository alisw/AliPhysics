// $Id$
// Category: sector
//
// Class AliMpVRowSegment
// ----------------------
// Class describing an interface for a row segment.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpVRowSegment.h"

ClassImp(AliMpVRowSegment)

//_____________________________________________________________________________
AliMpVRowSegment::AliMpVRowSegment()
  : AliMpVIndexed()
{
// 
}

//_____________________________________________________________________________
AliMpVRowSegment::~AliMpVRowSegment() {
//  
}

//_____________________________________________________________________________
AliMpVPadIterator* AliMpVRowSegment::CreateIterator() const
{
// Iterator is not yet implemented.
// ---

  Fatal("CreateIterator", "Iterator is not yet implemented.");
  
  return 0;
}  



