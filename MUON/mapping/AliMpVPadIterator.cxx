// $Id$
// Category: basic
//
// Class AliMpVPadIterator
// -----------------------
// Abstract base class, which defines an iterator over pads
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpVPadIterator.h"

ClassImp(AliMpVPadIterator)

//___________________________________________________________________
AliMpVPadIterator::AliMpVPadIterator():
    TObject()
{
// default constructor
}

//___________________________________________________________________
AliMpVPadIterator::AliMpVPadIterator(const AliMpVPadIterator& right)
  : TObject(right)
{
// copy constructor
}

//___________________________________________________________________
AliMpVPadIterator::~AliMpVPadIterator()
{
// destructor
}

//
// operators
//

//___________________________________________________________________
AliMpVPadIterator& 
AliMpVPadIterator::operator = (const AliMpVPadIterator& right)
{
// assignement operator

  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  TObject::operator=(right);

  return *this;
}  

