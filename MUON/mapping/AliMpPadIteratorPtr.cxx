// $Id$
// Category: basic
//
// Class AliMpPadIteratorPtr
// --------------------------
// Pointer to the virtual pad iterator;
// enables to allocate the virtual pad iterator on stack.
// Usage:
// MVIndexed* myIndexed = MyIndexed()
// MVIterator& it = *AliMpPadIteratorPtr(myIndexed->CreateIterator());
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpPadIteratorPtr.h"

ClassImp(AliMpPadIteratorPtr)

AliMpPadIteratorPtr::AliMpPadIteratorPtr(AliMpVPadIterator* it)
  : fIterator(it)
{}

AliMpPadIteratorPtr::~AliMpPadIteratorPtr() {
//
  delete fIterator;
}

