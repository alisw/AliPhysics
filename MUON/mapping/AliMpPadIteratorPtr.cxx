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
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpPadIteratorPtr.h"
#include "AliMpVPadIterator.h"

ClassImp(AliMpPadIteratorPtr)

//_____________________________________________________________________________
AliMpPadIteratorPtr::AliMpPadIteratorPtr(AliMpVPadIterator* it)
  : fIterator(it)
{}

//_____________________________________________________________________________
AliMpPadIteratorPtr::AliMpPadIteratorPtr(const AliMpPadIteratorPtr& right) 
  : TObject(right) {
// 
  Fatal("AliMpPadIteratorPtr", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpPadIteratorPtr::~AliMpPadIteratorPtr() {
//
  delete fIterator;
}

// operators

//_____________________________________________________________________________
AliMpPadIteratorPtr& 
AliMpPadIteratorPtr::operator=(const AliMpPadIteratorPtr& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

