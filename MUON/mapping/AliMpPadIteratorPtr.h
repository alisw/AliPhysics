// $Id$
// Category: basic
//
// Class AliMpPadIteratorPtr
// --------------------------
// Pointer to the virtual pad iterator;
// enables to allocate the virtual pad iterator on stack.
// Usage:
// AliMpVIndexed* myIndexed = MyIndexed()
// MVIterator& it = *AliMpPadIteratorPtr(myIndexed->CreateIterator());
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PAD_ITERATOR_PTR_H
#define ALI_MP_PAD_ITERATOR_PTR_H

#include <TObject.h>

#include "AliMpVPadIterator.h"

class AliMpPadIteratorPtr : public TObject
{
  public:
    AliMpPadIteratorPtr(AliMpVPadIterator* it);
    // AliMpPadIteratorPtr(const AliMpPadIteratorPtr& right); --> private
    virtual ~AliMpPadIteratorPtr();
  
    AliMpVPadIterator* operator->() { return  fIterator; }
    AliMpVPadIterator& operator*()  { return *fIterator; }

  private:   
    // disallow copy and assignment to avoid
    // multiple deletion of fIterator
    AliMpPadIteratorPtr(const AliMpPadIteratorPtr& right);
    AliMpPadIteratorPtr& operator=(const AliMpPadIteratorPtr& right);
     
    // data members
    AliMpVPadIterator*  fIterator; //The pad iterator
     
  ClassDef(AliMpPadIteratorPtr,1) // Pointer to abstract pad iterator
};

#endif // ALI_MP_PAD_ITERATOR_PTR_H
