// $Id$
// Category: basic
//
// Class AliMpTransformPadIterator
// -------------------------------
// Composite of iterator and transformer.
// Transforms returned pad. 
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_TRANSFORM_PAD_ITERATOR_H
#define ALI_MP_TRANSFORM_PAD_ITERATOR_H

#include "AliMpVPadIterator.h"
#include "AliMpPad.h"

class AliMpTransformer;

class AliMpTransformPadIterator : public AliMpVPadIterator
{
  public:
    AliMpTransformPadIterator(AliMpVPadIterator* it, 
                              const AliMpTransformer* transformer);
    AliMpTransformPadIterator(const AliMpTransformPadIterator& right);
    AliMpTransformPadIterator();
    virtual ~AliMpTransformPadIterator();     

    // operators
    AliMpTransformPadIterator& operator=(const AliMpTransformPadIterator& right);

    // methods
    virtual void First();
    virtual void Next();
    virtual Bool_t IsDone() const;
    virtual AliMpPad CurrentItem() const;
    virtual void Invalidate();
 
  private:
    AliMpVPadIterator*      fIterator;     // iterator
    const AliMpTransformer* fkTransformer; // transformer
 
  ClassDef(AliMpTransformPadIterator,1) // abstract pad iterator
};

#endif // ALI_MP_TRANSFORM_PAD_ITERATOR_H
