// $Id$
// Category: basic
//
// Class AliMpVPadIterator
// -----------------------
// Abstract base class, which defines an iterator over pads
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_VPAD_ITERATOR_H
#define ALI_MP_VPAD_ITERATOR_H

#include <TObject.h>

#include "AliMpPad.h"

class AliMpVPadIterator : public TObject
{
  public:
    AliMpVPadIterator();
    AliMpVPadIterator(const AliMpVPadIterator& right);
    virtual ~AliMpVPadIterator();     

    // operators
    AliMpVPadIterator& operator = (const AliMpVPadIterator& right);

    // methods
    virtual void First() = 0;
    virtual void Next() = 0;
    virtual Bool_t IsDone() const = 0;
    virtual AliMpPad CurrentItem() const = 0;
    virtual void Invalidate() = 0;
 
  ClassDef(AliMpVPadIterator,1) // abstract pad iterator
};

#endif // ALI_MP_V_PAD_ITERATOR_H
