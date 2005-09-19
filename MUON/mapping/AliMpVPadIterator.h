/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpVPadIterator.h,v 1.4 2005/08/26 15:43:36 ivana Exp $

/// \ingroup basic
/// \class AliMpVPadIterator
/// \brief An interface for an iterator over pads
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

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
