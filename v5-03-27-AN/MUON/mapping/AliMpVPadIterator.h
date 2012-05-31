/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpVPadIterator.h,v 1.6 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \class AliMpVPadIterator
/// \brief An interface for an iterator over pads
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_V_PAD_ITERATOR_H
#define ALI_MP_V_PAD_ITERATOR_H

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
                 /// Set iterator to the first pad
    virtual void First() = 0;
                 /// Set iterator to the next pad
    virtual void Next() = 0;
                 /// Is iterator done
    virtual Bool_t IsDone() const = 0;
                 /// Return current pad
    virtual AliMpPad CurrentItem() const = 0;
                 /// Invalidate iterator (
    virtual void Invalidate() = 0;
 
  ClassDef(AliMpVPadIterator,1) // abstract pad iterator
};

#endif // ALI_MP_V_PAD_ITERATOR_H
