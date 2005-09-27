/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifTypePadIterator.h,v 1.6 2005/09/26 16:10:46 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifTypePadIterator
/// \brief An iterator over the pads of a given motif type
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_TYPE_PAD_ITERATOR_H
#define ALI_MP_MOTIF_TYPE_PAD_ITERATOR_H

#include "AliMpVPadIterator.h"
#include "AliMpIntPair.h"

class AliMpMotifType;

class AliMpMotifTypePadIterator : public AliMpVPadIterator
{
  public:
    AliMpMotifTypePadIterator();
    AliMpMotifTypePadIterator(const AliMpMotifType* motifType);
    AliMpMotifTypePadIterator(const AliMpMotifTypePadIterator& right);
    virtual ~AliMpMotifTypePadIterator();     

    // operators
    AliMpMotifTypePadIterator& 
      operator = (const AliMpMotifTypePadIterator& right);

    virtual void First();
    virtual void Next();
    virtual Bool_t IsDone() const;
    virtual AliMpPad CurrentItem() const;
    virtual void Invalidate();

  private:
    // private methods
    AliMpIntPair FindFirstPadInLine(AliMpIntPair indices) const;
    Bool_t IsValid() const;

    // private data members
    const AliMpMotifType* fMotifType;// the motif type over which iterate
    AliMpIntPair fCurrentPosition;   // the current position inside the motif type

 ClassDef(AliMpMotifTypePadIterator,1) // iterator over motif's pads
};

#endif // ALI_MP_MOTIF_TYPE_PAD_ITERATOR_H
