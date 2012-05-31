/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifTypePadIterator.h,v 1.8 2006/05/24 13:58:18 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifTypePadIterator
/// \brief An iterator over the pads of a given motif type
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_TYPE_PAD_ITERATOR_H
#define ALI_MP_MOTIF_TYPE_PAD_ITERATOR_H

#include "AliMpVPadIterator.h"

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
    Bool_t  FindFirstPadInLine(Int_t ix, Int_t iy, 
                               Int_t& newIx, Int_t& newIy) const;
    Bool_t  IsValid() const;

    // private data members
    const AliMpMotifType* fkMotifType;///< the motif type over which iterate
    Int_t fCurrentIx;    ///< the current ix position inside the motif type
    Int_t fCurrentIy;    ///< the current iy position inside the motif type

 ClassDef(AliMpMotifTypePadIterator,2) // iterator over motif's pads
};

#endif // ALI_MP_MOTIF_TYPE_PAD_ITERATOR_H
