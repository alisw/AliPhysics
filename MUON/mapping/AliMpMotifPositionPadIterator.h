/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifPositionPadIterator.h,v 1.7 2006/05/24 13:58:18 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifPositionPadIterator
/// \brief An iterator over the pads of a given motif position
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_POSITION_PAD_ITERATOR_H
#define ALI_MP_MOTIF_POSITION_PAD_ITERATOR_H

#include "AliMpVPadIterator.h"
#include "AliMpMotifTypePadIterator.h"

class AliMpMotifPosition;

class AliMpMotifPositionPadIterator : public AliMpVPadIterator
{
  public:
    AliMpMotifPositionPadIterator();
    AliMpMotifPositionPadIterator(const AliMpMotifPosition* motifPos);
    AliMpMotifPositionPadIterator(const AliMpMotifPositionPadIterator& right);
    virtual ~AliMpMotifPositionPadIterator();     

    // operators
    AliMpMotifPositionPadIterator& 
      operator = (const AliMpMotifPositionPadIterator& right);

    virtual void First();
    virtual void Next();
    virtual Bool_t IsDone() const;
    virtual AliMpPad CurrentItem() const;
    virtual void Invalidate();

  private:
    // private methods
    Bool_t IsValid() const;

    // private data members
    const AliMpMotifPosition* fkMotifPos;///< the AliMpMotifPosition over which iterate
    AliMpMotifTypePadIterator fIterator; ///< Iterator over the motif type

 ClassDef(AliMpMotifPositionPadIterator,2) // iterator over motif's pads
};

#endif // ALI_MP_MOTIF_POSITION_PAD_ITERATOR_H
