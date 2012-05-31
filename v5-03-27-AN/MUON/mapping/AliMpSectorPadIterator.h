/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSectorPadIterator.h,v 1.7 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpSectorPadIterator
/// \brief An iterator over the pads of a sector
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SECTOR_PAD_ITERATOR_H
#define ALI_MP_SECTOR_PAD_ITERATOR_H

#include "AliMpVPadIterator.h"
#include "AliMpMotifPositionPadIterator.h"

class AliMpSector;
class AliMpMotifPosition;


class AliMpSectorPadIterator : public AliMpVPadIterator
{
  public:
    AliMpSectorPadIterator();
    AliMpSectorPadIterator(const AliMpSector* sector);
    AliMpSectorPadIterator(const AliMpSectorPadIterator& src);
    virtual ~AliMpSectorPadIterator();

    // operators
    AliMpSectorPadIterator& operator = (const AliMpSectorPadIterator& right);

    // methods
    virtual void First();
    virtual void Next();
    virtual Bool_t IsDone() const;
    virtual AliMpPad CurrentItem() const;
    virtual void Invalidate();

  private:
    // private methods
    AliMpMotifPosition* ResetToCurrentMotifPosition();
    Bool_t IsValid() const;

    // private data members
    const AliMpSector*  fkSector;  ///< the sector over which to iterate
    UInt_t              fCurrentIndex; ///< the current motif position index
    AliMpMotifPosition* fMotifPos; ///< the current motif position
    AliMpMotifPositionPadIterator  fIterator; ///< iterator over the current motif type

 ClassDef(AliMpSectorPadIterator,1) // iterator over motif's pads
};

#endif // ALI_MP_SECTOR_PAD_ITERATOR_H
