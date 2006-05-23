/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSectorAreaHPadIterator.h,v 1.6 2006/05/23 13:07:44 ivana Exp $

/// \ingroup sector
/// \class AliMpSectorAreaHPadIterator
/// \brief An iterator over the pads inside a given area in a sector
/// in horizontal direction.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SECTOR_AREA_H_PAD_ITERATOR_H
#define ALI_MP_SECTOR_AREA_H_PAD_ITERATOR_H

#include "AliMpVPadIterator.h"
#include "AliMpArea.h"
#include "AliMpPad.h"

class AliMpSectorSegmentation;

class AliMpSectorAreaHPadIterator : public AliMpVPadIterator
{
  public:
    AliMpSectorAreaHPadIterator(const AliMpSectorSegmentation* segmentation, 
                                const AliMpArea& area);
    AliMpSectorAreaHPadIterator(const AliMpSectorAreaHPadIterator& src);
    AliMpSectorAreaHPadIterator();
    virtual ~AliMpSectorAreaHPadIterator();

    // operators
    AliMpSectorAreaHPadIterator& 
      operator = (const AliMpSectorAreaHPadIterator& right);

    // methods
    virtual void First();
    virtual void Next();
    virtual Bool_t IsDone() const;
    virtual AliMpPad CurrentItem() const;
    virtual void Invalidate();

  private:
    // private methods
    Bool_t IsValid() const;
    void MoveUp();

    // private data members
    const AliMpSectorSegmentation*  fkSegmentation; ///< \brief the sector segmentation 
                                    /// over which we iterate
    //const AliMpArea  fkArea;      ///< \brief the area
                                    /// (const caused problem with CINT)
    AliMpArea  fkArea;              ///< the area
    AliMpPad   fCurrentPad;         ///< the current pad
    Double_t   fCurrentRowPosition; ///< the current row position

 ClassDef(AliMpSectorAreaHPadIterator,1) // iterator over motif's pads
};

#endif // ALI_MP_SECTOR_AREA_H_PAD_ITERATOR_H
