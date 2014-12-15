/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSectorAreaVPadIterator.h,v 1.7 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpSectorAreaVPadIterator
/// \brief An iterator over the pads inside a given area in a sector 
/// in vertical direction.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SECTOR_AREA_V_PAD_ITERATOR_H
#define ALI_MP_SECTOR_AREA_V_PAD_ITERATOR_H

#include "AliMpVPadIterator.h"
#include "AliMpArea.h"
#include "AliMpPad.h"

class AliMpSectorSegmentation;

class AliMpSectorAreaVPadIterator : public AliMpVPadIterator
{
  public:
    AliMpSectorAreaVPadIterator(const AliMpSectorSegmentation* segmentation, 
                                const AliMpArea& area);
    AliMpSectorAreaVPadIterator(const AliMpSectorAreaVPadIterator& src);
    AliMpSectorAreaVPadIterator();
    virtual ~AliMpSectorAreaVPadIterator();

    // operators
    AliMpSectorAreaVPadIterator& 
      operator = (const AliMpSectorAreaVPadIterator& right);

    // methods
    virtual void First();
    virtual void Next();
    virtual Bool_t IsDone() const;
    virtual AliMpPad CurrentItem() const;
    virtual void Invalidate();

  private:
    // private methods
    Bool_t IsValid() const;
    void MoveRight();

    // private data members
    const AliMpSectorSegmentation*  fkSegmentation; ///< \brief the sector segmentation 
                                       /// over which we iterate
    //const AliMpArea  fkArea;         ///< \brief the area
                                       /// (const caused problem with CINT)
    AliMpArea  fkArea;                 ///< the area
    AliMpPad   fCurrentPad;            ///< the current pad
    Double_t   fCurrentColumnPosition; ///< the current column position

 ClassDef(AliMpSectorAreaVPadIterator,1) // iterator over motif's pads
};
#endif // ALI_MP_SECTOR_AREA_V_PAD_ITERATOR_H
