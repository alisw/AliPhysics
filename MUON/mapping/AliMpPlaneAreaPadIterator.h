// $Id$
// Category: plane
//
// Class AliMpPlaneAreaPadIterator
// -------------------------------
// Class, which defines an iterator over the pads 
// inside a given area in a plane in horizontal direction.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PLANE_AREA_PAD_ITERATOR_H
#define ALI_MP_PLANE_AREA_PAD_ITERATOR_H

#include "AliMpPlaneTypes.h"
#include "AliMpVPadIterator.h"
#include "AliMpTransformPadIterator.h"
#include "AliMpArea.h"
#include "AliMpPad.h"

class AliMpPlaneSegmentation;

class AliMpPlaneAreaPadIterator : public AliMpVPadIterator
{
  public:
    AliMpPlaneAreaPadIterator(const AliMpPlaneSegmentation* segmentation, 
                              const AliMpArea& area);
    AliMpPlaneAreaPadIterator(const AliMpPlaneAreaPadIterator& src);
    AliMpPlaneAreaPadIterator();
    virtual ~AliMpPlaneAreaPadIterator();

    // operators
    AliMpPlaneAreaPadIterator& 
      operator = (const AliMpPlaneAreaPadIterator& right);

    // methods
    virtual void First();
    virtual void Next();
    virtual Bool_t IsDone() const;
    virtual AliMpPad CurrentItem() const;
    virtual void Invalidate();

  private:
    // private methods
    void DecomposeArea();

    // private data members
    const AliMpPlaneSegmentation*  fkPlaneSegmentation; // the plane segmentation 
                                               // over which we iterate
    //const AliMpArea  fkArea; // the area
                               // (const caused problem with CINT)
    AliMpArea                 fkArea;          // the area
    PadIteratorVector         fPadIterators;   // pad iterators
    PadIteratorVectorIterator fCurrentIterator;// the current iterator 
				               // in the vector of pad iterators

 ClassDef(AliMpPlaneAreaPadIterator,1) // iterator over motif's pads
};
#endif // ALI_MP_PLANE_AREA_PAD_ITERATOR_H
