// $Id$
// Category: sector
//
// Class AliMpPadRow
// -----------------
// Class describing a pad row composed of the pad row segments.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PAD_ROW_H
#define ALI_MP_PAD_ROW_H

#include <TObject.h>

#include "AliMpSectorTypes.h"

class AliMpPadRowSegment;

class AliMpPadRow : public TObject
{
  public:
    AliMpPadRow();
    virtual ~AliMpPadRow();
  
    // methods
    void  AddPadRowSegment(AliMpPadRowSegment* padRowSegment);
    AliMpPadRowSegment*  FindPadRowSegment(Double_t x) const;
    Double_t  HalfSizeY() const;
    
    // set methods
    void  SetID(Int_t id);
    void  SetOffsetX(Double_t offsetX);
    
    // get methods
    Int_t   GetID() const;
    Int_t   GetNofPadRowSegments() const;
    AliMpPadRowSegment*  GetPadRowSegment(Int_t i) const;
    Int_t   GetNofPads() const;

  private:
    // data members
    Int_t               fID;
    Double_t            fOffsetX; //the x position of the right border
    PadRowSegmentVector fSegments;

  ClassDef(AliMpPadRow,1)  //Pad row
};

#endif //ALI_MP_PAD_ROW_H

