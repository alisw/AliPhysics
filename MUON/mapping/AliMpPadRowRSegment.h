// $Id$
// Category: sector
//
// Class AliMpPadRowRSegment
// -------------------------
// Class describing a pad row segment composed of the 
// the identic pads;
// the pads are placed from the offset (defined in the base class)
// to the right.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PAD_ROW_R_SEGMENT_H
#define ALI_MP_PAD_ROW_R_SEGMENT_H

#include <TObject.h>

#include "AliMpVPadRowSegment.h"

class AliMpPadRow;
class AliMpMotif;

class AliMpPadRowRSegment : public AliMpVPadRowSegment
{
  public:
    AliMpPadRowRSegment(AliMpPadRow* padRow, AliMpMotif* motif, Int_t motifPositionId,
                   Int_t nofPads);
    AliMpPadRowRSegment();
    virtual ~AliMpPadRowRSegment();

    // methods
    virtual Double_t  LeftBorderX() const;
    virtual Double_t  RightBorderX() const;

  private:
    // methods
    Double_t  FirstPadCenterX() const;
    Double_t  LastPadCenterX() const;
    Double_t  FirstPadBorderX() const;
    Double_t  LastPadBorderX() const;
    
  ClassDef(AliMpPadRowRSegment,1)  //Row segment
};

#endif //ALI_MP_PAD_ROW_R_SEGMENT_H

