// $Id$
// Category: sector
//
// Class AliMpPadRowLSegment
// -------------------------
// Class describing a pad row segment composed of the 
// the identic pads;
// the pads are placed from the offset (defined in the base class)
// to the left.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PAD_ROW_L_SEGMENT_H
#define ALI_MP_PAD_ROW_L_SEGMENT_H

#include <TObject.h>

#include "AliMpVPadRowSegment.h"

class AliMpPadRow;
class AliMpMotif;

class AliMpPadRowLSegment : public AliMpVPadRowSegment
{
  public:
    AliMpPadRowLSegment(AliMpPadRow* padRow, AliMpMotif* motif, Int_t motifPositionId,
                   Int_t nofPads);
    AliMpPadRowLSegment();
    virtual ~AliMpPadRowLSegment();

    // methods
    virtual Double_t  LeftBorderX() const;
    virtual Double_t  RightBorderX() const;

  private:
    // methods
    Double_t  FirstPadCenterX() const;
    Double_t  LastPadCenterX() const;
    Double_t  FirstPadBorderX() const;
    Double_t  LastPadBorderX() const;
    
  ClassDef(AliMpPadRowLSegment,1)  //Row segment
};

#endif //ALI_MP_PAD_ROW_L_SEGMENT_H

