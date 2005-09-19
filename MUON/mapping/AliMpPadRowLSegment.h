/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPadRowLSegment.h,v 1.4 2005/08/26 15:43:36 ivana Exp $

/// \ingroup sector
/// \class AliMpPadRowLSegment
/// \brief A left pad row segment composed of the identic pads
///
/// A pad row segment composed of the identic pads;
/// the pads are placed from the offset (defined in the base class)
/// to the left.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

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

