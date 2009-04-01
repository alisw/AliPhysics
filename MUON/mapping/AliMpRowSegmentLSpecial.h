/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpRowSegmentLSpecial.h,v 1.9 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpRowSegmentLSpecial
/// \brief A special inner row segment composed of the pad rows.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_ROW_SEGMENT_L_SPECIAL_H
#define ALI_MP_ROW_SEGMENT_L_SPECIAL_H

#include "AliMpVRowSegmentSpecial.h"

class AliMpRow;
class AliMpPadRow;
class AliMpVPadRowSegment;

class AliMpRowSegmentLSpecial : public AliMpVRowSegmentSpecial
{
  public:
    AliMpRowSegmentLSpecial(AliMpRow* row, Double_t offsetX);
    AliMpRowSegmentLSpecial();
    virtual ~AliMpRowSegmentLSpecial();
    
    // methods
    virtual void  UpdatePadsOffset();
    virtual Double_t  LeftBorderX() const;
    virtual Double_t  RightBorderX() const;

    // geometry
    virtual Double_t  GetPositionX() const;
    virtual Double_t  GetPositionY() const;

    // set methods
    virtual void   SetGlobalIndices(AliMpRow* rowBefore);
    virtual Int_t  SetIndicesToMotifPosition(Int_t i, MpPair_t indices);

  protected:
    // methods
    virtual void  MotifCenterSlow(Int_t motifPositionId, 
                                  Double_t& x, Double_t& y) const;
    
  private:
    // methods
    AliMpVPadRowSegment* FindMostRightPadRowSegment(Int_t motifPositionId) const;
    
  ClassDef(AliMpRowSegmentLSpecial,1)  // Row segment
};

#endif //ALI_MP_ROW_SEGMENT_L_SPECIAL_H
