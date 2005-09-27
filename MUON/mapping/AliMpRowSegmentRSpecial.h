/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpRowSegmentRSpecial.h,v 1.5 2005/08/26 15:43:36 ivana Exp $

/// \ingroup sector
/// \class AliMpRowSegmentRSpecial
/// \brief A special outer row segment composed of the pad rows.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_ROW_SEGMENT_R_SPECIAL_H
#define ALI_MP_ROW_SEGMENT_R_SPECIAL_H

#include <TVector2.h>

#include "AliMpVRowSegmentSpecial.h"

class AliMpRow;
class AliMpPadRow;
class AliMpVPadRowSegment;
class AliMpIntPair;

class AliMpRowSegmentRSpecial : public AliMpVRowSegmentSpecial
{
  public:
    AliMpRowSegmentRSpecial(AliMpRow* row, Double_t offsetX);
    AliMpRowSegmentRSpecial();
    virtual ~AliMpRowSegmentRSpecial();
    
    // methods
     virtual void  UpdatePadsOffset() {}
    virtual Double_t  LeftBorderX() const;
    virtual Double_t  RightBorderX() const;

    // geometry
    virtual TVector2  Position() const;

    // set methods
    virtual void   SetGlobalIndices(AliMpRow* rowBefore);
    virtual Int_t  SetIndicesToMotifPosition(Int_t i, 
                             const AliMpIntPair& indices);

  protected:
    // methods
    virtual TVector2  MotifCenterSlow(Int_t motifPositionId) const;
    
  private:
    // methods
    AliMpVPadRowSegment* FindMostLeftPadRowSegment(Int_t motifPositionId) const;
    void SetGlobalIndicesLow();
    
  ClassDef(AliMpRowSegmentRSpecial,1)  //Row segment
};

#endif //ALI_MP_ROW_SEGMENT_R_SPECIAL_H
