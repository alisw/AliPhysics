/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpVPadRowSegment.h,v 1.7 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpVPadRowSegment
/// \brief The abstract base class for a pad row segment composed of the 
/// the identic pads.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_V_PAD_ROW_SEGMENT_H
#define ALI_MP_V_PAD_ROW_SEGMENT_H

#include <TObject.h>

class AliMpPadRow;
class AliMpMotif;

class AliMpVPadRowSegment : public TObject
{
  public:
    AliMpVPadRowSegment(AliMpPadRow* padRow, AliMpMotif* motif, 
                       Int_t motifPositionId, Int_t nofPads);
    AliMpVPadRowSegment();
    virtual ~AliMpVPadRowSegment();

    // methods
    virtual Double_t  LeftBorderX() const = 0;
    virtual Double_t  RightBorderX() const = 0;
    virtual Double_t  HalfSizeY() const;

    // get methods
    virtual AliMpPadRow*  GetPadRow() const;
    virtual AliMpMotif*   GetMotif() const;    
    virtual Int_t     GetMotifPositionId() const;
            Int_t     GetNofPads() const {return fNofPads;}     

    // set methods
    void  SetOffsetX(Double_t offsetX);  

  protected:
    AliMpVPadRowSegment(const AliMpVPadRowSegment& right);
    AliMpVPadRowSegment&  operator = (const AliMpVPadRowSegment& right);

    Double_t  GetOffsetX() const { return fOffsetX; }

  private:
    // data members
    Int_t         fNofPads;  ///< number of pads
    Double_t      fOffsetX;  ///< the x position of the right/left border
    AliMpPadRow*  fPadRow;   ///< the pad row containing this segment 
    AliMpMotif*   fMotif;    ///< the motif 
    Int_t         fMotifPositionId;  ///< the motif position id
    
  ClassDef(AliMpVPadRowSegment,1)  // Row segment
};

#endif //ALI_MP_V_PAD_ROW_SEGMENT_H

