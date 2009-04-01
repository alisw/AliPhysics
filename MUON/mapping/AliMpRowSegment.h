/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpRowSegment.h,v 1.10 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpRowSegment
/// \brief A row segment composed of the the identic motifs.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_ROW_SEGMENT_H
#define ALI_MP_ROW_SEGMENT_H

#include "AliMpVRowSegment.h"

class AliMpRow;
class AliMpVMotif;

class AliMpRowSegment : public AliMpVRowSegment
{
  public:
    AliMpRowSegment(AliMpRow* row, AliMpVMotif* motif, 
                Int_t padOffsetX, Int_t padOffsetY, 
                Int_t nofMotifs, Int_t motifPositionId, Int_t motifPositionDId);
    AliMpRowSegment();
    virtual ~AliMpRowSegment();

    // methods
    virtual Double_t  LeftBorderX() const;
    virtual Double_t  RightBorderX() const;
    virtual Double_t  HalfSizeY() const;

    // find methods
    virtual AliMpVMotif*  FindMotif(Double_t x, Double_t y) const;    
    virtual Int_t     FindMotifPositionId(Double_t x, Double_t y) const;
    virtual Bool_t    HasMotifPosition(Int_t motifPositionId) const;
    virtual void      MotifCenter(Int_t motifPositionId,
                                  Double_t& x, Double_t& y) const;

    // geometry
    virtual Double_t  GetPositionX() const;
    virtual Double_t  GetPositionY() const;
    virtual Double_t  GetDimensionX() const;
    virtual Double_t  GetDimensionY() const;

    // set methods
    virtual void      SetOffset(Double_t x, Double_t y);
    virtual void      SetGlobalIndices(AliMpRow* rowBefore);
    virtual Int_t     SetIndicesToMotifPosition(Int_t i, MpPair_t indices);

    // get methods
    virtual AliMpRow*     GetRow() const;
    virtual Int_t         GetNofMotifs() const;
    virtual AliMpVMotif*  GetMotif(Int_t /*i*/) const;
    virtual Int_t         GetMotifPositionId(Int_t i) const;

  private:
    /// Not implemented
    AliMpRowSegment(const AliMpRowSegment& right);
    /// Not implemented
    AliMpRowSegment&  operator = (const AliMpRowSegment& right);

    // methods
    Double_t  FirstMotifCenterX() const;
    Double_t  LastMotifCenterX() const;
    Double_t  MotifCenterX(Int_t motifPositionId) const;
    Double_t  MotifCenterY(Int_t motifPositionId) const;
    Bool_t    IsInside(Double_t x, Double_t y, Bool_t warn = true) const;

    // data members
    Int_t         fNofMotifs;  ///< number of motifs
    MpPair_t      fLPadOffset; ///< the offset in nof pads 
    Double_t      fOffsetX;    ///< \brief the x position of the centre of the first motif
                               /// wrt to left border
    Double_t      fOffsetY;    ///< \brief the y position of the centre of the first motif
                               /// wrt to row center
    AliMpRow*     fRow;        ///< the row containing this segment 
    AliMpVMotif*  fMotif;      ///< the motif 
    Int_t   fMotifPositionId;  ///< the first motif position id
    Int_t   fMotifPositionDId; ///< +1 if ids are increasing, -1 if decreasing
    
  ClassDef(AliMpRowSegment,2)  // Row segment
};

#endif //ALI_MP_ROW_SEGMENT_H

