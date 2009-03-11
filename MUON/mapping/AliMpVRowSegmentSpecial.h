/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpVRowSegmentSpecial.h,v 1.10 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpVRowSegmentSpecial
/// \brief Abstract base class for a special row segment composed of the 
/// pad rows.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_V_ROW_SEGMENT_SPECIAL_H
#define ALI_MP_V_ROW_SEGMENT_SPECIAL_H

#include "AliMpVRowSegment.h"
#include "AliMpVMotif.h"

#include <TVector2.h>
#include <TArrayI.h>
#include <TObjArray.h>

class AliMpRow;
class AliMpPadRow;
class AliMpVPadRowSegment;
class AliMpIntPair;

class AliMpVRowSegmentSpecial : public AliMpVRowSegment
{
  public:
    AliMpVRowSegmentSpecial(AliMpRow* row, Double_t offsetX);
    AliMpVRowSegmentSpecial();
    virtual ~AliMpVRowSegmentSpecial();
    
    //
    // methods
    //
    void  AddPadRow(AliMpPadRow* padRow);
    void  UpdateMotifVector();
    /// Update pads offset
    virtual void  UpdatePadsOffset() = 0;
    /// Return the x coordinate of the left border in the global coordinate system.
    virtual Double_t  LeftBorderX() const = 0;
    /// Return the x coordinate of the right border in the global coordinate system.
    virtual Double_t  RightBorderX() const= 0;
    /// Return the half size in y of this row segment.
    virtual Double_t  HalfSizeY() const;

    //
    // find methods
    //
    virtual AliMpVMotif*  FindMotif(const TVector2& position) const;    
    virtual Int_t     FindMotifPositionId(const TVector2& position) const;
    virtual Bool_t    HasMotifPosition(Int_t motifPositionId) const;
    virtual TVector2  MotifCenter(Int_t motifPositionId) const;

    //
    // geometry
    //
    /// Return the position of the row segment centre.
    virtual TVector2  Position() const = 0;
    virtual TVector2  Dimensions() const;

    //
    // set methods
    //
    /// Calculate offset
    virtual void   SetOffset(const TVector2& /*offset*/) {}
    /// Set global indices limits.
    virtual void   SetGlobalIndices(AliMpRow* rowBefore) = 0;
    /// Set global indices to i-th motif position and returns next index in x.
    virtual Int_t  SetIndicesToMotifPosition(Int_t i, 
                             const AliMpIntPair& indices) = 0;

    //
    // get methods
    //
    virtual AliMpRow*     GetRow() const;
    virtual Int_t         GetNofMotifs() const;
    virtual AliMpVMotif*  GetMotif(Int_t i) const;
    virtual Int_t         GetMotifPositionId(Int_t i) const;

  protected:
    // methods
    /// Return the coordinates of the motif specified with the given motif position Id                                           \n
    virtual TVector2  MotifCenterSlow(Int_t motifPositionId) const = 0;
    AliMpPadRow*         FindPadRow(Double_t y) const;
    AliMpVPadRowSegment* FindPadRowSegment(Int_t motifPositionId) const;
    AliMpIntPair         FindRelativeLowIndicesOf(Int_t motifPositionId) const;
    Int_t   MaxNofPadsInRow() const;
    Bool_t  HasMotif(const AliMpVMotif* motif) const;    

    // get methods
    Int_t         GetNofPadRows() const;
    AliMpPadRow*  GetPadRow(Int_t i) const;
    Double_t      GetOffsetX() const;

  private:
    /// Not implemented
    AliMpVRowSegmentSpecial(const AliMpVRowSegmentSpecial& right);
    /// Not implemented
    AliMpVRowSegmentSpecial&  operator = (const AliMpVRowSegmentSpecial& right);

    // static data members
    static const Int_t  fgkMaxNofMotifPositionIds; ///< dimension of fMotifPositionIds

    // data members
    AliMpRow*   fRow;     ///< the row containing this segment 
    Double_t    fOffsetX; ///< \brief the x position of the border that touches a standard
                          /// row segment
    TObjArray   fPadRows; ///< pad rows vector
    TObjArray   fMotifs;  ///< motifs vector
    TArrayI     fMotifPositionIds;    ///< motifs position Ids vector
    Int_t       fNofMotifPositionIds; ///< number of motif positions Ids
    
  ClassDef(AliMpVRowSegmentSpecial,1)  //Row segment
};

// inline functions

/// Return the x position of the border that touches a standard row segment
inline Double_t AliMpVRowSegmentSpecial::GetOffsetX() const
{ return fOffsetX; }    

#endif //ALI_MP_V_ROW_SEGMENT_SPECIAL_H
