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

#include "AliMpContainers.h"

#include "AliMpVRowSegment.h"
#include "AliMpVMotif.h"

#include <TVector2.h>
#ifdef WITH_ROOT
#include <TArrayI.h>
#include <TObjArray.h>
#endif

#ifdef WITH_STL
#include <vector>
#endif

class AliMpRow;
class AliMpPadRow;
class AliMpVPadRowSegment;
class AliMpIntPair;

class AliMpVRowSegmentSpecial : public AliMpVRowSegment
{
  public:
#ifdef WITH_STL
    typedef std::vector<AliMpPadRow*>  PadRowVector;
    typedef std::vector<AliMpVMotif*>  MotifVector;
    typedef std::vector<Int_t>         MotifPositionIdVector;
#endif
#ifdef WITH_ROOT
    typedef  TObjArray  PadRowVector;
    typedef  TObjArray  MotifVector;
    typedef  TArrayI    MotifPositionIdVector;
#endif

  public:
    AliMpVRowSegmentSpecial(AliMpRow* row, Double_t offsetX);
    AliMpVRowSegmentSpecial();
    virtual ~AliMpVRowSegmentSpecial();
    
    // methods
    void  AddPadRow(AliMpPadRow* padRow);
    void  UpdateMotifVector();
    virtual void  UpdatePadsOffset() = 0;
    virtual Double_t  LeftBorderX() const = 0;
    virtual Double_t  RightBorderX() const= 0;
    virtual Double_t  HalfSizeY() const;

    // find methods
    virtual AliMpVMotif*  FindMotif(const TVector2& position) const;    
    virtual Int_t     FindMotifPositionId(const TVector2& position) const;
    virtual Bool_t    HasMotifPosition(Int_t motifPositionId) const;
    virtual TVector2  MotifCenter(Int_t motifPositionId) const;

    // geometry
    virtual TVector2  Position() const = 0;
    virtual TVector2  Dimensions() const;

    // set methods
    virtual void   SetOffset(const TVector2& /*offset*/) {}
    virtual void   SetGlobalIndices(AliMpRow* rowBefore) = 0;
    virtual Int_t  SetIndicesToMotifPosition(Int_t i, 
                             const AliMpIntPair& indices) = 0;

    // get methods
    virtual AliMpRow*     GetRow() const;
    virtual Int_t         GetNofMotifs() const;
    virtual AliMpVMotif*  GetMotif(Int_t i) const;
    virtual Int_t         GetMotifPositionId(Int_t i) const;

  protected:
    // methods
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
    AliMpVRowSegmentSpecial(const AliMpVRowSegmentSpecial& right);
    AliMpVRowSegmentSpecial&  operator = (const AliMpVRowSegmentSpecial& right);

#ifdef WITH_ROOT
    // static data members
    static const Int_t  fgkMaxNofMotifPositionIds; // dimension of fMotifPositionIds
#endif    

    // data members
    AliMpRow*     fRow;     ///< the row containing this segment 
    Double_t      fOffsetX; ///< \brief the x position of the border that touches a standard
                            /// row segment
    PadRowVector  fPadRows; ///< pad rows vector
    MotifVector   fMotifs;  ///< motifs vector
    MotifPositionIdVector  fMotifPositionIds; ///< motifs position Ids vector

#ifdef WITH_ROOT
    Int_t  fNofMotifPositionIds; ///< number of motif positions Ids
#endif    
    
  ClassDef(AliMpVRowSegmentSpecial,1)  //Row segment
};

// inline functions

inline Double_t AliMpVRowSegmentSpecial::GetOffsetX() const
{ return fOffsetX; }    

#endif //ALI_MP_V_ROW_SEGMENT_SPECIAL_H
