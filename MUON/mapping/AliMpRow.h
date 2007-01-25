/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpRow.h,v 1.11 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpRow
/// \brief A row composed of the row segments.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_ROW_H
#define ALI_MP_ROW_H

#include "AliMpContainers.h"

#include "AliMpVIndexed.h"
#include "AliMpDirection.h"

#include <TVector2.h>
#ifdef WITH_ROOT
#include <TList.h>
#endif

#ifdef WITH_STL
#include <vector>
#endif

class AliMpVRowSegment;
class AliMpVPadIterator;
class AliMpMotifPosition;
class AliMpMotifMap;

class AliMpRow : public AliMpVIndexed
{
  public:
#ifdef WITH_STL
    typedef std::vector<AliMpVRowSegment*>  RowSegmentVector;
#endif
#ifdef WITH_ROOT
    typedef TList  RowSegmentVector;
#endif

  public:
    AliMpRow(Int_t id, AliMpMotifMap* motifMap);
    AliMpRow();
    virtual ~AliMpRow();
  
    // methods
    void  AddRowSegment(AliMpVRowSegment* rowSegment);
    void  AddRowSegmentInFront(AliMpVRowSegment* rowSegment);
    AliMpVRowSegment*  FindRowSegment(Double_t x) const;
    Double_t  LowBorderY() const;
    Double_t  UpperBorderY() const;
    virtual AliMpVPadIterator* CreateIterator() const;
    
    void      SetRowSegmentOffsets(const TVector2& offset);
    Double_t  SetOffsetY(Double_t offsetY);
    void      SetMotifPositions();
    void      SetGlobalIndices(AliMp::Direction constPadSizeDirection, 
                               AliMpRow* rowBefore);

    // geometry
    TVector2  Position() const;
    TVector2  Dimensions() const;    

    // get methods
    UInt_t   GetID() const;
    Int_t    GetNofRowSegments() const;
    AliMpVRowSegment*  GetRowSegment(Int_t i) const;
    AliMpMotifMap*     GetMotifMap() const;

  private:
    AliMpRow(const AliMpRow& right);
    AliMpRow&  operator = (const AliMpRow& right);

    // methods
    AliMpVRowSegment*    FindRowSegment(Int_t ix) const;
    AliMpMotifPosition*  FindMotifPosition(AliMpVRowSegment* segment, Int_t ix) const;
    void SetHighIndicesLimits(Int_t iy);
    void CheckEmpty() const;
  
    // data members
    UInt_t            fID;      ///< row ID
    Double_t          fOffsetY; ///< the y position of the centre of motifs
    RowSegmentVector  fSegments;///< row segments
    AliMpMotifMap*    fMotifMap;///< the motif map associated with its sector

  ClassDef(AliMpRow,1)  // Row
};

// inline functions

inline  UInt_t  AliMpRow::GetID() const { return fID; }
inline  AliMpMotifMap*  AliMpRow::GetMotifMap() const { return fMotifMap; }

#endif //ALI_MP_ROW_H

