// $Id$
// Category: sector
//
// Class AliMpRowSegmentSpecial
// ----------------------------
// Class describing a special row segment composed of the 
// pad rows.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_ROW_SEGMENT_SPECIAL_H
#define ALI_MP_ROW_SEGMENT_SPECIAL_H

#include <TVector2.h>

#include "AliMpSectorTypes.h"
#include "AliMpVRowSegment.h"
#include "AliMpVMotif.h"

class AliMpRow;
class AliMpPadRow;
class AliMpPadRowSegment;

class AliMpRowSegmentSpecial : public AliMpVRowSegment
{
  public:
    AliMpRowSegmentSpecial(AliMpRow* row, Double_t offsetX);
    AliMpRowSegmentSpecial();
    virtual ~AliMpRowSegmentSpecial();
    
    // methods
    void  AddPadRow(AliMpPadRow* padRow);
    void  UpdateMotifVector();
    void  UpdatePadsOffset();
    virtual Double_t  LeftBorderX() const;
    virtual Double_t  RightBorderX() const;
    virtual Double_t  HalfSizeY() const;

    // find methods
    virtual AliMpVMotif*  FindMotif(const TVector2& position) const;    
    virtual Int_t     FindMotifPositionId(const TVector2& position) const;
    virtual Bool_t    HasMotifPosition(Int_t motifPositionId) const;
    virtual TVector2  MotifCenter(Int_t motifPositionId) const;

    // geometry
    virtual TVector2  Position() const;
    virtual TVector2  Dimensions() const;

    // set methods
    virtual void      SetOffset(const TVector2& offset) {}
    virtual void      SetGlobalIndices();
    virtual Int_t     SetIndicesToMotifPosition(Int_t i, AliMpIntPair indices);

    // get methods
    virtual AliMpRow*     GetRow() const;
    virtual Int_t         GetNofMotifs() const;
    virtual AliMpVMotif*  GetMotif(Int_t i) const;
    virtual Int_t         GetMotifPositionId(Int_t i) const;

  private:
    // methods
    AliMpPadRow*         FindPadRow(Double_t y) const;
    AliMpPadRowSegment*  FindPadRowSegment(Int_t motifPositionId) const;
    AliMpPadRowSegment*  FindMostRightPadRowSegment(Int_t motifPositionId) const;
    AliMpIntPair         FindRelativeLowIndicesOf(Int_t motifPositionId) const;
    TVector2  MotifCenterSlow(Int_t motifPositionId) const;
    Int_t     GetNofPadRows() const;
    AliMpPadRow*  GetPadRow(Int_t i) const;
    Int_t     MaxNofPadsInRow() const;
    Bool_t    HasMotif(const AliMpVMotif* motif) const;
    
    // data members
    Double_t      fOffsetX; //the x position of the border that touches a standard
                            //row segment
    AliMpRow*     fRow;     //the row containing this segment 
    PadRowVector  fPadRows; //pad rows vector
    MotifVector   fMotifs;  //motifs vector
    MotifPositionIdVector  fMotifPositionIds; //motifs position Ids vector
    
  ClassDef(AliMpRowSegmentSpecial,1)  //Row segment
};

#endif //ALI_MP_ROW_SEGMENT_SPECIAL_H
