// $Id$
// Category: sector
//
// Class AliMpVRowSegment
// ----------------------
// Class describing an interface for a row segment.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_V_ROW_SEGMENT_H
#define ALI_MP_V_ROW_SEGMENT_H

#include <TVector2.h>

#include "AliMpVIndexed.h"

class AliMpRow;
class AliMpVMotif;
class AliMpIntPair;

class AliMpVRowSegment : public AliMpVIndexed
{
  public:
    AliMpVRowSegment();
    virtual ~AliMpVRowSegment();

    // methods  
    virtual Double_t  LeftBorderX() const = 0;
    virtual Double_t  RightBorderX() const = 0;
    virtual Double_t  HalfSizeY() const = 0;
    virtual AliMpVPadIterator* CreateIterator() const;

    // find methods
    virtual AliMpVMotif*  FindMotif(const TVector2& position) const = 0;    
    virtual Int_t     FindMotifPositionId(const TVector2& position) const = 0;
    virtual Bool_t    HasMotifPosition(Int_t motifPositionId) const = 0;
    virtual TVector2  MotifCenter(Int_t motifPositionId) const = 0;

    // geometry
    virtual TVector2  Position() const = 0;
    virtual TVector2  Dimensions() const = 0;
    
    // set methods
    virtual void      SetOffset(const TVector2& offset) = 0;
    virtual void      SetGlobalIndices() = 0;
    virtual Int_t     SetIndicesToMotifPosition(Int_t i, 
                               const AliMpIntPair& indices) = 0;
    
    // get methods
    virtual AliMpRow*  GetRow() const = 0;
    virtual Int_t      GetNofMotifs() const = 0;
    virtual AliMpVMotif*  GetMotif(Int_t i) const = 0;
    virtual Int_t      GetMotifPositionId(Int_t i) const = 0;
    
  ClassDef(AliMpVRowSegment,1)  //Row segment
};

#endif //ALI_MP_V_ROW_SEGMENT_H

