/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPadRow.h,v 1.10 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpPadRow
/// \brief A pad row composed of the pad row segments.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PAD_ROW_H
#define ALI_MP_PAD_ROW_H

#include <TObject.h>

#include "AliMpContainers.h"

#include "AliMpXDirection.h"

#ifdef WITH_ROOT
#include <TObjArray.h>
#endif

#ifdef WITH_STL
#include <vector>
#endif

class AliMpVPadRowSegment;
class AliMpMotif;

class AliMpPadRow : public TObject
{
  public:
#ifdef WITH_STL
    typedef std::vector<AliMpVPadRowSegment*>  PadRowSegmentVector;
#endif
#ifdef WITH_ROOT
    typedef TObjArray  PadRowSegmentVector;
#endif

  public:
    AliMpPadRow(AliMp::XDirection direction);
    AliMpPadRow();
    virtual ~AliMpPadRow();
  
    // methods
    AliMpVPadRowSegment*  AddPadRowSegment(AliMpMotif* motif, 
                                          Int_t motifPositionId, 
                                          Int_t nofPads);
    AliMpVPadRowSegment*  FindPadRowSegment(Double_t x) const;
    Double_t  HalfSizeY() const;
    
    // set methods
    void  SetID(Int_t id);
    void  SetOffsetX(Double_t offsetX);
    
    // get methods
    Int_t   GetID() const;
    Int_t   GetNofPadRowSegments() const;
    AliMpVPadRowSegment*  GetPadRowSegment(Int_t i) const;
    Int_t   GetNofPads() const;

  private:
    // methods
    Double_t CurrentBorderX() const;

    // data members
    AliMp::XDirection   fDirection; ///< the pad row x direction
    Int_t               fID;        ///< the pad row ID
    Double_t            fOffsetX;   ///< the x position of the border
    PadRowSegmentVector fSegments;  ///< the pad row segments

  ClassDef(AliMpPadRow,1)  // Pad row
};

#endif //ALI_MP_PAD_ROW_H

