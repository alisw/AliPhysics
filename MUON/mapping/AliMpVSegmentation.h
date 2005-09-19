/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpVSegmentation.h,v 1.6 2005/08/26 15:43:36 ivana Exp $

/// \ingroup basic
/// \class AliMpVSegmentation
/// \brief The abstract base class for the segmentation.
///
/// Provides methods related to pads:
/// conversion between pad indices, pad location, pad position;
/// finding pad neighbour.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_V_SEGMENTATION_H
#define ALI_MP_V_SEGMENTATION_H

#include <TObject.h>

#include "AliMpPadPair.h"
#include "AliMpPad.h"

class TVector2;

class AliMpVPadIterator;
class AliMpIntPair;
class AliMpArea;

class AliMpVSegmentation : public TObject
{
  public:
    AliMpVSegmentation();
    virtual ~AliMpVSegmentation();
  
    // factory method 
    virtual AliMpVPadIterator* CreateIterator(const AliMpArea& area) const = 0;

    // methods  
    virtual AliMpPad PadByLocation(const AliMpIntPair& location, 
                               Bool_t warning) const = 0;
    virtual AliMpPad PadByIndices (const AliMpIntPair& indices,  
                               Bool_t warning) const = 0;
    virtual AliMpPad PadByPosition(const TVector2& position,
                               Bool_t warning) const = 0;

    virtual AliMpPadPair PadsUp(const AliMpPad& pad) const;
    virtual AliMpPadPair PadsDown(const AliMpPad& pad) const;
    virtual AliMpPadPair PadsLeft(const AliMpPad& pad) const;
    virtual AliMpPadPair PadsRight(const AliMpPad& pad) const;

    virtual Int_t  MaxPadIndexX() = 0;
    virtual Int_t  MaxPadIndexY() = 0;

    virtual Bool_t HasPad(const AliMpIntPair& indices) const = 0;
    
  private:  
    // methods
    AliMpPadPair FindPads(const TVector2& position1, 
                          const TVector2& position2) const;

  ClassDef(AliMpVSegmentation,1)  // Segmentation
};

#endif //ALI_MP_V_SEGMENTATION_H

