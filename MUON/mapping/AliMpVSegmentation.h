/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpVSegmentation.h,v 1.12 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \class AliMpVSegmentation
/// \brief The abstract base class for the segmentation.
///
/// Provides methods related to pads:
/// conversion between pad indices, pad location, pad position;
/// finding pad neighbour.
///
/// \author David Guez, Ivana Hrivnacova, IPN Orsay;
///         Laurent Aphecetche, SUBATECH

#ifndef ALI_MP_V_SEGMENTATION_H
#define ALI_MP_V_SEGMENTATION_H

#include <TObject.h>

#include "AliMpPadPair.h"
#include "AliMpPad.h"
#include "AliMpPlaneType.h"

class AliMpVPadIterator;
class AliMpIntPair;
class AliMpArea;

class TArrayI;
class TVector2;
class TObjArray;

class AliMpVSegmentation : public TObject
{
  public:
    AliMpVSegmentation();
    virtual ~AliMpVSegmentation();
  
    // factory methods
    /// Create a pad iterator over the given area
    virtual AliMpVPadIterator* CreateIterator(const AliMpArea& area) const = 0;
    
    /// Create a pad iterator over the whole area
    virtual AliMpVPadIterator* CreateIterator() const = 0;
    
    /** Fills the array with the pads that are neighbours of pad. Returns
        the number of neighbours. */
    virtual Int_t GetNeighbours(const AliMpPad& pad, TObjArray& neighbours,
                                Bool_t includeSelf=kFALSE,
                                Bool_t includeVoid=kFALSE) const = 0;

    // methods  
    virtual AliMpPad PadByLocation(const AliMpIntPair& location, 
                               Bool_t warning = true) const = 0;
    virtual AliMpPad PadByIndices (const AliMpIntPair& indices,  
                               Bool_t warning = true) const = 0;
    virtual AliMpPad PadByPosition(const TVector2& position,
                               Bool_t warning = true) const = 0;

    virtual AliMpPadPair PadsUp(const AliMpPad& pad) const;
    virtual AliMpPadPair PadsDown(const AliMpPad& pad) const;
    virtual AliMpPadPair PadsLeft(const AliMpPad& pad) const;
    virtual AliMpPadPair PadsRight(const AliMpPad& pad) const;

    virtual Int_t  MaxPadIndexX() const = 0;
    virtual Int_t  MaxPadIndexY() const = 0;
    virtual Int_t  NofPads() const = 0;

    virtual Bool_t HasPad(const AliMpIntPair& indices) const = 0;
    
    virtual void GetAllElectronicCardIDs(TArrayI& ecn) const = 0;

    virtual AliMp::PlaneType PlaneType() const = 0;
    
    /// Gives the half-sizes (in cm) of the underlying detection element.
    virtual TVector2 Dimensions() const = 0;
    
  private:  
    // methods
    AliMpPadPair FindPads(const TVector2& position1, 
                          const TVector2& position2) const;

  ClassDef(AliMpVSegmentation,1)  // Segmentation
};

#endif //ALI_MP_V_SEGMENTATION_H

