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

#include "AliMpPad.h"
#include "AliMpPlaneType.h"
#include "AliMpStationType.h"

class AliMpVPadIterator;
class AliMpArea;
class AliMpMotifPosition;

class TArrayI;
class TVector2;
class TObjArray;

class AliMpVSegmentation : public TObject
{
  public:
    AliMpVSegmentation();
    virtual ~AliMpVSegmentation();
  
    //
    // methods 
    //

    // factory methods
    /// Create iterator over pads in the given area 
    virtual AliMpVPadIterator* CreateIterator(const AliMpArea& area) const = 0;

    /// Create a pad iterator over the whole area
    virtual AliMpVPadIterator* CreateIterator() const = 0;
    
    /// Fill the array with the pads that are neighbours of pad. Returns
    /// the number of neighbours.
    virtual Int_t GetNeighbours(const AliMpPad& pad, TObjArray& neighbours,
                                Bool_t includeSelf=kFALSE,
                                Bool_t includeVoid=kFALSE) const = 0;

            /// Find pad by location
    virtual AliMpPad PadByLocation(Int_t manuId, Int_t manuChannel, 
                                Bool_t warning = true) const = 0;
            /// Find pad by indices
    virtual AliMpPad PadByIndices (Int_t ix, Int_t iy,  
                               Bool_t warning = true) const = 0;
            /// Find pad by position
    virtual AliMpPad PadByPosition(const TVector2& position,
                               Bool_t warning = true) const = 0;
            /// Return true if the pad with given indices exists.
            /// Compared with the PadByIndices method, this one can generally be implemented
            /// faster, as one does not have to create an AliMpPad object... 
    virtual Bool_t HasPadByIndices(Int_t ix, Int_t iy) const;
  
            /// Return true if the pad with given location exists
    virtual Bool_t HasPadByLocation(Int_t manuId, Int_t manuChannel) const;
  
            /// Return maximum pad index in X direction
    virtual Int_t  MaxPadIndexX() const = 0;
            /// Return maximum pad index in Y direction
    virtual Int_t  MaxPadIndexY() const = 0;
            /// Return the number of pads in the detection element
    virtual Int_t  NofPads() const = 0;

            /// Fill the given array with the electronic card IDs
    virtual void GetAllElectronicCardIDs(TArrayI& ecn) const = 0;

            /// Get the number of electronic card IDs 
    virtual Int_t GetNofElectronicCards() const = 0;
    
            /// Whether or not we have a given manu
    virtual Bool_t HasMotifPosition(Int_t manuId) const = 0;
  
            /// Return the position of a given manu (aka motifPosition)
    virtual AliMpMotifPosition* MotifPosition(Int_t manuId) const = 0;

            /// Return the plane type
    virtual AliMp::PlaneType PlaneType() const = 0;
    
            /// Return the station type
    virtual AliMp::StationType StationType() const = 0;

            /// Return the half-sizes of the detection element
    virtual TVector2 Dimensions() const = 0;
    
            /// Return the position of the origine of the detection element
    virtual TVector2 Position() const = 0;

  
  ClassDef(AliMpVSegmentation,1)  // Segmentation
};

#endif //ALI_MP_V_SEGMENTATION_H

