/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSectorSegmentation.h,v 1.15 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpSectorSegmentation
/// \brief A segmentation of the sector.        
///
/// Provides methods related to pads:                                     \n
/// conversion between pad indices, pad location, pad position;
/// finding pad neighbour.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SECTOR_SEGMENTATION_H
#define ALI_MP_SECTOR_SEGMENTATION_H

#include "AliMpVSegmentation.h"
#include "AliMpPad.h"

#include <TExMap.h>

class AliMpSector;
class AliMpMotifPosition;
class AliMpVPadIterator;
class AliMpArea;

class AliMpSectorSegmentation : public AliMpVSegmentation
{
  public:
    AliMpSectorSegmentation(const AliMpSector* sector, Bool_t own = false);
    AliMpSectorSegmentation();
    virtual ~AliMpSectorSegmentation();
    
    // factory methods  
    virtual AliMpVPadIterator* CreateIterator(const AliMpArea& area) const;
    virtual AliMpVPadIterator* CreateIterator() const;

    virtual Int_t GetNeighbours(const AliMpPad& pad, TObjArray& neighbours,
                               Bool_t includeSelf = kFALSE,
                               Bool_t includeVoid = kFALSE) const;
    
    // methods  
    virtual AliMpPad PadByLocation(Int_t manuId, Int_t manuChannel,
                               Bool_t warning = kTRUE) const;
    virtual AliMpPad PadByIndices (Int_t ix, Int_t iy, 
                               Bool_t warning = kTRUE) const;
    virtual AliMpPad PadByPosition(Double_t x, Double_t y,
                               Bool_t warning = kTRUE) const;
    virtual AliMpPad PadByDirection(Double_t startx, Double_t starty, 
                               Double_t distance) const;
 
    virtual Bool_t HasPadByIndices(Int_t ix, Int_t iy) const;
    virtual Bool_t HasPadByLocation(Int_t manuId, Int_t manuChannel) const;
  
    virtual Int_t  MaxPadIndexX() const;
    virtual Int_t  MaxPadIndexY() const;
    virtual Int_t  NofPads() const;

    virtual void   GetAllElectronicCardIDs(TArrayI& ecn) const;
    virtual Int_t  GetNofElectronicCards() const;
    virtual Bool_t HasMotifPosition(Int_t motifPositionID) const;
    virtual AliMpMotifPosition* MotifPosition(Int_t manuId) const;

    virtual AliMp::PlaneType   PlaneType() const;
    virtual AliMp::StationType StationType() const;

    virtual Double_t  GetDimensionX() const;
    virtual Double_t  GetDimensionY() const;
 
    virtual Double_t  GetPositionX() const;
    virtual Double_t  GetPositionY() const;
  
    virtual void Print(Option_t* opt="") const;
    
    Double_t GetMinPadDimensionX() const;
    Double_t GetMinPadDimensionY() const;

    Bool_t CircleTest(Int_t ix, Int_t iy) const;
   
    const AliMpSector* GetSector() const;
  
  private:
    /// Not implemented
    AliMpSectorSegmentation(const AliMpSectorSegmentation& right);
    /// Not implemented
    AliMpSectorSegmentation&  operator = (const AliMpSectorSegmentation& right);

    // methods
    AliMpMotifPosition*  FindMotifPosition(Int_t ix, Int_t iy) const;
    virtual AliMpPad PadByXDirection(Double_t startx, Double_t starty, 
                                     Double_t maxX) const;
    virtual AliMpPad PadByYDirection(Double_t startx, Double_t starty, 
                                     Double_t maxY) const;
 
    // data members   
    const AliMpSector*  fkSector;     ///< Sector
    Bool_t              fIsOwner;     ///< Sector ownership     
    AliMpPad*           fPadBuffer;   ///< The pad buffer
    Int_t               fMaxIndexInX; ///< maximum pad index in x    
    Int_t               fMaxIndexInY; ///< maximum pad index in y  

  ClassDef(AliMpSectorSegmentation,3)  // Segmentation
};


// inline functions

/// Return the sector
inline const AliMpSector* AliMpSectorSegmentation::GetSector() const
{ return fkSector; }

/// Return station type
inline AliMp::StationType AliMpSectorSegmentation::StationType() const
{ return AliMp::kStation12; }


#endif //ALI_MP_SECTOR_SEGMENTATION_H

