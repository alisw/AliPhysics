// $Id$
// Category: sector
//
// Class AliMpSectorSegmentation
// -----------------------------
// Class describing the segmentation of the sector.        
// Provides methods related to pads:
// conversion between pad indices, pad location, pad position;
// finding pad neighbour.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SECTOR_SEGMENTATION_H
#define ALI_MP_SECTOR_SEGMENTATION_H

#include <TVector2.h>

#include "AliMpSectorTypes.h"
#include "AliMpVSegmentation.h"
#include "AliMpPad.h"

class AliMpSector;
class AliMpMotifPosition;
class AliMpVPadIterator;
class AliMpIntPair;
class AliMpArea;

class AliMpSectorSegmentation : public AliMpVSegmentation
{
  public:
    AliMpSectorSegmentation(const AliMpSector* sector);
    AliMpSectorSegmentation();
    virtual ~AliMpSectorSegmentation();
    
    // factory methods  
    virtual AliMpVPadIterator* CreateIterator(const AliMpArea& area) const;
    AliMpVPadIterator* CreateIterator(const AliMpPad& centerPad,
                                  Bool_t includeCenter=kFALSE) const;

    // methods  
    virtual AliMpPad PadByLocation(const AliMpIntPair& location,
                               Bool_t warning = kTRUE) const;
    virtual AliMpPad PadByIndices (const AliMpIntPair& indices, 
                               Bool_t warning = kTRUE) const;
    virtual AliMpPad PadByPosition(const TVector2& position ,
                               Bool_t warning = kTRUE) const;
    virtual AliMpPad PadByDirection(const TVector2& startPosition, 
                               Double_t distance) const;

    virtual Int_t    Zone(const AliMpPad& pad, Bool_t warning = kTRUE) const;
    virtual TVector2 PadDimensions(Int_t zone, Bool_t warning = kTRUE) const;

    virtual Bool_t HasPad(const AliMpIntPair& indices) const;
    Bool_t HasMotifPosition(Int_t motifPositionID) const;
    TVector2 GetMinPadDimensions() const;
    Bool_t CircleTest(const AliMpIntPair& indices) const;
    void   PrintZones() const;

  protected:
    AliMpSectorSegmentation(const AliMpSectorSegmentation& right);
    AliMpSectorSegmentation&  operator = (const AliMpSectorSegmentation& right);

  private:
#ifdef WITH_ROOT
    static const Double_t   fgkS1;  // the separators used for conversion
    static const Double_t   fgkS2;  // of TVector2 to Long_t
    
    // methods
    Long_t    GetIndex(const TVector2& vector2) const;
    TVector2  GetVector(Long_t index) const;
#endif

    // methods
    void  FillPadDimensionsMap();
    AliMpMotifPosition*  FindMotifPosition(const AliMpIntPair& indices) const;
    virtual AliMpPad PadByXDirection(const TVector2& startPosition, 
                                     Double_t maxX) const;
    virtual AliMpPad PadByYDirection(const TVector2& startPosition, 
                                     Double_t maxY) const;
    virtual AliMpVPadIterator* CreateIterator() const;
 
    // data members        
    const AliMpSector*  fkSector;   // Sector
    AliMpPad*           fPadBuffer; // The pad buffer
#ifdef WITH_STL
    PadDimensionsMap    fPadDimensionsMap; //! Map between zone IDs and pad dimensions
                              // EXCLUDED FOR CINT (does not compile on HP)    
#endif    
#ifdef WITH_ROOT
    mutable PadDimensionsMap  fPadDimensionsMap; //  Map between zone IDs and pad dimensions
#endif    

  ClassDef(AliMpSectorSegmentation,1)  // Segmentation
};

#endif //ALI_MP_SECTOR_SEGMENTATION_H

