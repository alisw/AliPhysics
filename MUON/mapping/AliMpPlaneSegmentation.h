// $Id$
// Category: plane
//
// Class AliMpPlaneSegmentation
// ----------------------------
// Class describing the segmentation of the plane.
//
// Transformation of pad characteristics according to sectors:
//
//   I.  ( posId,  Guassi ), ( i, j), ( x, y)         II. |  I.
//  II.  ( posId', Guassi'), (-i, j), (-x, y)       _____ | ____
// III.  (-posId,  Guassi),  (-i,-j), (-x,-y)             |
//  IV.  (-posId', Guassi'), ( i,-j), ( x,-y)        III. |  IV.
//   
// Where (posId', Guassi') is the location of the pad
// in the clipped sector.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PLANE_SEGMENTATION_H
#define ALI_MP_PLANE_SEGMENTATION_H

#include <TVector2.h>

#include "AliMpPlaneTypes.h"
#include "AliMpVSegmentation.h"
#include "AliMpTransformer.h"
#include "AliMpIntPair.h"
#include "AliMpPad.h"

class AliMpPlane;
class AliMpSectorSegmentation;
class AliMpVPadIterator;

class AliMpPlaneSegmentation : public AliMpVSegmentation
{
  public:
    AliMpPlaneSegmentation(const AliMpPlane* plane);
    AliMpPlaneSegmentation();
    virtual ~AliMpPlaneSegmentation();

    // factory method  
    virtual AliMpVPadIterator* CreateIterator(const AliMpArea& area) const;

    // methods  
    virtual AliMpPad PadByLocation(const AliMpIntPair& location, 
                                   Bool_t warning = kTRUE) const;
    virtual AliMpPad PadByIndices (const AliMpIntPair& indices,  
                                   Bool_t warning = kTRUE) const;
    virtual AliMpPad PadByPosition(const TVector2& position,
                                   Bool_t warning = kTRUE) const;

    virtual Int_t    Zone(const AliMpPad& pad, Bool_t warning = kTRUE) const;
    virtual TVector2 PadDimensions(Int_t zone, Bool_t warning = kTRUE) const;

    virtual Bool_t HasPad(const AliMpIntPair& indices) const;
    Bool_t CircleTest(const AliMpIntPair& indices) const;

    // get methods
    Int_t GetNofTransformers() const;
    AliMpTransformer* GetTransformer(Int_t i) const;
    AliMpSectorSegmentation* GetSectorSegmentation(
                                   const AliMpIntPair& scale) const;

  private:
    // methods    
    const AliMpTransformer* GetTransformer(const AliMpIntPair& scale) const;
    AliMpIntPair GetScale(const AliMpIntPair& pair) const;
    AliMpIntPair GetScale(const TVector2& vector) const;
    AliMpIntPair GetLocationScale(const AliMpIntPair& location) const;    
    AliMpSectorSegmentation* GetSectorSegmentation(Int_t motifPositionId) const;

    // data members        
    const AliMpPlane*         fkPlane;                 // plane
    AliMpSectorSegmentation*  fFrontSectorSegmentation;// front sector segmentation
    AliMpSectorSegmentation*  fBackSectorSegmentation; // back sector segmentation
    TransformerVector         fTransformers;    // transformer for each quadrant

  ClassDef(AliMpPlaneSegmentation,1)  // Plane segmentation
};

#endif //ALI_MP_PLANE_SEGMENTATION_H

