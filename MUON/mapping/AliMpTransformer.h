// $Id$
// Category: basic
//
// Class AliMpTransformer
// ----------------------
// Class contains definition of transformation and
// provides functions for transforming pads.
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

#ifndef ALI_MP_TRANSFORMER_H
#define ALI_MP_TRANSFORMER_H

#include <TVector2.h>
#include <TObject.h>

#include "AliMpIntPair.h"
#include "AliMpArea.h"
#include "AliMpPad.h"

class AliMpTransformer : public TObject
{
  public:
    AliMpTransformer(const TVector2& offset, const AliMpIntPair& scale);
    AliMpTransformer();
    virtual ~AliMpTransformer();
  
    // public methods
    AliMpIntPair Scale(const AliMpIntPair& pair) const;
    TVector2     Scale(const TVector2& vector) const;
    AliMpIntPair ScaleLocation(const AliMpIntPair& orig) const;
    AliMpPad     Scale(const AliMpPad& pad) const;    
    TVector2 Transform(const TVector2& vector) const;
    TVector2 ITransform(const TVector2& vector) const; 
    AliMpPad Transform(const AliMpPad& pad) const; 
    AliMpPad ITransform(const AliMpPad& pad) const;
    AliMpArea    CutArea(const AliMpArea& area) const;
    
    // get methods
    TVector2 GetOffset() const;
    AliMpIntPair GetScale() const;
 
  private:
    // methods
    void CutInterval(Double_t x1, Double_t x2, Double_t x0, Double_t s,
                     Double_t& pos, Double_t& dx) const;
  
    // data members
    TVector2       fOffset;  // translation transformation
    AliMpIntPair   fScale;   // reflection transformation    

  ClassDef(AliMpTransformer,1)  // Transformer
};

// inline functions

inline TVector2 AliMpTransformer::GetOffset() const { return fOffset; }
inline AliMpIntPair AliMpTransformer::GetScale() const  { return fScale; }

#endif //ALI_MP_TRANSFORMER_H

