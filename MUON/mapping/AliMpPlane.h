// $Id$
// Category: plane
//
// Class AliMpPlane
// ----------------
// Class represents the plane composed of 4 sector positions:
// 
//   I.  FS                             II. |  I.
//  II.  BS inverted in x             _____ | ____
// III.  FS inverted in x, y                |
//  IV.  BS inverted in y              III. |  IV.
//   
// FS - front sector
// BS - back sector    
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PLANE_H
#define ALI_MP_PLANE_H

#include <TObject.h>

#include "AliMpPlaneTypes.h"
#include "AliMpStationType.h"
#include "AliMpPlaneType.h"

class TVector2;

class AliMpSector;
class AliMpSectorPosition;
class AliMpIntPair;

class AliMpPlane : public TObject
{
  public:
    AliMpPlane(AliMpSector* frontSector, AliMpSector* backSector,
               const TVector2& q1Position, const TVector2& q2Position,
	       const TVector2& q3Position, const TVector2& q4Position);
    AliMpPlane();
    virtual ~AliMpPlane();
    
    // factory methods
    static AliMpPlane* Create(AliMpStationType station, AliMpPlaneType type,
               const TVector2& q1Position, const TVector2& q2Position,
               const TVector2& q3Position, const TVector2& q4Position);
    static AliMpPlane* Create(AliMpStationType station, AliMpPlaneType type);
    
    // methods
    const AliMpSectorPosition* SectorPosition(const AliMpIntPair& scale) const;

    // get methods
    const AliMpSector*   GetFrontSector() const;
    const AliMpSector*   GetBackSector() const;
    Int_t GetNofSectorPositions() const;
    AliMpSectorPosition* GetSectorPosition(Int_t i) const;

  private:
    // data members    
    const AliMpSector*    fkFrontSector;    // front sector in the 1st quadrant
    const AliMpSector*    fkBackSector;     // back sector in the 1st quadrant
    SectorPositionVector  fSectorPositions; // sector positions

  ClassDef(AliMpPlane,1)  //Plane
};

// inline functions

inline const AliMpSector* AliMpPlane::GetFrontSector() const
{ return fkFrontSector; }

inline const AliMpSector* AliMpPlane::GetBackSector() const
{ return fkBackSector; }

#endif //ALI_MP_PLANE_H

