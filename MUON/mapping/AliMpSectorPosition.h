// $Id$
// Category: plane
//
// Class AliMpSectorPosition
// -------------------------
// Class that represents a placed sector.
// Only translation + reflection transformations can
// be applied.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SECTOR_POSITION_H
#define ALI_MP_SECTOR_POSITION_H

#include <TObject.h>
#include <TVector2.h>

#include "AliMpIntPair.h"

class AliMpSector;

class AliMpSectorPosition : public TObject
{
  public:
    AliMpSectorPosition(const AliMpSector* sector, 
                        const TVector2& offset, const AliMpIntPair& scale);
    AliMpSectorPosition();
    virtual ~AliMpSectorPosition();
  
    // get methods
    const AliMpSector* GetSector() const;
    TVector2       GetOffset() const;
    AliMpIntPair   GetScale() const;

  protected:
    AliMpSectorPosition(const AliMpSectorPosition& right);

    // operators
    AliMpSectorPosition& operator=(const AliMpSectorPosition& right);
    
  private:
    // data members
    const AliMpSector* fkSector; // sector
    TVector2           fOffset;  // translation transformation
    AliMpIntPair       fScale;   // reflection transformation
    
  ClassDef(AliMpSectorPosition,1)  //Sector position
};

// inline functions

inline const AliMpSector* AliMpSectorPosition::GetSector() const
{ return fkSector; }

inline TVector2 AliMpSectorPosition::GetOffset() const
{ return fOffset; }

inline AliMpIntPair AliMpSectorPosition::GetScale() const
{ return fScale; }

#endif //ALI_MP_SECTOR_POSITION_H

