// $Id$
// Category: sector
//
// Class AliMpSector
// -----------------
// Class describing the sector of the MUON chamber of station 1.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SECTOR_H
#define ALI_MP_SECTOR_H

#include <TObject.h>
#include <TString.h>
#include <TVector2.h>

#include "AliMpSectorTypes.h"
#include "AliMpDirection.h"

class AliMpZone;
class AliMpRow;
class AliMpVRowSegment;
class AliMpVMotif;
class AliMpVPadIterator;
class AliMpMotifMap;

class AliMpSector : public TObject
{
  public:
    AliMpSector(const TString& id, Int_t nofZones, Int_t nofRows,
                AliMpDirection direction, const TVector2& offset);
    AliMpSector();
    virtual ~AliMpSector();
  
    // methods  
    virtual AliMpVPadIterator* CreateIterator() const;
    void  SetRowSegmentOffsets();
    void  Initialize(); 
    void  PrintGeometry() const;

    // find methods   
    AliMpRow*     FindRow(const TVector2& position) const;    
    AliMpVMotif*  FindMotif(const TVector2& position) const;
    Int_t         FindMotifPositionId(const TVector2& position) const;

    AliMpRow*          FindRow(Int_t motifPositionId) const;
    AliMpVRowSegment*  FindRowSegment(Int_t motifPositionId) const;
    TVector2           FindPosition(Int_t motifPositionId) const;

    AliMpZone*  FindZone(const TVector2& padDimensions) const;

    // geometry 
    TVector2  Position() const;
    TVector2  Dimensions() const;
    TVector2  Offset() const;
   
    // get methods
    Int_t       GetNofZones() const;
    AliMpZone*  GetZone(Int_t i) const;    
    Int_t       GetNofRows() const;
    AliMpRow*   GetRow(Int_t i) const;
    AliMpDirection  GetDirection() const;  
    TVector2        GetMinPadDimensions() const;
    AliMpMotifMap*  GetMotifMap() const;
    
  protected:
    AliMpSector(const AliMpSector& right);
    AliMpSector&  operator = (const AliMpSector& right);

  private:
    // methods
    AliMpVRowSegment* FindRowSegment(const TVector2& position) const;
    void SetRowOffsets();
    void SetMotifPositions();
    void SetGlobalIndices();
    void SetMinPadDimensions();

    // data members        
    TString    fID;       // sector ID
    TVector2   fOffset;   // sector position
    ZoneVector fZones;    // zones
    RowVector  fRows;     // rows
    AliMpMotifMap* fMotifMap; // motif map
    AliMpDirection fDirection;// the direction of constant pad size
    TVector2       fMinPadDimensions; // minimal pad dimensions

  ClassDef(AliMpSector,1)  //Sector
};

// inline functions

inline TVector2  AliMpSector::Offset() const
{ return fOffset; }

inline AliMpDirection AliMpSector::GetDirection() const 
{ return fDirection; }    

inline TVector2   AliMpSector::GetMinPadDimensions() const
{ return fMinPadDimensions; }

inline AliMpMotifMap* AliMpSector::GetMotifMap() const 
{ return fMotifMap; }    

#endif //ALI_MP_SECTOR_H

