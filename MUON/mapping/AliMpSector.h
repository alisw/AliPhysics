/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpSector.h,v 1.14 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpSector
/// \brief A sector (quadrant) of the MUON chamber of stations 1 and 2.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_SECTOR_H
#define ALI_MP_SECTOR_H

#include <TNamed.h>

#include "AliMpDirection.h"
#include "AliMpPlaneType.h"
#include "AliMpEncodePair.h"

#include <TString.h>
#include <TVector2.h>
#include <TObjArray.h>

class AliMpZone;
class AliMpRow;
class AliMpVRowSegment;
class AliMpVMotif;
class AliMpVPadIterator;
class AliMpMotifMap;

class TArrayI;

class AliMpSector : public TNamed
{
  public:
    AliMpSector(const TString& id, Int_t nofZones, Int_t nofRows,
                AliMp::Direction direction, const TVector2& offset);
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
   
    //
    // get methods

    Int_t       GetNofZones() const;
    AliMpZone*  GetZone(Int_t i) const;    

    Int_t       GetNofRows() const;
    AliMpRow*   GetRow(Int_t i) const;

    AliMp::Direction  GetDirection() const;  
    AliMp::PlaneType  GetPlaneType() const;  

    TVector2        GetMinPadDimensions() const;
    TVector2        GetMaxPadDimensions() const;
    MpPair_t        GetMaxPadIndices() const;
    Int_t           GetNofPads() const;

    AliMpMotifMap*  GetMotifMap() const;
    void            GetAllMotifPositionsIDs(TArrayI& ecn) const;
    
    virtual void Print(Option_t* opt="") const;
    
  Int_t GetNofMotifPositions() const;
    
  private:
    /// Not implemented
    AliMpSector(const AliMpSector& right);
    /// Not implemented
    AliMpSector&  operator = (const AliMpSector& right);

    // methods
    AliMpVRowSegment* FindRowSegment(const TVector2& position) const;
    void SetRowOffsets();
    void SetMotifPositions();
    void SetGlobalIndices();
    void SetMinMaxPadDimensions();
    void SetMaxPadIndices();
    void SetNofPads();

    // data members        
    TString    fID;       ///< sector ID
    TVector2   fOffset;   ///< sector position
    TObjArray  fZones;    ///< zones
    TObjArray  fRows;     ///< rows
    AliMpMotifMap*   fMotifMap;         ///< motif map
    AliMp::Direction fDirection;        ///< the direction of constant pad size
    TVector2         fMinPadDimensions; ///< minimum pad dimensions
    TVector2         fMaxPadDimensions; ///< miximum pad dimensions
    MpPair_t         fLMaxPadIndices;   ///< maximum pad indices    
    Int_t            fNofPads;          ///<  total number of pads

  ClassDef(AliMpSector,2)  // Sector
};

// inline functions

/// Return the direction of constant pad size
inline AliMp::Direction AliMpSector::GetDirection() const 
{ return fDirection; }    

/// Return minimum pad dimensions
inline TVector2   AliMpSector::GetMinPadDimensions() const
{ return fMinPadDimensions; }

/// Return maxmum pad dimensions
inline TVector2   AliMpSector::GetMaxPadDimensions() const
{ return fMaxPadDimensions; }

/// Return maximum pad indices
inline MpPair_t  AliMpSector::GetMaxPadIndices() const
{ return fLMaxPadIndices; }

/// Return total number of pads
inline Int_t  AliMpSector::GetNofPads() const
{ return fNofPads; }

/// Return the motif map
inline AliMpMotifMap* AliMpSector::GetMotifMap() const 
{ return fMotifMap; }    

#endif //ALI_MP_SECTOR_H

