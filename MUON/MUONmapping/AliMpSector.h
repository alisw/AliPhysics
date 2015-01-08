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
                AliMp::Direction direction, 
                Double_t offsetx, Double_t offsety);
    AliMpSector();
    virtual ~AliMpSector();
  
    // methods  
    virtual AliMpVPadIterator* CreateIterator() const;
    
    void  SetRowSegmentOffsets();
    void  Initialize(); 
    void  PrintGeometry() const;

    // find methods   
    Int_t  FindMotifPositionId(Double_t x, Double_t y) const;

    AliMpRow*          FindRow(Int_t motifPositionId) const;
    AliMpVRowSegment*  FindRowSegment(Int_t motifPositionId) const;


    // geometry 
    Double_t  GetPositionX() const;
    Double_t  GetPositionY() const;
    Double_t  GetDimensionX() const;
    Double_t  GetDimensionY() const;
   
    //
    // get methods

    Int_t       GetNofZones() const;
    AliMpZone*  GetZone(Int_t i) const;    

    Int_t       GetNofRows() const;
    AliMpRow*   GetRow(Int_t i) const;

    AliMp::Direction  GetDirection() const;  
    AliMp::PlaneType  GetPlaneType() const;  

    Double_t    GetMinPadDimensionX() const;
    Double_t    GetMinPadDimensionY() const;
    Double_t    GetMaxPadDimensionX() const;
    Double_t    GetMaxPadDimensionY() const;
    MpPair_t    GetMaxPadIndices() const;
    Int_t       GetNofPads() const;

    AliMpMotifMap*  GetMotifMap() const;

    Int_t  GetNofMotifPositions() const;
    void   GetAllMotifPositionsIDs(TArrayI& ecn) const;
    
    virtual void Print(Option_t* opt="") const;
    
    
  private:
    /// Not implemented
    AliMpSector(const AliMpSector& right);
    /// Not implemented
    AliMpSector&  operator = (const AliMpSector& right);

    // methods
    AliMpRow*         FindRow(Double_t y) const;    
    AliMpVRowSegment* FindRowSegment(Double_t x, Double_t y) const;

    void SetRowOffsets();
    void SetMotifPositions();
    void SetGlobalIndices();
    void SetMinMaxPadDimensions();
    void SetMaxPadIndices();
    void SetNofPads();
    void SetDimensions();

    // data members        
    TString    fID;       ///< sector ID
    Double_t   fOffsetX;  ///< sector x position
    Double_t   fOffsetY;  ///< sector y position
    Double_t   fDimensionX;  ///< sector x dimension
    Double_t   fDimensionY;  ///< sector y dimension
    TObjArray  fZones;    ///< zones
    TObjArray  fRows;     ///< rows
    AliMpMotifMap*   fMotifMap;         ///< motif map
    AliMp::Direction fDirection;        ///< the direction of constant pad size
    Double_t         fMinPadDimensionX; ///< minimum pad x dimensions
    Double_t         fMinPadDimensionY; ///< minimum pad y dimensions
    Double_t         fMaxPadDimensionX; ///< miximum pad x dimensions
    Double_t         fMaxPadDimensionY; ///< miximum pad y dimensions
    MpPair_t         fLMaxPadIndices;   ///< maximum pad indices    
    Int_t            fNofPads;          ///<  total number of pads

  ClassDef(AliMpSector,3)  // Sector
};

// inline functions

/// Return the direction of constant pad size
inline AliMp::Direction AliMpSector::GetDirection() const 
{ return fDirection; }    

/// Return minimum x pad dimensions
inline Double_t  AliMpSector::GetMinPadDimensionX() const
{ return fMinPadDimensionX; }

/// Return maximum y pad dimensions
inline Double_t  AliMpSector::GetMinPadDimensionY() const
{ return fMinPadDimensionY; }

/// Return maximum x pad dimensions
inline Double_t  AliMpSector::GetMaxPadDimensionX() const
{ return fMaxPadDimensionX; }

/// Return minimum y pad dimensions
inline Double_t  AliMpSector::GetMaxPadDimensionY() const
{ return fMaxPadDimensionY; }

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

