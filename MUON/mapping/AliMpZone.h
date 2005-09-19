/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpZone.h,v 1.6 2005/08/26 15:43:36 ivana Exp $

/// \ingroup sector
/// \class AliMpZone
/// \brief A region of pads of the same dimensions composed of subzones.
///
/// The zone contains pads of the same dimensions,
/// it is composed of the subzones.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_ZONE_H
#define ALI_MP_ZONE_H

#include <TObject.h>
#include <TVector2.h>

#include "AliMpSectorTypes.h"

class AliMpSubZone;

class AliMpZone : public TObject
{
  public:
    AliMpZone(Int_t id);
    AliMpZone();
    virtual ~AliMpZone();
  
    // methods
    void AddSubZone(AliMpSubZone* subZone);

    // find methods
    AliMpSubZone* FindSubZone(AliMpVMotif* motif) const;
    
    // set methods
    void SetPadDimensions(const TVector2& padDimensions);
    
    // access methods
    UInt_t    GetID() const;
    Int_t     GetNofSubZones() const;
    AliMpSubZone*  GetSubZone(Int_t i) const;
    TVector2  GetPadDimensions() const;

  private:
    // data members
    UInt_t        fID;           // ID
    SubZoneVector fSubZones;     // subzones
    TVector2      fPadDimensions;// pad dimensions

  ClassDef(AliMpZone,1)  //Zone
};

// inline functions

inline  void AliMpZone::SetPadDimensions(const TVector2& padDimensions)
{ fPadDimensions = padDimensions; }

inline  UInt_t  AliMpZone::GetID() const 
{ return fID; }

inline  TVector2  AliMpZone::GetPadDimensions() const 
{ return fPadDimensions;}

#endif //ALI_MP_ZONE_H
