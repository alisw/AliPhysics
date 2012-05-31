/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpZone.h,v 1.10 2006/05/24 13:58:21 ivana Exp $

/// \ingroup sector
/// \class AliMpZone
/// \brief A region of pads of the same dimensions composed of subzones.
///
/// The zone contains pads of the same dimensions,
/// it is composed of the subzones.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_ZONE_H
#define ALI_MP_ZONE_H

#include <TObject.h>
#include <TObjArray.h>

class AliMpSubZone;
class AliMpVMotif;

class AliMpZone : public TObject
{
  public:
    AliMpZone(Int_t id);
    AliMpZone();
    virtual ~AliMpZone();
  
    // methods
    void AddSubZone(AliMpSubZone* subZone);

    // find methods
    AliMpSubZone* FindSubZone(const AliMpVMotif* motif) const;
    
    // set methods
    void SetPadDimensions(Double_t dx, Double_t dy);
    
    // access methods
    UInt_t    GetID() const;
    Int_t     GetNofSubZones() const;
    AliMpSubZone*  GetSubZone(Int_t i) const;

    Double_t  GetPadDimensionX() const;
    Double_t  GetPadDimensionY() const;

  private:
    // data members
    UInt_t        fID;           ///< ID
    TObjArray     fSubZones;     ///< subzones
    Double_t      fPadDimensionX;///< pad x dimension
    Double_t      fPadDimensionY;///< pad y dimension

  ClassDef(AliMpZone,2)  // Zone
};

// inline functions

/// Return ID
inline  UInt_t  AliMpZone::GetID() const 
{ return fID; }

/// Return pad x dimensions
inline  Double_t AliMpZone::GetPadDimensionX() const 
{ return fPadDimensionX; }

/// Return pad y dimensions
inline  Double_t AliMpZone::GetPadDimensionY() const 
{ return fPadDimensionY; }

#endif //ALI_MP_ZONE_H
