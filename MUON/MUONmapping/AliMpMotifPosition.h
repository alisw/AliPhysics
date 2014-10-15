/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifPosition.h,v 1.9 2006/05/24 13:58:18 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifPosition
/// \brief A placed motif.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_POSITION_H
#define ALI_MP_MOTIF_POSITION_H

#include "AliMpVIndexed.h"
#include "AliMpVMotif.h"

class AliMpVPadIterator;

class AliMpMotifPosition : public AliMpVIndexed
{
 public:
  AliMpMotifPosition(Int_t id, AliMpVMotif* motif, Double_t x, Double_t y);
  AliMpMotifPosition();
  virtual ~AliMpMotifPosition();

  // methods
  virtual AliMpVPadIterator* CreateIterator() const;

  // get methods
  Int_t        GetID() const;
  AliMpVMotif* GetMotif() const;
  Bool_t       HasPadByIndices(MpPair_t indices) const;
  Bool_t       HasPadByManuChannel(Int_t manuChannel) const;

  // Geometry
  Double_t  GetPositionX() const;
  Double_t  GetPositionY() const;
  Double_t  GetDimensionX() const;
  Double_t  GetDimensionY() const;
  
  // Sets the ID (which is the MANU ID)
  void SetID(Int_t id); 
  // Sets the position.
  void SetPosition(Double_t x, Double_t y);

  void Print(Option_t* option="") const;

 private:
  /// Not implemented
  AliMpMotifPosition(const AliMpMotifPosition& right);
  /// Not implemented
  AliMpMotifPosition&  operator = (const AliMpMotifPosition& right);

  // methods
  // data members 
  Int_t         fID;        ///< identifier=manu id
  AliMpVMotif*  fMotif;     ///< motif
  Double_t      fPositionX; ///< x position
  Double_t      fPositionY; ///< y position

  ClassDef(AliMpMotifPosition,2) // A motif position
};

// inline functions

/// Return motif position ID = manu id
inline Int_t  AliMpMotifPosition::GetID() const 
{ return fID; }

/// Return motif 
inline AliMpVMotif*  AliMpMotifPosition::GetMotif() const
{ return fMotif; }
 
/// Return x position
inline Double_t AliMpMotifPosition::GetPositionX() const
{ return fPositionX; }

/// Return y position
inline Double_t AliMpMotifPosition::GetPositionY() const
{ return fPositionY; }

/// Return x dimension
inline Double_t AliMpMotifPosition::GetDimensionX() const
{ return fMotif->DimensionX(); }

/// Return y dimension
inline Double_t AliMpMotifPosition::GetDimensionY() const
{ return fMotif->DimensionY(); }

#endif //ALI_MP_MOTIF_POSITION_H
