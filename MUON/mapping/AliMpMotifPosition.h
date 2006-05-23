/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifPosition.h,v 1.8 2006/05/23 13:07:42 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifPosition
/// \brief A placed motif.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_POSITION_H
#define ALI_MP_MOTIF_POSITION_H

#include "AliMpVIndexed.h"
#include "AliMpVMotif.h"

#include <TVector2.h>

class AliMpVPadIterator;

class AliMpMotifPosition : public AliMpVIndexed
{
 public:
  AliMpMotifPosition(Int_t id, AliMpVMotif* motif, TVector2 position);
  AliMpMotifPosition();
  virtual ~AliMpMotifPosition();

  // methods
  virtual AliMpVPadIterator* CreateIterator() const;

  // get methods
  Int_t        GetID() const;
  AliMpVMotif* GetMotif() const;
  Bool_t       HasPad(const AliMpIntPair& indices) const;

  // Geometry
  TVector2 Position() const;
  TVector2 Dimensions() const;
  
  // Sets the ID (which is the MANU ID)
  void SetID(Int_t id); 
  // Sets the position.
  void SetPosition(const TVector2& pos);

  void Print(Option_t* option="") const;

 protected:
  AliMpMotifPosition(const AliMpMotifPosition& right);
  AliMpMotifPosition&  operator = (const AliMpMotifPosition& right);

 private:
  // methods
  // data members 
  Int_t         fID;       ///< identifier=manu id
  AliMpVMotif*  fMotif;    ///< motif
  TVector2      fPosition; ///< position

  ClassDef(AliMpMotifPosition,1) // A motif position
};

// inline functions

inline Int_t  AliMpMotifPosition::GetID() const 
{ return fID; }

inline AliMpVMotif*  AliMpMotifPosition::GetMotif() const
{ return fMotif; }
 
inline TVector2 AliMpMotifPosition::Position() const
{ return fPosition; }

inline TVector2 AliMpMotifPosition::Dimensions() const
{ return fMotif->Dimensions(); }

#endif //ALI_MP_MOTIF_POSITION_H
