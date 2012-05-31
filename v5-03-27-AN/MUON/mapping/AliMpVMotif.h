/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpVMotif.h,v 1.8 2006/05/24 13:58:18 ivana Exp $

/// \ingroup motif
/// \class AliMpVMotif
/// \brief Abstract base class for a motif with its unique ID and the motif type.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_V_MOTIF_H
#define ALI_MP_V_MOTIF_H

#include "AliMpEncodePair.h"

#include <TObject.h>
#include <TString.h>

class AliMpMotifType;
class AliMpConnection;

class AliMpVMotif : public TObject
{
 public:
  AliMpVMotif(const TString &id, AliMpMotifType *motifType);
  AliMpVMotif();
  virtual ~AliMpVMotif();

  // Access methods
  AliMpMotifType  *GetMotifType() const;
  TString          GetID() const;

                   /// Return the number of pad dimensions
  virtual Int_t    GetNofPadDimensions() const=0;
                   /// Return the i-th pad x dimension
  virtual Double_t GetPadDimensionX(Int_t i) const=0;
                   /// Return the i-th pad y dimension
  virtual Double_t GetPadDimensionY(Int_t i) const=0;

                   /// Return the dimensions of the pad specified by localIndices
  virtual void GetPadDimensionsByIndices(MpPair_t localIndices,
                   Double_t& dx, Double_t& dy) const=0;
                   /// Return the dimensions of the pad specified by localIndices
  virtual void GetPadDimensionsByIndices(Int_t ixLocal, Int_t iyLocal,
                   Double_t& dx, Double_t& dy) const=0;

  // Geometry
                   /// Return x dimensions
  virtual Double_t DimensionX() const=0;
                   /// Return y dimensions
  virtual Double_t DimensionY() const=0;

  // Other methods
                   /// Fill local position of the pad specified by local indices
  virtual void PadPositionLocal(MpPair_t localIndices,
                      Double_t& posx, Double_t& posy  ) const=0;
                   /// Fill local position of the pad specified by local indices
  virtual void PadPositionLocal(Int_t ixLocal, Int_t iyLocal,
                      Double_t& posx, Double_t& posy  ) const=0;

  AliMpConnection *FindConnectionByLocalPos(
                      Double_t localPosX, Double_t localPosY) const;

                   /// Return local indices of the pad specified by local position
  virtual MpPair_t PadIndicesLocal(Double_t localPosX, Double_t localPosY) const=0;

  virtual void Print(Option_t *option) const;

 private:
  /// Not implemented
  AliMpVMotif(const AliMpVMotif& right);
  /// Not implemented
  AliMpVMotif&  operator = (const AliMpVMotif& right);

  // data members 
  TString         fID;        ///< identifier
  AliMpMotifType *fMotifType; ///< the motif type

  ClassDef(AliMpVMotif,1) // A motif with its ID
};

// inline functions

/// Return the motif type
inline  AliMpMotifType* AliMpVMotif::GetMotifType() const {return fMotifType;}

/// Return the motif identifier
inline  TString  AliMpVMotif::GetID() const {return fID;}

#endif //ALI_MP_V_MOTIF_H
