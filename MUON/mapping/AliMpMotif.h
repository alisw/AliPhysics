/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotif.h,v 1.8 2006/05/24 13:58:18 ivana Exp $

/// \ingroup motif
/// \class AliMpMotif
/// \brief A motif with its unique ID and the motif type.
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_H
#define ALI_MP_MOTIF_H

#include "AliMpVMotif.h"
#include "AliMpEncodePair.h"

#include <TObject.h>
#include <TVector2.h>

class TString;

class AliMpMotif : public AliMpVMotif
{
 public:
  AliMpMotif(const TString &id,AliMpMotifType *motifType, const TVector2& padDimension);
  AliMpMotif();
  virtual ~AliMpMotif();

  // Access methods
  virtual Int_t    GetNofPadDimensions() const;
  virtual TVector2 GetPadDimensions(Int_t /*i*/ = 0) const;
  virtual TVector2 GetPadDimensionsByIndices(MpPair_t localIndices) const;
  virtual TVector2 GetPadDimensionsByIndices(Int_t ixLocal, Int_t iyLocal) const;

  // Geometry
  virtual TVector2 Dimensions() const;

  // Other methods
  virtual TVector2 PadPositionLocal(MpPair_t localIndices) const;
  virtual TVector2 PadPositionLocal(Int_t ixLocal, Int_t iyLocal) const;
  virtual MpPair_t PadIndicesLocal(const TVector2& localPos) const;

 private:
  // methods

  // data members 
  TVector2  fPadDimensions; ///< pad dimensions (halflength x, y size) 

  ClassDef(AliMpMotif,1) // A motif with its ID
};

// inline functions

                            /// Return 1 as the number of pad dimensions 
inline Int_t    AliMpMotif::GetNofPadDimensions() const 
{ return 1; }

                            /// Return the pad dimensions 
inline TVector2 AliMpMotif::GetPadDimensions(Int_t /*i*/) const 
{ return fPadDimensions; }  

#endif //ALI_MP_MOTIF_H
