/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotif.h,v 1.7 2006/05/23 13:07:41 ivana Exp $

/// \ingroup motif
/// \class AliMpMotif
/// \brief A motif with its unique ID and the motif type.
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_H
#define ALI_MP_MOTIF_H

#include "AliMpVMotif.h"

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
  virtual TVector2 GetPadDimensions(const AliMpIntPair& localIndices) const;

  // Geometry
  virtual TVector2 Dimensions() const;

  // Other methods
  virtual TVector2 PadPositionLocal(const AliMpIntPair& localIndices) const;
  virtual AliMpIntPair PadIndicesLocal(const TVector2& localPos) const;

 private:
  // methods

  // data members 
  TVector2  fPadDimensions; ///< pad dimensions (halflength x, y size) 

  ClassDef(AliMpMotif,1) // A motif with its ID
};

// inline functions

inline Int_t    AliMpMotif::GetNofPadDimensions() const 
{ return 1; }

inline TVector2 AliMpMotif::GetPadDimensions(Int_t /*i*/) const 
{ return fPadDimensions; }  

#endif //ALI_MP_MOTIF_H
