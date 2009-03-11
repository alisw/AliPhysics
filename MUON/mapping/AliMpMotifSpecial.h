/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifSpecial.h,v 1.11 2006/05/24 13:58:18 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifSpecial
/// \brief A special motif with varying pad dimensions
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_SPECIAL_H
#define ALI_MP_MOTIF_SPECIAL_H

#include "AliMpVMotif.h"
#include "AliMpExMap.h"

#include <TVector2.h>
#include <TObjArray.h>

class TString;

class AliMpMotifSpecial : public AliMpVMotif
{
 public:
  AliMpMotifSpecial(const TString &id, AliMpMotifType *motifType);
  AliMpMotifSpecial(TRootIOCtor* ioCtor);
  virtual ~AliMpMotifSpecial();

  // Access methods
  virtual TVector2 GetPadDimensions(const AliMpIntPair& localIndices) const;
  virtual Int_t    GetNofPadDimensions() const;
  virtual TVector2 GetPadDimensions(Int_t i) const;

  // Set methods
  void SetPadDimensions(const AliMpIntPair& localIndices,
                        const TVector2& dimensions);
  
  // Geometry
  void CalculateDimensions();
  virtual TVector2 Dimensions() const;

  // Other methods
  virtual TVector2     PadPositionLocal(const AliMpIntPair& localIndices) const;
  virtual AliMpIntPair PadIndicesLocal(const TVector2& localPos) const;

 private:
  /// Not implemented
  AliMpMotifSpecial();
  // methods
  Int_t VectorIndex(const AliMpIntPair& indices) const;

  // data members
  TVector2     fDimensions;           ///< motif dimensions
  AliMpExMap   fPadDimensionsVector;  ///< the vector of pad dimensions
  TObjArray    fPadDimensionsVector2; ///< the vector of different pad dimensions

  ClassDef(AliMpMotifSpecial,2) // A motif with its ID
};

#endif //ALI_MP_MOTIF_SPECIAL_H
