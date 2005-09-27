/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifSpecial.h,v 1.7 2005/09/26 16:10:46 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifSpecial
/// \brief A special motif with varying pad dimensions
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_SPECIAL_H
#define ALI_MP_MOTIF_SPECIAL_H

#include "AliMpContainers.h"

#ifdef WITH_STL
#include <vector>
#endif

#ifdef WITH_ROOT
#include <TObjArray.h>
#include "AliMpExMap.h"
#endif

#include <TVector2.h>

#include "AliMpVMotif.h"

class TString;

class AliMpMotifSpecial : public AliMpVMotif
{
 public:
#ifdef WITH_STL
  typedef std::vector< TVector2 > DimensionsMap;
  typedef std::vector< TVector2 > DimensionsMap2;
#endif    
#ifdef WITH_ROOT
  typedef AliMpExMap DimensionsMap;
  typedef TObjArray  DimensionsMap2;
#endif    

 public:
  AliMpMotifSpecial(const TString &id, AliMpMotifType *motifType);
  AliMpMotifSpecial();
  virtual ~AliMpMotifSpecial();

  // Access methods
  virtual TVector2 GetPadDimensions(const AliMpIntPair& localIndices) const;
  virtual Int_t    GetNofPadDimensions() const;
  virtual TVector2 GetPadDimensions(Int_t i) const;

  // Set methods
  void SetPadDimensions(const AliMpIntPair& localIndices,
                        const TVector2& dimensions);
  
  // Geometry
  virtual TVector2 Dimensions() const;

  // Other methods
  virtual TVector2     PadPositionLocal(const AliMpIntPair& localIndices) const;
  virtual AliMpIntPair PadIndicesLocal(const TVector2& localPos) const;

 private:
  // methods
  Int_t VectorIndex(const AliMpIntPair& indices) const;

  // data members
  DimensionsMap   fPadDimensionsVector;  // the vector of pad dimensions
  DimensionsMap2  fPadDimensionsVector2; // the vector of different pad dimensions

  ClassDef(AliMpMotifSpecial,1) // A motif with its ID
};

#endif //ALI_MP_MOTIF_SPECIAL_H
