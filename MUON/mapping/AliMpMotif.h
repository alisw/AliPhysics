// $Id$
// Category: motif
//
// Class AliMpMotif
// ----------------
// Class that defines a motif with its unique ID
// and the motif type.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef M_MOTIF_H
#define M_MOTIF_H

#include <TObject.h>
#include <TString.h>
#include <TVector2.h>

#include "AliMpVMotif.h"

class AliMpMotif : public AliMpVMotif
{
 public:
  AliMpMotif(const TString &id,AliMpMotifType *motifType, const TVector2& padDimension);
  AliMpMotif();

  // Access methods
  virtual Int_t    GetNofPadDimensions() const;
  virtual TVector2 GetPadDimensions(Int_t i = 0) const;
  virtual TVector2 GetPadDimensions(const AliMpIntPair& localIndices) const;

  // Geometry
  virtual TVector2 Dimensions() const;

  // Other methods
  virtual TVector2 PadPositionLocal(const AliMpIntPair& localIndices) const;
  virtual AliMpIntPair PadIndicesLocal(const TVector2& localPos) const;

 private:
  // methods

  // data members 
  TVector2    fPadDimensions; //pad dimensions (halflength x, y size) 

  ClassDef(AliMpMotif,1) // A motif with its ID
};

// inline functions

inline Int_t    AliMpMotif::GetNofPadDimensions() const { return 1; }
inline TVector2 AliMpMotif::GetPadDimensions(Int_t i) const { return fPadDimensions; }  

#endif //M_MOTIF_H
