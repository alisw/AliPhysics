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

class TString;

class AliMpMotif : public AliMpVMotif
{
 public:
  AliMpMotif(const TString& id,
             AliMpMotifType* motifType, 
             Double_t dx, Double_t dy);
  AliMpMotif();
  virtual ~AliMpMotif();

  // Access methods
  virtual Int_t    GetNofPadDimensions() const;
  virtual Double_t GetPadDimensionX(Int_t /*i*/ = 0) const;
  virtual Double_t GetPadDimensionY(Int_t /*i*/ = 0) const;

  virtual void GetPadDimensionsByIndices(MpPair_t localIndices,
                   Double_t& dx, Double_t& dy) const;
  virtual void GetPadDimensionsByIndices(Int_t ixLocal, Int_t iyLocal,
                   Double_t& dx, Double_t& dy) const;

  // Geometry
  virtual Double_t DimensionX() const;
  virtual Double_t DimensionY() const;

  // Other methods
  virtual void PadPositionLocal(MpPair_t localIndices,
                      Double_t& posx, Double_t& posy  ) const;
  virtual void PadPositionLocal(Int_t ixLocal, Int_t iyLocal,
                      Double_t& posx, Double_t& posy  ) const;

  virtual MpPair_t PadIndicesLocal(
                      Double_t localPosX, Double_t localPosY) const;

 private:
  // methods

  // data members 
  Double_t  fPadDimensionX; ///< pad dimensions (halflength x, y size) 
  Double_t  fPadDimensionY; ///< pad dimensions (halflength x, y size) 

  ClassDef(AliMpMotif,2) // A motif with its ID
};

// inline functions

                            /// Return 1 as the number of pad dimensions 
inline Int_t    AliMpMotif::GetNofPadDimensions() const 
{ return 1; }

                            /// Return the pad x dimension
inline Double_t AliMpMotif::GetPadDimensionX(Int_t /*i*/) const 
{ return fPadDimensionX; }  

                            /// Return the pad y dimension
inline Double_t AliMpMotif::GetPadDimensionY(Int_t /*i*/) const 
{ return fPadDimensionY; }  

#endif //ALI_MP_MOTIF_H
