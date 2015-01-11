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

#include <TObjArray.h>
#include <TArrayD.h>

class TString;

class AliMpMotifSpecial : public AliMpVMotif
{
 public:
  AliMpMotifSpecial(const TString &id, AliMpMotifType *motifType);
  AliMpMotifSpecial(TRootIOCtor* ioCtor);
  virtual ~AliMpMotifSpecial();

  // Access methods
  virtual void GetPadDimensionsByIndices(MpPair_t localIndices,
                      Double_t& dx, Double_t& dy) const;
  virtual void GetPadDimensionsByIndices(Int_t ixLocal, Int_t iyLocal,
                      Double_t& dx, Double_t& dy) const;

  virtual Int_t    GetNofPadDimensions() const;
  virtual Double_t GetPadDimensionX(Int_t i) const;
  virtual Double_t GetPadDimensionY(Int_t i) const;

  // Set methods
  void SetPadDimensions(MpPair_t localIndices, 
                        Double_t dx, Double_t dy);
  void SetPadDimensions(Int_t ixLocal, Int_t iyLocal, 
                        Double_t dx, Double_t dy);
  
  // Geometry
  void CalculateDimensions();

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
  /// Not implemented
  AliMpMotifSpecial();

  // static data members
  static Int_t fgkPadDimensions2Size; ///< The fPadDimensionsX/Y2 array size

  // data members
  Double_t     fDimensionX; ///< motif x dimensions
  Double_t     fDimensionY; ///< motif y dimensions
  AliMpExMap   fPadDimensionsVector;  ///< the vector of pad dimensions
  Int_t        fNofPadDimensions2; ///< number of different pad dimensions
  TArrayD      fPadDimensions2X; ///< the vector of x of different pad dimensions
  TArrayD      fPadDimensions2Y; ///< the vector of y of different pad dimensions

  ClassDef(AliMpMotifSpecial,3) // A motif with its ID
};

#endif //ALI_MP_MOTIF_SPECIAL_H
