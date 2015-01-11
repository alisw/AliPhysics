/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
// $MpId: AliMpMotif.cxx,v 1.8 2006/05/24 13:58:41 ivana Exp $
// Category: motif

//-----------------------------------------------------------------------------
// Class AliMpMotif
// ----------------
// Class that defines a motif with its unique ID
// and the motif type.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpMotif.h"

#include "AliMpConstants.h"
#include "AliMpEncodePair.h"
#include "AliMpMotifType.h"

/// \cond CLASSIMP
ClassImp(AliMpMotif)
/// \endcond

//_____________________________________________________________________________
AliMpMotif::AliMpMotif()
  : AliMpVMotif(),
    fPadDimensionX(0.),
    fPadDimensionY(0.)
{
  /// Default constructor
}

//_____________________________________________________________________________
AliMpMotif::AliMpMotif(const TString &id, AliMpMotifType *motifType,
                       Double_t dx, Double_t dy)
  : AliMpVMotif(id,motifType),
    fPadDimensionX(dx),
    fPadDimensionY(dy)
{
  /// Standard constructor.                                                \n
  /// The dimension in a given direction is calculated by
  /// multiplying the total dimension by the number of pads

}
//_____________________________________________________________________________
AliMpMotif::~AliMpMotif()
{
  /// Destructor
}


//_____________________________________________________________________________
void AliMpMotif::GetPadDimensionsByIndices(MpPair_t localIndices,
                                           Double_t& dx, Double_t& dy) const
{
  /// Give the dimension of the specified pad in the motif

  if ( GetMotifType()->HasPadByLocalIndices(localIndices) ) {
    dx = fPadDimensionX;
    dy = fPadDimensionY;
  }   
  else {
    Warning("GetPadDimensionsByIndices","indices outside range");
    dx = 0.;
    dy = 0.;
  }
}

//_____________________________________________________________________________
void AliMpMotif::GetPadDimensionsByIndices(Int_t ixLocal, Int_t iyLocal,
                                           Double_t& dx, Double_t& dy) const
{
  /// Give the dimension of the specified pad in the motif

  GetPadDimensionsByIndices(AliMp::Pair(ixLocal, iyLocal), dx, dy);
}

//_____________________________________________________________________________
Double_t AliMpMotif::DimensionX() const
{
  /// Give the x dimension of the motif

  return GetMotifType()->GetNofPadsX()*fPadDimensionX;
}

//_____________________________________________________________________________
Double_t AliMpMotif::DimensionY() const
{
  /// Give the y dimension of the motif

  return GetMotifType()->GetNofPadsY()*fPadDimensionY;
}

//_____________________________________________________________________________
void AliMpMotif::PadPositionLocal(MpPair_t localIndices,
                         Double_t& posx, Double_t& posy ) const 
{
  /// Give the local position of the pad number (ix,iy)
  /// (0,0 is the center of the motif)

  PadPositionLocal(AliMp::PairFirst(localIndices), 
                   AliMp::PairSecond(localIndices),
                   posx, posy);
}

//_____________________________________________________________________________
void AliMpMotif::PadPositionLocal(Int_t ixLocal, Int_t iyLocal,
                         Double_t& posx, Double_t& posy) const 
{
  /// Give the local position of the pad number (ix,iy)
  /// (0,0 is the center of the motif)

  posx = (2.*ixLocal+1.)*fPadDimensionX - DimensionX();
  posy = (2.*iyLocal+1.)*fPadDimensionY - DimensionY();
}

//_____________________________________________________________________________
MpPair_t AliMpMotif::PadIndicesLocal(Double_t localPosX, Double_t localPosY) const
{
  /// Return the pad indices from a given local position
  /// or (-1,-1) if this position doesn't correspond to any valid
  /// connection

  Double_t lowerLeftX = localPosX;
  Double_t lowerLeftY = localPosY;
  
  lowerLeftX += DimensionX();
  lowerLeftY += DimensionY();

  if ( lowerLeftX < - AliMpConstants::LengthTolerance() || 
       lowerLeftY < - AliMpConstants::LengthTolerance() ) 
    {
      return -1;
    }

  Int_t ix = (Int_t)(lowerLeftX/(2.*fPadDimensionX));
  Int_t iy = (Int_t)(lowerLeftY/(2.*fPadDimensionY));
  
  if ( ! GetMotifType()->FindConnectionByLocalIndices(ix,iy) )
  {
    return -1;
  }
  
  return AliMp::Pair(ix,iy);
}
