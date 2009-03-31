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
    fPadDimensions(TVector2(0.,0.))
{
  /// Default constructor
}

//_____________________________________________________________________________
AliMpMotif::AliMpMotif(const TString &id, AliMpMotifType *motifType,
	               const TVector2& padDimension)
  : AliMpVMotif(id,motifType),
    fPadDimensions(padDimension)
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
TVector2 AliMpMotif::GetPadDimensionsByIndices(MpPair_t localIndices) const
{
  /// Give the dimension of the specified pad in the motif

  if ( GetMotifType()->HasPadByLocalIndices(localIndices) )
    return fPadDimensions;
  else {
    Warning("GetPadDimensionsByIndices","indices outside range");
    return TVector2(0.,0.);
  }
}

//_____________________________________________________________________________
TVector2 AliMpMotif::GetPadDimensionsByIndices(Int_t ixLocal, Int_t iyLocal) const
{
  /// Give the dimension of the specified pad in the motif

  return GetPadDimensionsByIndices(AliMp::Pair(ixLocal, iyLocal));
}

//_____________________________________________________________________________
TVector2 AliMpMotif::Dimensions() const
{
  /// Give the dimension of the motif

  return TVector2(GetMotifType()->GetNofPadsX()*fPadDimensions.X(),
		GetMotifType()->GetNofPadsY()*fPadDimensions.Y());
}

//_____________________________________________________________________________
TVector2 AliMpMotif::PadPositionLocal(MpPair_t localIndices) const 
{
  /// Give the local position of the pad number (ix,iy)
  /// (0,0 is the center of the motif)

  return PadPositionLocal(AliMp::PairFirst(localIndices), 
                          AliMp::PairSecond(localIndices));
}

//_____________________________________________________________________________
TVector2 AliMpMotif::PadPositionLocal(Int_t ixLocal, Int_t iyLocal) const 
{
  /// Give the local position of the pad number (ix,iy)
  /// (0,0 is the center of the motif)

  TVector2 dim=Dimensions();
  return TVector2((2.*ixLocal+1.)*fPadDimensions.X()-dim.X(),
		  (2.*iyLocal+1.)*fPadDimensions.Y()-dim.Y());
}

//_____________________________________________________________________________
MpPair_t AliMpMotif::PadIndicesLocal(const TVector2& localPos) const
{
  /// Return the pad indices from a given local position
  /// or (-1,-1) if this position doesn't correspond to any valid
  /// connection

  TVector2 lowerLeft(localPos);
  
  lowerLeft += Dimensions();

  if ( lowerLeft.X() < - AliMpConstants::LengthTolerance() || 
       lowerLeft.Y() < - AliMpConstants::LengthTolerance() ) 
    {
      return -1;
    }

  Int_t ix = (Int_t)(lowerLeft.X()/(2.*fPadDimensions.X()));
  Int_t iy = (Int_t)(lowerLeft.Y()/(2.*fPadDimensions.Y()));
  
  if ( ! GetMotifType()->FindConnectionByLocalIndices(ix,iy) )
  {
    return -1;
  }
  
  return AliMp::Pair(ix,iy);
}
