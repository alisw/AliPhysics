// $Id$
// Category: motif
//
// Class AliMpMotif
// ----------------
// Class that defines a motif with its unique ID
// and the motif type.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpMotif.h"
#include "AliMpMotifType.h"
#include "AliMpIntPair.h"

ClassImp(AliMpMotif)

//_____________________________________________________________________________
AliMpMotif::AliMpMotif()
  : AliMpVMotif(),
    fPadDimensions(TVector2(0.,0.))
{
  //default dummy constructor
}

//_____________________________________________________________________________
AliMpMotif::AliMpMotif(const TString &id, AliMpMotifType *motifType,
	               const TVector2& padDimension)
  : AliMpVMotif(id,motifType),
    fPadDimensions(padDimension)
{
  // Normal constructor.
  // The dimension in a given direction is calculated by
  // multiplying the total dimension by the number of pads

}
//_____________________________________________________________________________
AliMpMotif::~AliMpMotif()
{
  // destructor
}


//_____________________________________________________________________________
TVector2 AliMpMotif::GetPadDimensions(const AliMpIntPair& localIndices) const
{
  // gives the dimension of the specified pad in the motif

  if (GetMotifType()->HasPad(localIndices))
    return fPadDimensions;
  else {
    Warning("GetPadDimensions","indices outside range");
    return TVector2(0.,0.);
  }
}

//_____________________________________________________________________________
TVector2 AliMpMotif::Dimensions() const
{
  // gives the dimension of the motif

  return TVector2(GetMotifType()->GetNofPadsX()*fPadDimensions.X(),
		GetMotifType()->GetNofPadsY()*fPadDimensions.Y());
}

//_____________________________________________________________________________
TVector2 AliMpMotif::PadPositionLocal(const AliMpIntPair& localIndices) const 
{
  // gives the local position of the pad number (ix,iy)
  // (0,0 is the center of the motif)

  TVector2 dim=Dimensions();
  return TVector2((2.*localIndices.GetFirst()+1.)*fPadDimensions.X()-dim.X(),
		  (2.*localIndices.GetSecond()+1.)*fPadDimensions.Y()-dim.Y());
}

//_____________________________________________________________________________
AliMpIntPair AliMpMotif::PadIndicesLocal(const TVector2& localPos) const
{
  // return the pad indices from a given local position
  // or (-1,-1) if this position doesn't correspond to any valid
  // connection

  TVector2 lowerLeft = localPos+Dimensions();
  Int_t ix = (Int_t)(lowerLeft.X()/(2.*fPadDimensions.X()));
  Int_t iy = (Int_t)(lowerLeft.Y()/(2.*fPadDimensions.Y()));
  
  if (!GetMotifType()->FindConnectionByLocalIndices(AliMpIntPair(ix,iy))) {
    //Warning("PadIndicesLocal","Position outside motif");
    return AliMpIntPair::Invalid();
  }
  return AliMpIntPair(ix,iy);
}
