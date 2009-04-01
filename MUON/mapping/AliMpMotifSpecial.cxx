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
// $MpId: AliMpMotifSpecial.cxx,v 1.12 2006/05/24 13:58:41 ivana Exp $
// Category: motif

//-----------------------------------------------------------------------------
// Class AliMpMotifSpecial
// -----------------------
// Class that defines a motif with its unique ID
// and the motif type.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpMotifSpecial.h"
#include "AliMpMotifType.h"
#include "AliMpEncodePair.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <TString.h>
#include <Riostream.h>
#include <TVector2.h>

/// \cond CLASSIMP
ClassImp(AliMpMotifSpecial)
/// \endcond

Int_t AliMpMotifSpecial::fgkPadDimensions2Size = 5;

//______________________________________________________________________________
AliMpMotifSpecial::AliMpMotifSpecial(const TString &id, 
                                     AliMpMotifType *motifType)
  : AliMpVMotif(id,motifType),
    fDimensionX(0.),
    fDimensionY(0.),
    fPadDimensionsVector(),
    fNofPadDimensions2(0),
    fPadDimensions2X(fgkPadDimensions2Size),
    fPadDimensions2Y(fgkPadDimensions2Size)
  
{
  /// Standard constructor.
}

//______________________________________________________________________________
AliMpMotifSpecial::AliMpMotifSpecial(TRootIOCtor* ioCtor):
  AliMpVMotif(),
  fDimensionX(0.),
  fDimensionY(0.),
  fPadDimensionsVector(ioCtor),
  fNofPadDimensions2(),
  fPadDimensions2X(),
  fPadDimensions2Y()
{
  /// Root IO constructor
}


//______________________________________________________________________________
AliMpMotifSpecial::~AliMpMotifSpecial()
{
  /// Destructor
}


//
// public methods
//

//______________________________________________________________________________
void 
AliMpMotifSpecial::GetPadDimensionsByIndices(MpPair_t localIndices,
                                    Double_t& dx, Double_t& dy) const
{
/// Return the dimensions of pad located at the given indices

  GetPadDimensionsByIndices(AliMp::PairFirst(localIndices), 
                            AliMp::PairSecond(localIndices),
                            dx, dy);
}

//______________________________________________________________________________
void 
AliMpMotifSpecial::GetPadDimensionsByIndices(Int_t ixLocal, Int_t iyLocal,
                                    Double_t& dx, Double_t& dy) const
{
/// Return the dimensions of pad located at the given indices

  if ( GetMotifType()->HasPadByLocalIndices(ixLocal, iyLocal) ) {
    if (!fPadDimensionsVector.GetValue(ixLocal, iyLocal)) {
      Warning("GetPadDimensionsByIndices","Indices outside limits");
      dx = 0.;
      dy = 0.;
    }
    else {     
      dx = ((TVector2*)fPadDimensionsVector.GetValue(ixLocal, iyLocal))->X();
      dy = ((TVector2*)fPadDimensionsVector.GetValue(ixLocal, iyLocal))->Y();
    }  
  } 
  else {
    Warning("GetPadDimensionsByIndices","Indices outside limits");
    dx = 0.;
    dy = 0.;
  }
}

//______________________________________________________________________________
Int_t AliMpMotifSpecial::GetNofPadDimensions() const
{
/// Return number of different pad dimensions in this motif

  return fNofPadDimensions2;
}  

//______________________________________________________________________________
Double_t AliMpMotifSpecial::GetPadDimensionX(Int_t i) const
{
/// Returns the i-th different pad dimensions 

  if ( i < 0 || i > fNofPadDimensions2 ) {
    AliFatal("Index outside limits.");
    return 0;
  }  

  return fPadDimensions2X[i];
}  

//______________________________________________________________________________
Double_t AliMpMotifSpecial::GetPadDimensionY(Int_t i) const
{
/// Returns the i-th different pad dimensions 

  if ( i < 0 || i > fNofPadDimensions2 ) {
    AliFatal("Index outside limits.");
    return 0;
  }  

  return fPadDimensions2Y[i];
}  

//______________________________________________________________________________
void AliMpMotifSpecial::CalculateDimensions()
{
  /// Calculate motif dimensions and keep them in fDimensionX/Y data

  Int_t i,j;
  fDimensionY = 0.;
  fDimensionX = 0.;
  
  Double_t* tabSizeX = new Double_t[GetMotifType()->GetNofPadsY()];
  
  for ( j=0; j<GetMotifType()->GetNofPadsY(); ++j ) tabSizeX[j]=0.0;
  
  for ( i=0; i<GetMotifType()->GetNofPadsX(); ++i ) {
    Double_t trSizeY=0.;
    for ( j=0; j<GetMotifType()->GetNofPadsY(); ++j ) {
      Double_t dimx, dimy;
      GetPadDimensionsByIndices(i,j, dimx, dimy);
      trSizeY += dimy;
      tabSizeX[j] += dimx;
    }
    if ( trSizeY > fDimensionY ) fDimensionY = trSizeY;
  }

  for ( j=0; j<GetMotifType()->GetNofPadsY(); ++j ) {
    if ( tabSizeX[j] > fDimensionX ) fDimensionX = tabSizeX[j];
  }

  delete [] tabSizeX;
}  

//______________________________________________________________________________
Double_t AliMpMotifSpecial::DimensionX() const
{
  /// Give the dimension of the motif

  return fDimensionX;
}

//______________________________________________________________________________
Double_t AliMpMotifSpecial::DimensionY() const
{
  /// Give the dimension of the motif

  return fDimensionY;
}

//______________________________________________________________________________
void
AliMpMotifSpecial::PadPositionLocal(MpPair_t localIndices,
                                    Double_t& posx, Double_t& posy) const 
{
  /// Give the local position of the pad number (ix,iy)
  /// (0,0 is the center of the motif)

  return PadPositionLocal(AliMp::PairFirst(localIndices),
                          AliMp::PairSecond(localIndices),
                          posx, posy);

}

//______________________________________________________________________________
void
AliMpMotifSpecial::PadPositionLocal(Int_t ixLocal, Int_t iyLocal,
                                    Double_t& posx, Double_t& posy) const 
{
  /// Give the local position of the pad number (ix,iy)
  /// (0,0 is the center of the motif)

  Double_t dx0, dy0;
  GetPadDimensionsByIndices(ixLocal, iyLocal, dx0, dy0);
  
  posx = dx0;
  for ( Int_t i=0 ;i<ixLocal; ++i ) {
    Double_t dxi, dyi;
    GetPadDimensionsByIndices(i, iyLocal, dxi, dyi);
    posx += 2.*dxi;
  }
  
  posy = dy0;
  for ( Int_t j=0; j<iyLocal; ++j ) {
    Double_t dxi, dyi;
    GetPadDimensionsByIndices(ixLocal, j, dxi, dyi);
    posy += 2.*dyi;
  }
  
  posx -= DimensionX();
  posy -= DimensionY();
}

//______________________________________________________________________________
MpPair_t 
AliMpMotifSpecial::PadIndicesLocal(Double_t localPosX, Double_t localPosY) const
{
  /// Return the pad indices from a given local position
  /// or -1 if this position doesn't correspond to any valid
  /// connection
  ///
  /// *SOLEIL* : This code suppose that
  /// - 1) all cells have the same size along the Y direction
  /// - 2) the column 0 is entierly filled
    

  // First : find the j index
  Int_t j=0;
  Double_t y = localPosY + DimensionY();
  
  while (j<GetMotifType()->GetNofPadsY()) {
    Double_t padDimX, padDimY;
    GetPadDimensionsByIndices(0, j, padDimX, padDimY);
    y -= 2.*padDimY;
    if ( y < 0. ) break;
    j++;
  }

  // Test if it's outside limits
  if (j==GetMotifType()->GetNofPadsY()){
    Warning("PadIndicesLocal","The position is outside the motif");
    return -1;
  }
  
  
  // now find the i index, in the j_th row
  Int_t i=0;
  Double_t x = localPosX + DimensionX();
  
  while (i<GetMotifType()->GetNofPadsX()) {
    Double_t padDimX, padDimY;
    GetPadDimensionsByIndices(i, j, padDimX, padDimY);
    x -= 2.*padDimX;
    if ( x < 0. ) break;
    i++;
  }
  
  
  // Test if it's outside limits

  if (i==GetMotifType()->GetNofPadsX()){
    Warning("PadIndicesLocal","The position is outside the motif");
    return -1;
  }
   
  // then return the found (i,j)
  return AliMp::Pair(i,j);  
}

//______________________________________________________________________________
void AliMpMotifSpecial::SetPadDimensions(MpPair_t localIndices,
                                         Double_t dx, Double_t dy)
{
  /// Set the dimensions of the pad located at \a localIndices to the given
  /// \a dimensions

  SetPadDimensions(AliMp::PairFirst(localIndices), 
                   AliMp::PairSecond(localIndices), dx, dy);  
}

//______________________________________________________________________________
void AliMpMotifSpecial::SetPadDimensions(Int_t ixLocal, Int_t iyLocal,
                                         Double_t dx, Double_t dy)
{
  /// Set the dimensions of the pad located at \a localIndices to the given
  /// \a dimensions
  
  if ( ! GetMotifType()->HasPadByLocalIndices(ixLocal, iyLocal) ) {
    Warning("SetPadDimensions","Pad indices outside limits");
    return;
  }  

  // fill the dimensions map vector
  TVector2* dimensionsObj = new TVector2(dx, dy);
  fPadDimensionsVector.Add(ixLocal, iyLocal, dimensionsObj);

  // fill the vector of different pad dimensions
  // only if these dimensions are not yet present
  Bool_t isPresent = false;
  for (Int_t i=0; i<GetNofPadDimensions(); i++) {
    if ( AliMpConstants::IsEqual(
            fPadDimensions2X[i], fPadDimensions2Y[i], dx, dy) ) 
      isPresent = true;    
  }    
  
  if (!isPresent) {
    fPadDimensions2X.AddAt(dx, fNofPadDimensions2);
    fPadDimensions2Y.AddAt(dy, fNofPadDimensions2++);
  }  
}
