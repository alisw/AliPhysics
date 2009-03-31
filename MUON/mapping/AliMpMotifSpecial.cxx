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

/// \cond CLASSIMP
ClassImp(AliMpMotifSpecial)
/// \endcond


//______________________________________________________________________________
AliMpMotifSpecial::AliMpMotifSpecial(const TString &id, 
                                     AliMpMotifType *motifType)
  : AliMpVMotif(id,motifType),
    fDimensions(),
    fPadDimensionsVector(),
    fPadDimensionsVector2()
  
{
  /// Standard constructor.
}

//______________________________________________________________________________
AliMpMotifSpecial::AliMpMotifSpecial(TRootIOCtor* ioCtor):
  AliMpVMotif(),
  fDimensions(),
  fPadDimensionsVector(ioCtor),
  fPadDimensionsVector2()
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
TVector2 
AliMpMotifSpecial::GetPadDimensionsByIndices(MpPair_t localIndices) const
{
/// Return the dimensions of pad located at the given indices

  return GetPadDimensionsByIndices(AliMp::PairFirst(localIndices), 
                                   AliMp::PairSecond(localIndices));
}

//______________________________________________________________________________
TVector2 
AliMpMotifSpecial::GetPadDimensionsByIndices(Int_t ixLocal, Int_t iyLocal) const
{
/// Return the dimensions of pad located at the given indices

  if ( GetMotifType()->HasPadByLocalIndices(ixLocal, iyLocal) ) {
    if (!fPadDimensionsVector.GetValue(ixLocal, iyLocal)) {
      Warning("GetPadDimensionsByIndices","Indices outside limits");
      return TVector2(0.,0.);
    }
    else      
      return  *((TVector2*)fPadDimensionsVector.GetValue(ixLocal, iyLocal));
  } 
  else {
    Warning("GetPadDimensionsByIndices","Indices outside limits");
    return TVector2(0.,0.);
  }
}

//______________________________________________________________________________
Int_t AliMpMotifSpecial::GetNofPadDimensions() const
{
/// Return number of different pad dimensions in this motif

  return fPadDimensionsVector2.GetEntriesFast();
}  

//______________________________________________________________________________
TVector2 AliMpMotifSpecial::GetPadDimensions(Int_t i) const
{
/// Returns the i-th different pad dimensions 

  if (i<0 || i>GetNofPadDimensions()) {
    AliFatal("Index outside limits.");
    return TVector2();
  }  

  return *((TVector2*) fPadDimensionsVector2[i]);
}  

//______________________________________________________________________________
void AliMpMotifSpecial::CalculateDimensions()
{
  /// Calculate motif dimensions and keep them in fDimensions data

  Int_t i,j;
  Double_t sizeY=0.;
  Double_t sizeX=0.;
  
  Double_t* tabSizeX = new Double_t[GetMotifType()->GetNofPadsY()];
  for (j=0;j<GetMotifType()->GetNofPadsY();++j) tabSizeX[j]=0.0;
  
  for (i=0;i<GetMotifType()->GetNofPadsX();++i) {
    Double_t trSizeY=0.;
    for (j=0;j<GetMotifType()->GetNofPadsY();++j) {
      TVector2 dim = GetPadDimensionsByIndices(i,j);
      trSizeY+=dim.Y();
      tabSizeX[j]+=dim.X();
    }
    if (trSizeY>sizeY) sizeY=trSizeY;
  }
  for (j=0;j<GetMotifType()->GetNofPadsY();++j) {
    if (tabSizeX[j]>sizeX) sizeX = tabSizeX[j];
  }

  delete [] tabSizeX;
  
  fDimensions = TVector2(sizeX,sizeY);
}  

//______________________________________________________________________________
TVector2 AliMpMotifSpecial::Dimensions() const
{
  /// Give the dimension of the motif

  return fDimensions;
}

//______________________________________________________________________________
TVector2 
AliMpMotifSpecial::PadPositionLocal(MpPair_t localIndices) const 
{
  /// Give the local position of the pad number (ix,iy)
  /// (0,0 is the center of the motif)

  return PadPositionLocal(AliMp::PairFirst(localIndices),
                          AliMp::PairSecond(localIndices));

}

//______________________________________________________________________________
TVector2 
AliMpMotifSpecial::PadPositionLocal(Int_t ixLocal, Int_t iyLocal) const 
{
  /// Give the local position of the pad number (ix,iy)
  /// (0,0 is the center of the motif)

  TVector2 dim = GetPadDimensionsByIndices(ixLocal, iyLocal);
  
  Double_t posX= dim.X();
  for (Int_t i=0;i<ixLocal;++i) {
    posX+=2.*GetPadDimensionsByIndices(i,iyLocal).X();
  }
  
  Double_t posY= dim.Y();
  for (Int_t j=0;j<iyLocal;++j) {
    posY+=2.*GetPadDimensionsByIndices(ixLocal,j).Y();
  }

  return TVector2(posX,posY)-Dimensions();

}
//______________________________________________________________________________
MpPair_t AliMpMotifSpecial::PadIndicesLocal(const TVector2& localPos) const
{
  /// Return the pad indices from a given local position
  /// or -1 if this position doesn't correspond to any valid
  /// connection
  ///
  /// *SOLEIL* : This code suppose that
  /// - 1) all cells have the same size along the Y direction
  /// - 2) the column 0 is entierly filled
    

  // First : find the j index
  TVector2 pos = localPos + Dimensions();
  Int_t j=0;
  Double_t y=pos.Y();
  
  while (j<GetMotifType()->GetNofPadsY()) {
    TVector2 padDim = GetPadDimensionsByIndices(0,j);
    y-=2.*padDim.Y();
    if (y<0.) break;
    j++;
  }

  // Test if it's outside limits
  if (j==GetMotifType()->GetNofPadsY()){
    Warning("PadIndicesLocal","The position is outside the motif");
    return -1;
  }
  
  
  // now find the i index, in the j_th row
  Int_t i=0;
  Double_t x=pos.X();
  
  while (i<GetMotifType()->GetNofPadsX()) {
    TVector2 padDim = GetPadDimensionsByIndices(i,j);
    x-=2.*padDim.X();
    if (x<0.) break;
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
                                         const TVector2& dimensions)
{
  /// Set the dimensions of the pad located at \a localIndices to the given
  /// \a dimensions

  SetPadDimensions(AliMp::PairFirst(localIndices), 
                   AliMp::PairSecond(localIndices), 
                   dimensions);  
}

//______________________________________________________________________________
void AliMpMotifSpecial::SetPadDimensions(Int_t ixLocal, Int_t iyLocal,
                                         const TVector2& dimensions)
{
  /// Set the dimensions of the pad located at \a localIndices to the given
  /// \a dimensions
  
  if ( ! GetMotifType()->HasPadByLocalIndices(ixLocal, iyLocal) ) {
    Warning("SetPadDimensions","Pad indices outside limits");
    return;
  }  

  // fill the dimensions map vector
  TVector2* dimensionsObj = new TVector2(dimensions);
  fPadDimensionsVector.Add(ixLocal, iyLocal, dimensionsObj);

  // fill the vector of different pad dimensions
  // only if these dimensions are not yet present
  Bool_t isPresent = false;
  for (Int_t i=0; i<GetNofPadDimensions(); i++) {
    if (AliMpConstants::IsEqual(*((TVector2*) fPadDimensionsVector2[i]), dimensions)) 
      isPresent = true;    
  }    
  
  if (!isPresent) fPadDimensionsVector2.Add(dimensionsObj);
}
