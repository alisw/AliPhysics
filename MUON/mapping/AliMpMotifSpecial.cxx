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
//
// Class AliMpMotifSpecial
// -----------------------
// Class that defines a motif with its unique ID
// and the motif type.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpMotifSpecial.h"
#include "AliMpMotifType.h"
#include "AliMpIntPair.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <TString.h>

/// \cond CLASSIMP
ClassImp(AliMpMotifSpecial)
/// \endcond


//______________________________________________________________________________
AliMpMotifSpecial::AliMpMotifSpecial():
  AliMpVMotif(),
  fDimensions(),
  fPadDimensionsVector(),
  fPadDimensionsVector2()
{
  /// Default constructor
}


//______________________________________________________________________________
AliMpMotifSpecial::AliMpMotifSpecial(const TString &id, 
                                     AliMpMotifType *motifType)
  : AliMpVMotif(id,motifType),
    fDimensions(),
#ifdef WITH_STL
    fPadDimensionsVector(),
#endif    
#ifdef WITH_ROOT
    fPadDimensionsVector(true),
#endif    
    fPadDimensionsVector2()
  
{
  /// Standard constructor.

#ifdef WITH_STL
  fPadDimensionsVector.resize(motifType->GetNofPadsX()*motifType->GetNofPadsY());
#endif  
}

//______________________________________________________________________________
AliMpMotifSpecial::~AliMpMotifSpecial()
{
  /// Destructor
}


//
// private methods
//

//______________________________________________________________________________
Int_t AliMpMotifSpecial::VectorIndex(const AliMpIntPair& indices) const
{
/// Transform indices to linear vector index

  return indices.GetFirst()*GetMotifType()->GetNofPadsY() + indices.GetSecond();
}


//
// public methods
//

#include <Riostream.h>
//______________________________________________________________________________
TVector2 
AliMpMotifSpecial::GetPadDimensions(const AliMpIntPair& localIndices) const
{
/// Return the dimensions of pad located at the given indices

  if (GetMotifType()->HasPad(localIndices)) {
#ifdef WITH_STL
    return fPadDimensionsVector[VectorIndex(localIndices)];
#endif  
#ifdef WITH_ROOT
    if (!fPadDimensionsVector.GetValue(localIndices)) {
      Warning("GetPadDimensions","Indices outside limits");
      return TVector2(0.,0.);
    }
    else      
      return  *((TVector2*)fPadDimensionsVector.GetValue(localIndices));
#endif 
  } 
  else {
    Warning("GetPadDimensions","Indices outside limits");
    return TVector2(0.,0.);
  }
}

//______________________________________________________________________________
Int_t AliMpMotifSpecial::GetNofPadDimensions() const
{
/// Return number of different pad dimensions in this motif

#ifdef WITH_STL
  return fPadDimensionsVector2.size();
#endif  

#ifdef WITH_ROOT
  return fPadDimensionsVector2.GetEntriesFast();
#endif  
}  

//______________________________________________________________________________
TVector2 AliMpMotifSpecial::GetPadDimensions(Int_t i) const
{
/// Returns the i-th different pad dimensions 

  if (i<0 || i>GetNofPadDimensions()) {
    AliFatal("Index outside limits.");
    return TVector2();
  }  

#ifdef WITH_STL
  return fPadDimensionsVector2[i];
#endif  

#ifdef WITH_ROOT
  return *((TVector2*) fPadDimensionsVector2[i]);
#endif  
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
      TVector2 dim = GetPadDimensions(AliMpIntPair(i,j));
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
AliMpMotifSpecial::PadPositionLocal(const AliMpIntPair& localIndices) const 
{
  /// Give the local position of the pad number (ix,iy)
  /// (0,0 is the center of the motif)

  TVector2 dim = GetPadDimensions(localIndices);
  
  Double_t posX= dim.X();
  for (Int_t i=0;i<localIndices.GetFirst();++i) {
    posX+=2.*GetPadDimensions(AliMpIntPair(i,localIndices.GetSecond())).X();
  }
  
  Double_t posY= dim.Y();
  for (Int_t j=0;j<localIndices.GetSecond();++j) {
    posY+=2.*GetPadDimensions(AliMpIntPair(localIndices.GetFirst(),j)).Y();
  }

  return TVector2(posX,posY)-Dimensions();

}
//______________________________________________________________________________
AliMpIntPair AliMpMotifSpecial::PadIndicesLocal(const TVector2& localPos) const
{
  /// Return the pad indices from a given local position
  /// or AliMpIntPair::Invalid() if this position doesn't correspond to any valid
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
    TVector2 padDim = GetPadDimensions(AliMpIntPair(0,j));
    y-=2.*padDim.Y();
    if (y<0.) break;
    j++;
  }

  // Test if it's outside limits
  if (j==GetMotifType()->GetNofPadsY()){
    Warning("PadIndicesLocal","The position is outside the motif");
    return AliMpIntPair::Invalid();
  }
  
  
  // now find the i index, in the j_th row
  Int_t i=0;
  Double_t x=pos.X();
  
  while (i<GetMotifType()->GetNofPadsX()) {
    TVector2 padDim = GetPadDimensions(AliMpIntPair(i,j));
    x-=2.*padDim.X();
    if (x<0.) break;
    i++;
  }
  
  
  // Test if it's outside limits

  if (i==GetMotifType()->GetNofPadsX()){
    Warning("PadIndicesLocal","The position is outside the motif");
    return AliMpIntPair::Invalid();
  }
   
  // then return the found (i,j)
  return AliMpIntPair(i,j);  
}
//______________________________________________________________________________
void AliMpMotifSpecial::SetPadDimensions(const AliMpIntPair& localIndices,
                                         const TVector2& dimensions)
{
  /// Set the dimensions of the pad located at \a localIndices to the given
  /// \a dimensions
  
  if ( !GetMotifType()->HasPad(localIndices)){
    Warning("SetPadDimensions","Pad indices outside limits");
    return;
  }  

  // fill the dimensions map vector
#ifdef WITH_STL
  fPadDimensionsVector[VectorIndex(localIndices)]=dimensions;
  
  // fill the vector of different pad dimensions
  // only if these dimensions are not yet present
  Bool_t isPresent = false;
  for (Int_t i=0; i<GetNofPadDimensions(); i++) {
    if (AliMpConstants::IsEqual(fPadDimensionsVector2[i], dimensions)) 
      isPresent = true;    
  }    
  
  if (!isPresent) fPadDimensionsVector2.push_back(dimensions);
#endif  

#ifdef WITH_ROOT
  TVector2* dimensionsObj = new TVector2(dimensions);
  fPadDimensionsVector.Add(localIndices, dimensionsObj);

  // fill the vector of different pad dimensions
  // only if these dimensions are not yet present
  Bool_t isPresent = false;
  for (Int_t i=0; i<GetNofPadDimensions(); i++) {
    if (AliMpConstants::IsEqual(*((TVector2*) fPadDimensionsVector2[i]), dimensions)) 
      isPresent = true;    
  }    
  
  if (!isPresent) fPadDimensionsVector2.Add(dimensionsObj);
#endif  
  
}
