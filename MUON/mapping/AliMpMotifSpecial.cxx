// $Id$
// Category: motif
//
// Class AliMpMotifSpecial
// -----------------------
// Class that defines a motif with its unique ID
// and the motif type.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TString.h>

#include "AliMpMotifSpecial.h"
#include "AliMpMotifType.h"
#include "AliMpIntPair.h"
#include "AliMpConstants.h"

ClassImp(AliMpMotifSpecial)


// private methods
//______________________________________________________________________________
Int_t AliMpMotifSpecial::VectorIndex(const AliMpIntPair& indices) const
{
// transform indices to linear vector index
  return indices.GetFirst()*GetMotifType()->GetNofPadsY() + indices.GetSecond();
}


//public methods

//______________________________________________________________________________
AliMpMotifSpecial::AliMpMotifSpecial():
  AliMpVMotif(),
  fPadDimensionsVector(),
  fPadDimensionsVector2()
{
  //default dummy constructor
}


//______________________________________________________________________________
AliMpMotifSpecial::AliMpMotifSpecial(const TString &id, 
                                     AliMpMotifType *motifType)
  : AliMpVMotif(id,motifType),
    fPadDimensionsVector(),
    fPadDimensionsVector2()
  
{
  // Normal constructor.

#ifdef WITH_STL
  fPadDimensionsVector.resize(motifType->GetNofPadsX()*motifType->GetNofPadsY());
#endif  

#ifdef WITH_ROOT
  fPadDimensionsVector.Expand(motifType->GetNofPadsX()*motifType->GetNofPadsY());
#endif  
}

//______________________________________________________________________________
AliMpMotifSpecial::~AliMpMotifSpecial()
{
  //destructor

#ifdef WITH_ROOT
  fPadDimensionsVector.Delete();
#endif  
}


//______________________________________________________________________________
TVector2 
AliMpMotifSpecial::GetPadDimensions(const AliMpIntPair& localIndices) const
{
// returns the dimensions of pad located at the given indices
  if (GetMotifType()->HasPad(localIndices))
#ifdef WITH_STL
    return fPadDimensionsVector[VectorIndex(localIndices)];
#endif  
#ifdef WITH_ROOT
    return  *((TVector2*)fPadDimensionsVector[VectorIndex(localIndices)]);
#endif  
  else {
    Warning("GetPadDimensions","Indices outside limits");
    return TVector2(0.,0.);
  }
}

//______________________________________________________________________________
Int_t AliMpMotifSpecial::GetNofPadDimensions() const
{
// returns number of different pad dimensions in this motif

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
// returns the i-th different pad dimensions 

  if (i<0 || i>GetNofPadDimensions()) {
    Fatal("GetPadDimensions(i)", "Index outside limits.");
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
TVector2 AliMpMotifSpecial::Dimensions() const
{
  // gives the dimension of the motif


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

  delete tabSizeX;
  
  return TVector2(sizeX,sizeY);
}

//______________________________________________________________________________
TVector2 
AliMpMotifSpecial::PadPositionLocal(const AliMpIntPair& localIndices) const 
{
  // gives the local position of the pad number (ix,iy)
  // (0,0 is the center of the motif)

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
  // return the pad indices from a given local position
  // or AliMpIntPair::Invalid() if this position doesn't correspond to any valid
  // connection

  // *SOLEIL* : This code suppose that
  // 1) all cells have the same size along the Y direction
  // 2) the column 0 is entierly filled
    

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
  // set the dimensions of the pad located at <localIndices> to the given
  // <dimensions>
  
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
  fPadDimensionsVector[VectorIndex(localIndices)]= dimensionsObj;

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
