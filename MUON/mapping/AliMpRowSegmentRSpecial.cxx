// $Id$
// Category: sector
//
// Class AliMpRowSegmentRSpecial
// -----------------------------
// Class describing a special outer row segment composed of the 
// pad rows.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>
#include <Riostream.h>

#include "AliMpRowSegmentRSpecial.h"
#include "AliMpRow.h"
#include "AliMpPadRow.h"
#include "AliMpVPadRowSegment.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"

ClassImp(AliMpRowSegmentRSpecial)

//______________________________________________________________________________
AliMpRowSegmentRSpecial::AliMpRowSegmentRSpecial(AliMpRow* row, Double_t offsetX)
  : AliMpVRowSegmentSpecial(row, offsetX)
{
// 
}

//______________________________________________________________________________
AliMpRowSegmentRSpecial::AliMpRowSegmentRSpecial() 
  : AliMpVRowSegmentSpecial()
{
//
}

//______________________________________________________________________________
AliMpRowSegmentRSpecial::~AliMpRowSegmentRSpecial() 
{
//  
}

//
// private methods  
//

//______________________________________________________________________________
AliMpVPadRowSegment*  
AliMpRowSegmentRSpecial::FindMostLeftPadRowSegment(Int_t motifPositionId) const
{
// Find the most left pad row segment with this motifPositionId.
// ---

  AliMpVPadRowSegment* found = 0;

  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);    

    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) { 
      AliMpVPadRowSegment* padRowSegment = padRow->GetPadRowSegment(j);

      if ( padRowSegment->GetMotifPositionId() == motifPositionId &&
           (!found || padRowSegment->LeftBorderX() < found->LeftBorderX()))
	   
        found = padRowSegment;  
    }
  }

  return found;   	
}

//______________________________________________________________________________
void AliMpRowSegmentRSpecial::SetGlobalIndicesLow()
{
// ...
   
  // Last normal row segment in the row
  // (preceding this special row segment)
  AliMpVRowSegment* rowSegment 
    =  GetRow()->GetRowSegment(GetRow()->GetNofRowSegments()-2);
    
  // Set low indices limit to continue indices of the
  // preceding row segment  
  Int_t ix = rowSegment->GetHighIndicesLimit().GetFirst() + 1;
  Int_t iy = rowSegment->GetLowIndicesLimit().GetSecond();
  
  SetLowIndicesLimit(AliMpIntPair(ix, iy));
}

//
// protected methods  
//

//______________________________________________________________________________
TVector2 AliMpRowSegmentRSpecial::MotifCenterSlow(Int_t motifPositionId) const
{
// Returns the coordinates of the motif specified with
// the given position identifier.
// !! Applicable only for motifs that have their most down pad in
// this row segment.
// ---

  // Find the first (left, down) pad row segment with this motifPositionId.
  AliMpVPadRowSegment* downPadRowSegment 
    = FindPadRowSegment(motifPositionId);
  AliMpVPadRowSegment* leftPadRowSegment 
    = FindMostLeftPadRowSegment(motifPositionId);
  
  // Check if the motifPositionId is present 
  if (!downPadRowSegment || !leftPadRowSegment) {
    Error("MotifCenter", "Outside row segment region");
    return 0;
  }

  // Check if both pad row segments have the same motif 
  if (downPadRowSegment->GetMotif() != leftPadRowSegment->GetMotif()) {
    Fatal("MotifCenter", "Outside row segment region");
    return 0;
  }

  // Get position of found row segment
  Double_t x = leftPadRowSegment->LeftBorderX();       
  Double_t y = GetRow()->LowBorderY()  ;   
  
  for (Int_t i=0; i<downPadRowSegment->GetPadRow()->GetID(); i++)
    y += GetPadRow(i)->HalfSizeY()*2.;
    
  // Add motifs dimensions
  x += downPadRowSegment->GetMotif()->Dimensions().X();
  y += downPadRowSegment->GetMotif()->Dimensions().Y();
  
  return TVector2(x, y);
}

//
// public methods  
//

//______________________________________________________________________________
Double_t  AliMpRowSegmentRSpecial::LeftBorderX() const
{
// Returns the x coordinate of the left row segment border
// in global coordinate system.
// ---

  // The right edge of the last normal segment
  Double_t sameBorder = GetOffsetX();

  // Consistence check  
  Double_t leftBorder = DBL_MAX;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
    Double_t border = padRow->GetPadRowSegment(0)->LeftBorderX();
    if (border < leftBorder) leftBorder =  border;
  }  

  if (TMath::Abs(sameBorder - leftBorder) > 1.e-04)  {
    Error("LeftBorderX", "WrongBorder");
    return sameBorder;
  }  
  
  return leftBorder;

}

//______________________________________________________________________________
Double_t  AliMpRowSegmentRSpecial::RightBorderX() const
{
// Returns the x coordinate of the right row segment border
// in global coordinate system.
// ---

  Double_t rightBorder = -DBL_MAX;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
    Double_t border 
      = padRow->GetPadRowSegment(padRow->GetNofPadRowSegments()-1)
        ->RightBorderX();
      
    if (border > rightBorder) rightBorder =  border;
  }  
  
  return rightBorder;
}

//______________________________________________________________________________
TVector2 AliMpRowSegmentRSpecial::Position() const
{
// Returns the position of the row segment centre.
// The centre is defined as the centre of the rectangular
// row segment envelope.
// ---

  // The right edge of the last normal segment
  Double_t x = GetOffsetX()  + Dimensions().X();
  Double_t y = GetRow()->Position().Y();  
    
  return TVector2(x, y);   
}

//______________________________________________________________________________
Int_t AliMpRowSegmentRSpecial::SetIndicesToMotifPosition(Int_t i, 
                                         const AliMpIntPair& indices)
{
// Sets global indices to i-th motif position and returns next index in x.
// ---

  // Update low indices limit for this row segment
  SetGlobalIndicesLow();
  
  // Check for consistence
  if (GetLowIndicesLimit().GetFirst() != indices.GetFirst()) 
    Fatal("SetIndicesToMotifPosition", "Inconsistent indices");

  // Get motif position
  AliMpMotifPosition* motifPosition
    = GetRow()->GetMotifMap()->FindMotifPosition(GetMotifPositionId(i));

  // Low limit
  AliMpIntPair low = GetLowIndicesLimit();
	  
  if (! motifPosition->GetHighIndicesLimit().IsValid()) {   
     motifPosition->SetLowIndicesLimit(low);
  } 
  else {
    if (motifPosition->GetLowIndicesLimit().GetFirst() > low.GetFirst())
      motifPosition->SetLowIndicesLimit(
                        AliMpIntPair(low.GetFirst(),
                                 motifPosition->GetLowIndicesLimit().GetSecond()));

    if (motifPosition->GetLowIndicesLimit().GetSecond() > low.GetSecond())
       motifPosition->SetLowIndicesLimit(
                         AliMpIntPair(motifPosition->GetLowIndicesLimit().GetFirst(),
                                  low.GetSecond()));
  }

  // High limit	     
  AliMpMotifType* motifType = motifPosition->GetMotif()->GetMotifType();  
  AliMpIntPair high 
    = motifPosition->GetLowIndicesLimit()
      + AliMpIntPair(motifType->GetNofPadsX()-1, motifType->GetNofPadsY()-1);            
  motifPosition->SetHighIndicesLimit(high);
  
  // No increment index needed (this is always the last element)
  return indices.GetFirst();
}


