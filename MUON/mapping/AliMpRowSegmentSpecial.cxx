// $Id$
// Category: sector
//
// Class AliMpRowSegmentSpecial
// ----------------------------
// Class describing a special row segment composed of the 
// pad rows.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpRowSegmentSpecial.h"
#include "AliMpRow.h"
#include "AliMpPadRow.h"
#include "AliMpPadRowSegment.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpConstants.h"

ClassImp(AliMpRowSegmentSpecial)

//______________________________________________________________________________
AliMpRowSegmentSpecial::AliMpRowSegmentSpecial(AliMpRow* row, Double_t offsetX)
  : AliMpVRowSegment(),
    fOffsetX(offsetX),
    fRow(row),
    fPadRows(),
    fMotifs(),
    fMotifPositionIds()
{
// 
}

//______________________________________________________________________________
AliMpRowSegmentSpecial::AliMpRowSegmentSpecial() 
  : AliMpVRowSegment(),
    fOffsetX(0.),
    fRow(0),
    fPadRows(),
    fMotifs(),
    fMotifPositionIds()
{
//
}

//______________________________________________________________________________
AliMpRowSegmentSpecial::~AliMpRowSegmentSpecial() 
{
//  
  for (Int_t i=0; i<GetNofPadRows(); i++)
    delete fPadRows[i];
}

//
// private methods  
//

//______________________________________________________________________________
Int_t AliMpRowSegmentSpecial::GetNofPadRows() const
{
// Returns number of pad rows.
// ---

  return fPadRows.size();
}  

//______________________________________________________________________________
AliMpPadRow* AliMpRowSegmentSpecial::GetPadRow(Int_t i) const
{
// Returns number of pad rows.
// ---

  return fPadRows[i];
}  

//______________________________________________________________________________
Int_t  AliMpRowSegmentSpecial::MaxNofPadsInRow() const 
{ 
// Returns the maximum number of pads in this row segment along the X direction
// ---

  Int_t maxNofPads = 0;    

  for (Int_t i=0; i<GetNofPadRows(); i++){
    Int_t nofPads = GetPadRow(i)->GetNofPads(); 

    // Find maximum
    if (nofPads > maxNofPads) maxNofPads = nofPads;
  }
    
  return maxNofPads;
}

//______________________________________________________________________________
AliMpPadRow*  AliMpRowSegmentSpecial::FindPadRow(Double_t y) const
{
// Finds the pad row in the given y coordinate.
// ---

  Double_t lowBorder =  fRow->LowBorderY();
  Double_t highBorder = fRow->LowBorderY();
  
  for (Int_t i=0; i<GetNofPadRows(); i++) {    

    AliMpPadRow* padRow = GetPadRow(i);
    highBorder += 2.*padRow->HalfSizeY();

    if ( y >= lowBorder &&  y <= highBorder)
      return padRow;

    lowBorder = highBorder;
  }
  
  return 0;   	
}

//______________________________________________________________________________
AliMpPadRowSegment*  
AliMpRowSegmentSpecial::FindPadRowSegment(Int_t motifPositionId) const
{
// Find the most down pad row segment with this motifPositionId.
// ---

  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);    

    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) { 
      AliMpPadRowSegment* padRowSegment = padRow->GetPadRowSegment(j);

      if (padRowSegment->GetMotifPositionId() == motifPositionId) 
        return padRowSegment;
    }
  }
  return 0;   	
}

//______________________________________________________________________________
AliMpPadRowSegment*  
AliMpRowSegmentSpecial::FindMostRightPadRowSegment(Int_t motifPositionId) const
{
// Find the most left pad row segment with this motifPositionId.
// ---

  AliMpPadRowSegment* found = 0;

  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);    

    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) { 
      AliMpPadRowSegment* padRowSegment = padRow->GetPadRowSegment(j);

      if ( padRowSegment->GetMotifPositionId() == motifPositionId &&
           (!found || padRowSegment->RightBorderX() > found->RightBorderX()))
	   
        found = padRowSegment;  
    }
  }

  return found;   	
}

//______________________________________________________________________________
AliMpIntPair 
AliMpRowSegmentSpecial::FindRelativeLowIndicesOf(Int_t motifPositionId) const 
{ 
// Returns the lowest pad indices where the motif of the given position ID
// exists in this segment.
// ---

  AliMpIntPair ans(0,1000);
  AliMpIntPair ans0 = ans;
  Int_t maxNofPadsX=0;
  
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);

    Int_t nofPadsX=0;
    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) {
      AliMpPadRowSegment* padRowSegment = padRow->GetPadRowSegment(j);
      nofPadsX += padRowSegment->GetNofPads();
      if (padRowSegment->GetMotifPositionId() == motifPositionId) {
         if (ans.GetFirst() < nofPadsX) ans.SetFirst(nofPadsX);
         if (ans.GetSecond()>i) ans.SetSecond(i);
                  // ans.First = max (nof pads of this pos ID)
                  // ans.Second = min of pad row number
      }
    }  
    if (nofPadsX > maxNofPadsX) maxNofPadsX = nofPadsX;
  }    
  if (ans == ans0) return AliMpIntPair::Invalid();
  
  return AliMpIntPair(maxNofPadsX-ans.GetFirst(), ans.GetSecond());
}
 
//______________________________________________________________________________
TVector2 AliMpRowSegmentSpecial::MotifCenterSlow(Int_t motifPositionId) const
{
// Returns the coordinates of the motif specified with
// the given position identifier.
// !! Applicable only for motifs that have their most down pad in
// this row segment.
// ---

  // Find the first (left, down) pad row segment with this motifPositionId.
  AliMpPadRowSegment* downPadRowSegment 
    = FindPadRowSegment(motifPositionId);
  AliMpPadRowSegment* rightPadRowSegment 
    = FindMostRightPadRowSegment(motifPositionId);
  
  // Check if the motifPositionId is present 
  if (!downPadRowSegment || !rightPadRowSegment) {
    Error("MotifCenter", "Outside row segment region");
    return 0;
  }

  // Check if both pad row segments have the same motif 
  if (downPadRowSegment->GetMotif() != rightPadRowSegment->GetMotif()) {
    Fatal("MotifCenter", "Outside row segment region");
    return 0;
  }

  // Get position of found row segment
  Double_t x = rightPadRowSegment->RightBorderX();       
  Double_t y = fRow->LowBorderY()  ;   
  
  for (Int_t i=0; i<downPadRowSegment->GetPadRow()->GetID(); i++)
    y += GetPadRow(i)->HalfSizeY()*2.;
    
  // Add motifs dimensions
  x -= downPadRowSegment->GetMotif()->Dimensions().X();
  y += downPadRowSegment->GetMotif()->Dimensions().Y();
  
  return TVector2(x, y);
}

//______________________________________________________________________________
Bool_t AliMpRowSegmentSpecial::HasMotif(const AliMpVMotif* motif) const
{
// Returns true if the specified motif is already in fMotifs vector,
// returns false otherwise.
// ---

  for (UInt_t i=0; i<fMotifs.size(); i++)
    if (fMotifs[i] == motif) return true;

  return false;	 
}

//
// public methods  
//

//______________________________________________________________________________
void  AliMpRowSegmentSpecial::AddPadRow(AliMpPadRow* padRow)
{
// Adds a pad row.
// ---

  padRow->SetID(GetNofPadRows());
  padRow->SetOffsetX(fOffsetX);

  fPadRows.push_back(padRow);
}  

//______________________________________________________________________________
void AliMpRowSegmentSpecial::UpdateMotifVector()
{
// Add motifs associated with the pad row segments in the specified
// pad row in the fMotifs vector.
// ---

  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
 
    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) {
      AliMpVMotif* motif = padRow->GetPadRowSegment(j)->GetMotif();            

      if (!HasMotif(motif)) {
        fMotifs.push_back(motif);	 
        fMotifPositionIds.push_back(
          padRow->GetPadRowSegment(j)->GetMotifPositionId());
      }
    }  
  }
}

//______________________________________________________________________________
void AliMpRowSegmentSpecial::UpdatePadsOffset()
{
// Sets low indices limit to the pad offset calculated
// from the neighbour normal segment.
// ---

  // Get the neighbour row segment
  // (the first normal segment)
  AliMpVRowSegment* neighbour = fRow->GetRowSegment(1);

  // Get the the pads offset of the neighbour row segment
  // (the first normal segment)
  AliMpIntPair offset = neighbour->GetLowIndicesLimit();
  
  // Find max nof pads in a row
  Int_t maxNofPads = MaxNofPadsInRow();

  // Set limits
  SetLowIndicesLimit(offset - AliMpIntPair(maxNofPads, 0));

  // Reset limits in the neighbour row segment
  // (pad offset is now included in the special segment)  
  neighbour->SetLowIndicesLimit(
    AliMpIntPair(0, neighbour->GetLowIndicesLimit().GetSecond()));
}

//______________________________________________________________________________
Double_t  AliMpRowSegmentSpecial::LeftBorderX() const
{
// Returns the x coordinate of the left row segment border
// in global coordinate system.
// ---

  Double_t leftBorder = DBL_MAX;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
    Double_t border 
      = padRow->GetPadRowSegment(padRow->GetNofPadRowSegments()-1)->LeftBorderX();
      
    if (border < leftBorder) leftBorder =  border;
  }  
  
  return leftBorder;
}

//______________________________________________________________________________
Double_t  AliMpRowSegmentSpecial::RightBorderX() const
{
// Returns the x coordinate of the right row segment border
// in global coordinate system.
// ---

  Double_t sameBorder = fOffsetX;

  // Consistence check  
  Double_t rightBorder = -DBL_MAX;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
    Double_t border = padRow->GetPadRowSegment(0)->RightBorderX();
    if (border > rightBorder) rightBorder =  border;
  }  

  if (TMath::Abs(fOffsetX - rightBorder) > 1.e-04)  {
    Error("RightBorderX", "WrongBorder");
    return sameBorder;
  }  
  
  return rightBorder;
}

//______________________________________________________________________________
Double_t  AliMpRowSegmentSpecial::HalfSizeY() const
{
// Returns the size in y of this row segment.
// ---

  Double_t halfSizeY = 0.;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    halfSizeY += GetPadRow(i)->HalfSizeY();
  }  
  
  return halfSizeY;
}

//______________________________________________________________________________
AliMpVMotif*  AliMpRowSegmentSpecial::FindMotif(const TVector2& position) const
{
// Returns the motif of this row; 
// ---

  AliMpPadRow* padRow 
    = FindPadRow(position.Y());
  
  if (!padRow) return 0;

  AliMpPadRowSegment* padRowSegment 
    = padRow->FindPadRowSegment(position.X());
    
  if (!padRowSegment) return 0;

  return padRowSegment->GetMotif();
}  

//______________________________________________________________________________
Int_t AliMpRowSegmentSpecial::FindMotifPositionId(const TVector2& position) const
{
// Returns the motif position identified for the given
// geometric position.
// ---

  AliMpPadRow* padRow 
    = FindPadRow(position.Y());
  
  if (!padRow) return 0;

  AliMpPadRowSegment* padRowSegment 
    = padRow->FindPadRowSegment(position.X());
    
  if (!padRowSegment) return 0;

  return padRowSegment->GetMotifPositionId();
}

//______________________________________________________________________________
Bool_t AliMpRowSegmentSpecial::HasMotifPosition(Int_t motifPositionId) const
{
// Returns true if the motif specified with the given position identifier
// is in this segment.
// ---

  if (FindPadRowSegment(motifPositionId))
    return true;
  else  
    return false;   	
}

//______________________________________________________________________________
TVector2 AliMpRowSegmentSpecial::MotifCenter(Int_t motifPositionId) const
{
// Returns the coordinates of the motif specified with
// the given position identifier.
// ---

  // Try to get the motif position from the motif map first
  AliMpMotifPosition* motifPosition
    = GetRow()->GetMotifMap()->FindMotifPosition(motifPositionId);
  if (motifPosition) return motifPosition->Position();

  // Use slow method otherwise
  return MotifCenterSlow(motifPositionId);
}

//______________________________________________________________________________
TVector2 AliMpRowSegmentSpecial::Position() const
{
// Returns the position of the row segment centre.
// The centre is defined as the centre of the rectangular
// row segment envelope.
// ---

  Double_t x = fOffsetX - Dimensions().X();		    
  Double_t y = fRow->Position().Y();  
    
  return TVector2(x, y);   
}


//______________________________________________________________________________
TVector2 AliMpRowSegmentSpecial::Dimensions() const
{
// Returns the halflengths in x, y of the row segment rectangular envelope.
// ---

  Double_t x = 0.;		    
  Double_t y = 0.;  
  for (Int_t i=0; i<GetNofPadRows(); i++) {    
    AliMpPadRow* padRow = GetPadRow(i); 
    
    // Add all pad rows y halfsizes   
    y += padRow->HalfSizeY();

    // Find the biggest pad rows x halfsize
    Double_t xx 
      = (padRow->GetPadRowSegment(0)->RightBorderX() -
         padRow->GetPadRowSegment(padRow->GetNofPadRowSegments()-1)->LeftBorderX())/2.;
    if (xx > x) x = xx;		   
  }                  
    
  return TVector2(x, y);   
}

//______________________________________________________________________________
void AliMpRowSegmentSpecial::SetGlobalIndices()
{
// Sets indices limits.
// ---

  AliMpMotifPosition* firstPos = 0;
  AliMpMotifPosition* lastPos = 0;
	
  for (Int_t i=0;i<GetNofMotifs();i++) {
    AliMpMotifPosition* mPos
      = GetRow()->GetMotifMap()
        ->FindMotifPosition(GetMotifPositionId(i));
	
    if (!firstPos || 
        mPos->GetLowIndicesLimit().GetFirst()
	< firstPos->GetLowIndicesLimit().GetFirst())
       firstPos = mPos;
       
    if (!lastPos || 
        mPos->GetHighIndicesLimit().GetFirst()
	>lastPos->GetHighIndicesLimit().GetFirst())
       lastPos = mPos;
  }

  // Check if the motif positions has the limits set
  if ( !firstPos->HasValidIndices() || !lastPos->HasValidIndices())
    Fatal("SetGlobalIndices", "Indices of motif positions have to be set first.");

  SetLowIndicesLimit(firstPos->GetLowIndicesLimit());
  SetHighIndicesLimit(lastPos->GetHighIndicesLimit());
}  

//______________________________________________________________________________
Int_t AliMpRowSegmentSpecial::SetIndicesToMotifPosition(
                                          Int_t i, AliMpIntPair indices)
{
// Sets global indices to i-th motif position and returns next index in x.
// ---

  // Get motif position
  AliMpMotifPosition* motifPosition
    = GetRow()->GetMotifMap()->FindMotifPosition(GetMotifPositionId(i));

  // Low limit
  AliMpIntPair low 
    = AliMpIntPair(GetLowIndicesLimit().GetFirst() + AliMpConstants::StartPadIndex(), 
	       indices.GetSecond())
      + FindRelativeLowIndicesOf(GetMotifPositionId(i));
            
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
  
  // Increment index only if last motif position is processed 
  if (i != GetNofMotifs()-1) 
    return indices.GetFirst();
    //return 0;
  else
    return indices.GetFirst() + MaxNofPadsInRow();  
    //return MaxNofPadsInRow();  
}

//______________________________________________________________________________
AliMpRow*  AliMpRowSegmentSpecial::GetRow() const
{
// Returns the row.which this row segment belongs to.
// ---

  return fRow;
}  

//______________________________________________________________________________
Int_t  AliMpRowSegmentSpecial::GetNofMotifs() const 
{ 
// Returns the number of different motifs present in this row segment.
// ---

  return fMotifs.size();
}  

//______________________________________________________________________________
AliMpVMotif* AliMpRowSegmentSpecial::GetMotif(Int_t i) const  
{
// Returns the i-th motif present in this row segment.
// ---

   return fMotifs[i]; 
}

//______________________________________________________________________________
Int_t  AliMpRowSegmentSpecial::GetMotifPositionId(Int_t i) const 
{ 
// Returns the i-th motif position Id present in this row segment.
// ---

   return fMotifPositionIds[i]; 
} 
