// $Id$
// Category: sector
//
// Class AliMpRow
// --------------
// Class describing a row composed of the row segments.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>
#include <TError.h>
#include <TMath.h>

#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpVRowSegmentSpecial.h"
#include "AliMpRowSegmentRSpecial.h"
#include "AliMpVMotif.h"
#include "AliMpMotifType.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifMap.h"
#include "AliMpConstants.h"

ClassImp(AliMpRow)

//_____________________________________________________________________________
AliMpRow::AliMpRow(Int_t id, AliMpMotifMap* motifMap) 
  : AliMpVIndexed(),
    fID(id),
    fOffsetY(0.),
    fSegments(),
    fMotifMap(motifMap)
{
//
}

//_____________________________________________________________________________
AliMpRow::AliMpRow() 
  : AliMpVIndexed(),
    fID(0),
    fOffsetY(0.),
    fSegments(),
    fMotifMap(0)
{
//
}

//_____________________________________________________________________________
AliMpRow::AliMpRow(const AliMpRow& right) 
  : AliMpVIndexed(right) {
// 
  Fatal("AliMpRow", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpRow::~AliMpRow() {
// 

#ifdef WITH_STL
  for (Int_t i=0; i<GetNofRowSegments(); i++)
    delete fSegments[i]; 
#endif

#ifdef WITH_ROOT
  fSegments.Delete();
#endif
}

//
// operators
//

//_____________________________________________________________________________
AliMpRow&  AliMpRow::operator=(const AliMpRow& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//
// private methods
//

//_____________________________________________________________________________
AliMpVRowSegment*  AliMpRow::FindRowSegment(Int_t ix) const
{    
// Finds first normal row segment with low indices limit >= ix.
// --- 

  for (Int_t i=0; i<GetNofRowSegments(); i++) {
    AliMpVRowSegment* segment = GetRowSegment(i);

    if (!dynamic_cast<AliMpVRowSegmentSpecial*>(segment) &&
         segment->GetHighIndicesLimit().GetFirst() >= ix)
	 
     return segment;	 
  }   

  return 0;	 
}

//_____________________________________________________________________________
AliMpMotifPosition*  
AliMpRow::FindMotifPosition(AliMpVRowSegment* segment, Int_t ix) const
{
// Finds first motif position in the specified row segment 
// with high indices limit >= ix.
// --- 

  if (!segment) return 0;

  for (Int_t i=0; i<segment->GetNofMotifs(); i++){
     AliMpMotifPosition* motifPosition 
       = GetMotifMap()->FindMotifPosition(segment->GetMotifPositionId(i));
       
     if(!motifPosition) {
       Fatal("FindMotifPosition", "Not found.");
       return 0;
     }  
     
     if (motifPosition->GetHighIndicesLimit().GetFirst()>=ix) 
       return motifPosition;
  }
  
  return 0;     
}


//_____________________________________________________________________________
void AliMpRow::SetHighIndicesLimits(Int_t iy)
{
// Sets the global indices high limit to its row segments,
// motif positions with a given value.
// Keeps ix unmodified.
// --- 

  for (Int_t j=0; j<GetNofRowSegments(); j++) {
     AliMpVRowSegment* rowSegment = GetRowSegment(j);       
     rowSegment
       ->SetHighIndicesLimit(
	    AliMpIntPair(rowSegment->GetHighIndicesLimit().GetFirst(),iy));

    for (Int_t k=0; k<rowSegment->GetNofMotifs(); k++) {

      Int_t motifPositionId = rowSegment->GetMotifPositionId(k);
      AliMpMotifPosition* motifPosition 
	= GetMotifMap()->FindMotifPosition(motifPositionId);

      motifPosition
	->SetHighIndicesLimit(
	     AliMpIntPair(motifPosition->GetHighIndicesLimit().GetFirst(), iy));
     
    }
  }  
}

//_____________________________________________________________________________
void  AliMpRow::CheckEmpty() const
{
// Give a fatal if row is empty.
// ---

  if (GetNofRowSegments() == 0) 
    Fatal("CheckEmpty", "Empty row");
}

//
// public methods
//

//_____________________________________________________________________________
void AliMpRow::AddRowSegment(AliMpVRowSegment* rowSegment)
{
// Adds row segment at the end.
// ---

#ifdef WITH_STL
  fSegments.push_back(rowSegment);
#endif

#ifdef WITH_ROOT
  fSegments.Add(rowSegment);
#endif
}  
  
//_____________________________________________________________________________
void AliMpRow::AddRowSegmentInFront(AliMpVRowSegment* rowSegment)
{
// Inserts row segment in the first vector position.
// ---

#ifdef WITH_STL
  fSegments.insert(fSegments.begin(), rowSegment);
#endif

#ifdef WITH_ROOT
  fSegments.AddFirst(rowSegment);
#endif
}  
  
//_____________________________________________________________________________
AliMpVRowSegment* AliMpRow::FindRowSegment(Double_t x) const
{
// Finds the row segment for the specified x position;
// returns 0 if no row segment is found.
// ---

  for (Int_t i=0; i<GetNofRowSegments(); i++) {

#ifdef WITH_STL
    AliMpVRowSegment* rs = fSegments[i];
#endif
#ifdef WITH_ROOT
    AliMpVRowSegment* rs = (AliMpVRowSegment*)fSegments.At(i);
#endif

    if (x >= rs->LeftBorderX() && x <= rs->RightBorderX())
      return rs;
  }
  
  return 0;    
}    

//_____________________________________________________________________________
Double_t AliMpRow::LowBorderY() const
{
// Returns the lowest row offset (the Y coordinate of the position of the
// low border of motif).
// ---

  CheckEmpty();

  return fOffsetY - GetRowSegment(0)->HalfSizeY();
}  

//_____________________________________________________________________________
Double_t AliMpRow::UpperBorderY() const
{
// Returns the uppermost row offset (the Y coordinate of the position of the
// upper border of motif).
// ---

  CheckEmpty();

  return fOffsetY + GetRowSegment(0)->HalfSizeY();
}  

//_____________________________________________________________________________
AliMpVPadIterator* AliMpRow::CreateIterator() const
{
// Iterator is not yet implemented.
// ---

  Fatal("CreateIterator", "Iterator is not yet implemented.");
  
  return 0;
}  

//_____________________________________________________________________________
void AliMpRow::SetMotifPositions()
{
// Creates motif positions objects and fills them in the motif map.
// ---

  CheckEmpty();

  for (Int_t j=0; j<GetNofRowSegments(); j++) {
     AliMpVRowSegment* rowSegment = GetRowSegment(j);

     for (Int_t k=0; k<rowSegment->GetNofMotifs(); k++) {
        // Get values 
	Int_t motifPositionId = rowSegment->GetMotifPositionId(k);
	AliMpVMotif* motif = rowSegment->GetMotif(k);
	TVector2 position = rowSegment->MotifCenter(motifPositionId);
       
        AliMpMotifPosition* motifPosition 
	  = new AliMpMotifPosition(motifPositionId, motif, position);
        // set the initial value to of HighIndicesLimit() Invalid()
        // (this is used for calculation of indices in case of
        // special row segments)
        motifPosition->SetHighIndicesLimit(AliMpIntPair::Invalid());

        //Bool_t warn = (rowSegment->GetNofMotifs()==1); 
        Bool_t warn = true;
	if (dynamic_cast<AliMpVRowSegmentSpecial*>(rowSegment)) warn = false; 
               // supress warnings for special row segments
	       // which motifs can overlap the row borders
	       
        Bool_t added = GetMotifMap()->AddMotifPosition(motifPosition, warn);
	
	if (!added) delete motifPosition;	
     }  
  }
}    

//_____________________________________________________________________________
void AliMpRow::SetGlobalIndices(AliMpDirection constPadSizeDirection, 
                                AliMpRow* rowBefore)
{
// Sets the global indices limits to its row segments,
// motif positions.
// ---

  Int_t ix = AliMpConstants::StartPadIndex();
  Int_t iy = AliMpConstants::StartPadIndex();

  for (Int_t j=0; j<GetNofRowSegments(); j++) {
     AliMpVRowSegment* rowSegment = GetRowSegment(j);
     
     ix += rowSegment->GetLowIndicesLimit().GetFirst();

     for (Int_t k=0; k<rowSegment->GetNofMotifs(); k++) {
     
       // Find the y index value of the low edge
       if (rowBefore) {
         if (constPadSizeDirection == kY) {
           iy = rowBefore->GetHighIndicesLimit().GetSecond()+1;
         } 
	 else {
           AliMpVRowSegment* seg = rowBefore->FindRowSegment(ix);	
	   AliMpMotifPosition* motPos =  FindMotifPosition(seg, ix);
	   if (!dynamic_cast<AliMpRowSegmentRSpecial*>(rowSegment)) {
             if (!motPos) 
	       Fatal("SetGlobalIndices", "Motif position in rowBefore not found.");
	   
             iy = motPos->GetHighIndicesLimit().GetSecond()+1;
	   }  
         }
       } 

       // Set (ix, iy) to k-th motif position and update ix
       ix = rowSegment->SetIndicesToMotifPosition(k, AliMpIntPair(ix, iy));		   
    }
    rowSegment->SetGlobalIndices();    
  }

  SetLowIndicesLimit(GetRowSegment(0)->GetLowIndicesLimit());
  SetHighIndicesLimit(GetRowSegment(GetNofRowSegments()-1)->GetHighIndicesLimit());

  return ;
}

//_____________________________________________________________________________
TVector2  AliMpRow::Position() const
{
// Returns the position of the row centre.
// ---

  Double_t x = (GetRowSegment(0)->LeftBorderX() +
                GetRowSegment(GetNofRowSegments()-1)->RightBorderX())/2.;
		    
  Double_t y = fOffsetY;  
    
  return TVector2(x, y);   
}

//_____________________________________________________________________________
TVector2  AliMpRow::Dimensions() const
{
// Returns the maximum halflengths of the row in x, y.
// ---

  Double_t x = (GetRowSegment(GetNofRowSegments()-1)->RightBorderX() -
                GetRowSegment(0)->LeftBorderX())/2.;
                  
  Double_t y = GetRowSegment(0)->HalfSizeY();  
    
  return TVector2(x, y);   
}

//_____________________________________________________________________________
void AliMpRow::SetRowSegmentOffsets(const TVector2& offset)
{
// Sets the row segments offsets in X .
// ---

  CheckEmpty();
  
  AliMpVRowSegment* previous = 0;

  for (Int_t j=0; j<GetNofRowSegments(); j++) {
     AliMpVRowSegment* rowSegment = GetRowSegment(j);

     Double_t offsetX;
     if (previous) 
      offsetX = previous->RightBorderX();
    else
      offsetX = offset.X();  
  
    rowSegment->SetOffset(TVector2(offsetX, 0.));
    previous = rowSegment;  
  }
}


//_____________________________________________________________________________
Double_t AliMpRow::SetOffsetY(Double_t offsetY)
{
// Sets the row offset (the Y coordinate of the position of the
// center of motif) and returns the offset of the top border.
// ---

  CheckEmpty();

  AliMpVRowSegment* first = GetRowSegment(0);
  Double_t rowSizeY = first->HalfSizeY();
  
  // Check if all next row segments have motif of
  // the same size in y
  for (Int_t i=1; i<GetNofRowSegments(); i++) {
     Double_t sizeY = GetRowSegment(i)->HalfSizeY();
     
     if (TMath::Abs(sizeY - rowSizeY) >= AliMpConstants::LengthTolerance()) {
       //cout << GetID() << "th row " << i << "th segment " 
       //     << sizeY << "  " << rowSizeY  << endl;
       Fatal("SetOffsetY", "Motif with different Y size in one row");
       return 0.;
     }  
  }

  offsetY += rowSizeY ;
    
  fOffsetY = offsetY;
    
  return offsetY += rowSizeY;
}  

//_____________________________________________________________________________
Int_t AliMpRow::GetNofRowSegments() const 
{
// Returns number of row segments.
// ---

#ifdef WITH_STL
  return fSegments.size();
#endif

#ifdef WITH_ROOT
  return fSegments.GetSize();
#endif
}  

//_____________________________________________________________________________
AliMpVRowSegment* AliMpRow::GetRowSegment(Int_t i) const 
{
// Returns i-th row segment.
// ---

  if (i<0 || i>=GetNofRowSegments()) {
    Warning("GetRowSegment", "Index outside range");
    return 0;
  }
  
#ifdef WITH_STL
  return fSegments[i];  
#endif

#ifdef WITH_ROOT
  return (AliMpVRowSegment*)fSegments.At(i);  
#endif
}
 
