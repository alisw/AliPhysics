// $Id$
// --------------------------------------------------------
// Category: sector
//
// Class AliMpPadRow
// -----------------
// Class describing a pad row composed of the pad row segments.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpPadRow.h"
#include "AliMpPadRowSegment.h"

ClassImp(AliMpPadRow)

//_____________________________________________________________________________
AliMpPadRow::AliMpPadRow() 
  : TObject(),
    fID(0) 
{
//
}

//_____________________________________________________________________________
AliMpPadRow::~AliMpPadRow() {
//  

  for (Int_t i=0; i<GetNofPadRowSegments() ; i++)
    delete fSegments[i];
}

//
// public methods
//

//_____________________________________________________________________________
void AliMpPadRow::AddPadRowSegment(AliMpPadRowSegment* padRowSegment)
{
// Adds row segment.
// ---

  // Set pad row segment offset
  if (GetNofPadRowSegments() == 0)
    padRowSegment
      ->SetOffsetX(fOffsetX);
  else 
    padRowSegment
      ->SetOffsetX(GetPadRowSegment(GetNofPadRowSegments()-1)->LeftBorderX());

  // Adds the pad row segment
  fSegments.push_back(padRowSegment);
}  
  
//_____________________________________________________________________________
AliMpPadRowSegment* AliMpPadRow::FindPadRowSegment(Double_t x) const
{
// Finds the row segment for the specified x position;
// returns 0 if no row segment is found.
// ---

  for (Int_t i=0; i<GetNofPadRowSegments(); i++) {
    AliMpPadRowSegment* rs = GetPadRowSegment(i);
    if (x >= rs->LeftBorderX() && x <= rs->RightBorderX())
      return rs;
  }
  
  return 0;    
}    

//_____________________________________________________________________________
Double_t  AliMpPadRow::HalfSizeY() const
{
  return GetPadRowSegment(0)->HalfSizeY();
}

//_____________________________________________________________________________
void  AliMpPadRow::SetID(Int_t id)
{
// Sets the ID.
// ---

  fID = id;
}    

//_____________________________________________________________________________
void  AliMpPadRow::SetOffsetX(Double_t offsetX)
{
// Sets the x offset.
// ---

  fOffsetX = offsetX;
}    

//_____________________________________________________________________________
Int_t AliMpPadRow::GetID() const 
{
// Returns the row ID.
// ---

  return fID;
}  

//_____________________________________________________________________________
Int_t AliMpPadRow::GetNofPadRowSegments() const 
{
// Returns number of row segments.
// ---

  return fSegments.size();
}  

//_____________________________________________________________________________
AliMpPadRowSegment* AliMpPadRow::GetPadRowSegment(Int_t i) const 
{
  if (i<0 || i>=GetNofPadRowSegments()) {
    Warning("GetRowSegment", "Index outside range");
    return 0;
  }
  
  return fSegments[i];  
}

//_____________________________________________________________________________
Int_t AliMpPadRow::GetNofPads() const 
{
// Returns number of pads in this pad row.
// ---

  Int_t nofPads=0;
  for (Int_t i=0; i<GetNofPadRowSegments(); i++)
    nofPads += GetPadRowSegment(i)->GetNofPads();

  return nofPads;
}  

