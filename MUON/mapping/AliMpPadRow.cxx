// $Id$
// --------------------------------------------------------
// Category: sector
//
// Class AliMpPadRow
// ------------------
// Class describing a pad row composed of the pad row segments.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpPadRow.h"
#include "AliMpPadRowLSegment.h"
#include "AliMpPadRowRSegment.h"

ClassImp(AliMpPadRow)

//_____________________________________________________________________________
AliMpPadRow::AliMpPadRow(AliMpXDirection direction) 
  : TObject(),
    fDirection(direction), 
    fID(0)
{
//
}

//_____________________________________________________________________________
AliMpPadRow::AliMpPadRow() 
  : TObject(),
    fDirection(kLeft), 
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
// private methods
//

//_____________________________________________________________________________
Double_t AliMpPadRow::CurrentBorderX() const
{
// Returns the left/right x border 
// (depending on the direction which the row segments are filled in).
// ---

  if (GetNofPadRowSegments() == 0)
      return fOffsetX;
  else 
    if (fDirection == kLeft)
      return GetPadRowSegment(GetNofPadRowSegments()-1)->LeftBorderX();
    else  
      return GetPadRowSegment(GetNofPadRowSegments()-1)->RightBorderX();
}

//
// public methods
//

//_____________________________________________________________________________
AliMpVPadRowSegment* 
AliMpPadRow::AddPadRowSegment(AliMpMotif* motif, Int_t motifPositionId,
                              Int_t nofPads)
{
// Adds pad row segment.
// ---

  AliMpVPadRowSegment* padRowSegment = 0;

  if (fDirection == kLeft) {
    padRowSegment 
      = new AliMpPadRowLSegment(this, motif, motifPositionId, nofPads);
  }    
  else  {
    padRowSegment 
      = new AliMpPadRowRSegment(this, motif, motifPositionId, nofPads);
  }     

  // Set pad row segment offset
  padRowSegment->SetOffsetX(CurrentBorderX());

  // Adds the pad row segment
#ifdef WITH_STL
  fSegments.push_back(padRowSegment);
#endif

#ifdef WITH_ROOT
  fSegments.Add(padRowSegment);
#endif
  
  return padRowSegment;
}  
  
//_____________________________________________________________________________
AliMpVPadRowSegment* AliMpPadRow::FindPadRowSegment(Double_t x) const
{
// Finds the row segment for the specified x position;
// returns 0 if no row segment is found.
// ---

  for (Int_t i=0; i<GetNofPadRowSegments(); i++) {
    AliMpVPadRowSegment* rs = GetPadRowSegment(i);
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

#ifdef WITH_STL
  return fSegments.size();
#endif

#ifdef WITH_ROOT
  return fSegments.GetEntriesFast();
#endif
}  

//_____________________________________________________________________________
AliMpVPadRowSegment* AliMpPadRow::GetPadRowSegment(Int_t i) const 
{
// Returns pad row segment with specified number.
// ---

  if (i<0 || i>=GetNofPadRowSegments()) {
    Warning("GetRowSegment", "Index outside range");
    return 0;
  }
  
#ifdef WITH_STL
  return fSegments[i];  
#endif

#ifdef WITH_ROOT
  return (AliMpVPadRowSegment*)fSegments[i];  
#endif
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

