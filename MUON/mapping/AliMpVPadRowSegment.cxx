// $Id$
// Category: sector
//
// Class AliMpVPadRowSegment
// --------------------
// The abstract base class for a pad row segment composed of the 
// the identic pads.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpVPadRowSegment.h"
#include "AliMpPadRow.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"

ClassImp(AliMpVPadRowSegment)

//_____________________________________________________________________________
AliMpVPadRowSegment::AliMpVPadRowSegment(AliMpPadRow* padRow, AliMpMotif* motif, 
                                       Int_t motifPositionId, Int_t nofPads)
  : TObject(),
    fNofPads(nofPads),
    fOffsetX(0.),
    fPadRow(padRow),
    fMotif(motif),
    fMotifPositionId(motifPositionId)
{
// 
}

//_____________________________________________________________________________
AliMpVPadRowSegment::AliMpVPadRowSegment() 
  : TObject(),
    fNofPads(0),
    fOffsetX(0.),
    fPadRow(0),
    fMotif(0),
    fMotifPositionId(0)
{
//
}

//_____________________________________________________________________________
AliMpVPadRowSegment::AliMpVPadRowSegment(const AliMpVPadRowSegment& right) 
  : TObject(right) {
// 
  Fatal("AliMpVPadRowSegment", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpVPadRowSegment::~AliMpVPadRowSegment() {
//  
}

//
// operators
//

//_____________________________________________________________________________
AliMpVPadRowSegment& 
AliMpVPadRowSegment::operator=(const AliMpVPadRowSegment& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//
// public methods  
//

//_____________________________________________________________________________
Double_t  AliMpVPadRowSegment::HalfSizeY() const
{
// Returns the size in y of this row segment.
// ---

  return fMotif->GetPadDimensions().Y();
}

//_____________________________________________________________________________
AliMpPadRow*  AliMpVPadRowSegment::GetPadRow() const
{
// Returns the pad row.which this pad row segment belongs to.
// ---

  return fPadRow;
}  

//_____________________________________________________________________________
AliMpMotif*  AliMpVPadRowSegment::GetMotif() const
{
// Returns the motif of this pad row segment. 
// ---

  return fMotif;
}  

//_____________________________________________________________________________
Int_t  AliMpVPadRowSegment::GetMotifPositionId() const
{
// Returns the motif of this pad row segment. 
// ---

  return fMotifPositionId;
}  

//_____________________________________________________________________________
void  AliMpVPadRowSegment::SetOffsetX(Double_t offsetX)
{
// Sets the x offset.
// ---

  fOffsetX = offsetX;
}    

