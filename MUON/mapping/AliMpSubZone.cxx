// $Id$
// Category: sector
//
// Class AliMpSubZone
// ------------------
// Class describing a zone segment composed of the 
// line segments with the same motif type.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>

#include "AliMpSubZone.h"
#include "AliMpVRowSegment.h"
#include "AliMpVMotif.h"

ClassImp(AliMpSubZone)

//_____________________________________________________________________________
AliMpSubZone::AliMpSubZone(AliMpVMotif* motif) 
  : TObject(),
    fMotif(motif)
{
//
}

//_____________________________________________________________________________
AliMpSubZone::AliMpSubZone() 
  : TObject(),
    fMotif(0)
{
//
}

//_____________________________________________________________________________
AliMpSubZone::~AliMpSubZone() {
//  
}

//_____________________________________________________________________________
void AliMpSubZone::AddRowSegment(AliMpVRowSegment* rowSegment)
{
// Adds row segment.
// ---

  fSegments.push_back(rowSegment);
} 


//_____________________________________________________________________________
void AliMpSubZone::Print(const char* /*option*/) const
{
// Prints motif position Ids for all row segments.
// --
 
  for (Int_t i=0; i<GetNofRowSegments(); i++) {
    AliMpVRowSegment* rowSegment = GetRowSegment(i);
    
    cout << rowSegment->GetNofMotifs() << " ";

    for (Int_t j=0; j<rowSegment->GetNofMotifs(); j++)
      cout << rowSegment->GetMotifPositionId(j) << " ";
    
    cout << endl;    
  }    
}
  
//_____________________________________________________________________________
Int_t AliMpSubZone::GetNofRowSegments() const 
{
// Returns number of row segments.

  return fSegments.size();
}  

//_____________________________________________________________________________
AliMpVRowSegment* AliMpSubZone::GetRowSegment(Int_t i) const 
{
  if (i<0 || i>=GetNofRowSegments()) {
    Warning("GetRowSegment", "Index outside range");
    return 0;
  }
  
  return fSegments[i];  
}

//_____________________________________________________________________________
AliMpVMotif*  AliMpSubZone:: GetMotif() const
{
// Returns the motif.
// ---

  return fMotif;
}  
