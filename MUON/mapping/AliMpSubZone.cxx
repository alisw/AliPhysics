// $Id$
// Category: sector
//
// Class AliMpSubZone
// ------------------
// Class describing a zone segment composed of the 
// line segments with the same motif type.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>
#include <TError.h>

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
AliMpSubZone::AliMpSubZone(const AliMpSubZone& right) 
  : TObject(right) {
// 
  Fatal("AliMpSubZone", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpSubZone::~AliMpSubZone() {
//  
}

//
// operators
//

//_____________________________________________________________________________
AliMpSubZone& AliMpSubZone::operator=(const AliMpSubZone& right)
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
void AliMpSubZone::AddRowSegment(AliMpVRowSegment* rowSegment)
{
// Adds row segment.
// ---

#ifdef WITH_STL
  fSegments.push_back(rowSegment);
#endif

#ifdef WITH_ROOT
  fSegments.Add(rowSegment);
#endif
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
// ---

#ifdef WITH_STL
  return fSegments.size();
#endif

#ifdef WITH_ROOT
  return fSegments.GetSize();
#endif
}  

//_____________________________________________________________________________
AliMpVRowSegment* AliMpSubZone::GetRowSegment(Int_t i) const 
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

//_____________________________________________________________________________
AliMpVMotif*  AliMpSubZone:: GetMotif() const
{
// Returns the motif.
// ---

  return fMotif;
}  
