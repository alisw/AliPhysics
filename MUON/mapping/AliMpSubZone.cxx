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
// $MpId: AliMpSubZone.cxx,v 1.6 2005/08/26 15:43:36 ivana Exp $
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
/// Standard constructor
}

//_____________________________________________________________________________
AliMpSubZone::AliMpSubZone() 
  : TObject(),
    fMotif(0)
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpSubZone::AliMpSubZone(const AliMpSubZone& right) 
  : TObject(right) 
{
/// Protected copy constructor (not provided) 

  Fatal("AliMpSubZone", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpSubZone::~AliMpSubZone() 
{
// Destructor 
}

//
// operators
//

//_____________________________________________________________________________
AliMpSubZone& AliMpSubZone::operator=(const AliMpSubZone& right)
{
/// Protected assignment operator (not provided)

  // check assignment to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignment operator not provided.");
    
  return *this;  
}    

//
// public methods
//

//_____________________________________________________________________________
void AliMpSubZone::AddRowSegment(AliMpVRowSegment* rowSegment)
{
/// Add row segment.

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
/// Print motif position Ids for all row segments.
 
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
/// Return number of row segments.

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
/// Return i-th row segment.

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
/// Return the motif.

  return fMotif;
}  
