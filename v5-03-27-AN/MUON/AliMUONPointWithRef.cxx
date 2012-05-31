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

///\class AliMUONPointWithRef
///
/// A class used to represent a point with an external integer reference
/// and with a specific sorting method (see AliMUONContourMaker)
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliMUONPointWithRef.h"

#include "AliMUONSegment.h"
#include "Riostream.h"
#include "TString.h"

//\cond CLASSIMP
ClassImp(AliMUONPointWithRef)
//\endcond

//_____________________________________________________________________________
AliMUONPointWithRef::AliMUONPointWithRef() : fX(), fY(), fRef(-1)
{
  /// default ctor
}

//_____________________________________________________________________________
AliMUONPointWithRef::AliMUONPointWithRef(Double_t x, Double_t y, Int_t ref)
: fX(x), fY(y), fRef(ref)
{
  /// ctor
}

//_____________________________________________________________________________
Int_t	
AliMUONPointWithRef::Compare(const TObject* obj) const
{
  /// Should serve to sort the vertical edges in ascending order, first on absissa, 
  /// then on ordinate
  
  if ( this == obj ) return 0;
  
  const AliMUONPointWithRef* rhs = static_cast<const AliMUONPointWithRef*>(obj);

  if ( AliMUONSegment::AreEqual(Y(),rhs->Y()) )
  {
    if ( AliMUONSegment::AreEqual(X(),rhs->X()) )
    {
      return 0;
    }
    else if ( X() > rhs->X() )
    {
      return 1;
    }
    else 
      return -1;
  }
  else if ( Y() < rhs->Y() )
  {
    return -1;
  }
  else
  {
    return 1;
  }
}

//_____________________________________________________________________________
void 
AliMUONPointWithRef::Print(Option_t*) const
{
  /// Printout
  cout << Form("(%10.5f,%10.5f) [%4d]",X(),Y(),Ref()) << endl;
}
