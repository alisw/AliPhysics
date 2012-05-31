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

/// \class AliMUONSegment
///
/// A basic line segment, to be used in contour making algorithms. 
///
/// In particular, this class defines what a left or right edge is.
///
/// Also, please note that, due to the way Root collections are sorted (relying
/// on TObject::Compare method), the way the AliMUONSegment::Compare method
/// is implemented below is really important when it comes to understand
/// contour making algorithm. Keep that in mind. 
///
/// \author Laurent Aphecetche, Subatech
///

#include "AliMUONSegment.h"

#include "TMath.h"
#include "Riostream.h"
#include "AliMpConstants.h"

/// \cond CLASSIMP
ClassImp(AliMUONSegment)
/// \endcond

const Double_t AliMUONSegment::fgkPrecision(AliMpConstants::LengthTolerance());
  
//_____________________________________________________________________________
AliMUONSegment::AliMUONSegment() : 
TObject(),
fStartX(), fStartY(), fEndX(), fEndY(), fSmallerY(), fIsHorizontal(), fIsVertical(),
fIsLeftEdge(), fIsRightEdge(), fIsAPoint(kTRUE)
{
  /// Ctor
  Set(fStartX,fStartY,fEndX,fEndY);
}

//_____________________________________________________________________________
AliMUONSegment::AliMUONSegment(Double_t xstart, Double_t ystart, Double_t xend, Double_t yend)
: TObject(),
fStartX(xstart), fStartY(ystart), fEndX(xend), fEndY(yend), fSmallerY(), fIsHorizontal(), fIsVertical(),
fIsLeftEdge(), fIsRightEdge(), fIsAPoint(kTRUE)
{
  /// Ctor
  Set(xstart,ystart,xend,yend);
}

//_____________________________________________________________________________
Bool_t
AliMUONSegment::AreEqual(double a, double b)
{
  /// Whether the two floats are equal within the given precision
  return (TMath::Abs(b-a) < fgkPrecision);
}

//_____________________________________________________________________________
Int_t	
AliMUONSegment::Compare(const TObject* obj) const
{
  /// Compare method, which sort segments in ascending x order
  /// if same x, insure that left edges are before right edges
  /// within same x, order by increasing bottommost y
  /// Mind your steps ! This method is critical to the contour merging algorithm !
  
  const AliMUONSegment* rhs = static_cast<const AliMUONSegment*>(obj);
  
  if ( AreEqual(StartX(),rhs->StartX()) )
  {
    if ( IsLeftEdge() && rhs->IsRightEdge() ) return -1;
    if ( IsRightEdge() && rhs->IsLeftEdge() ) return 1;
    if ( SmallerY() < rhs->SmallerY() ) return -1;
    if ( SmallerY() > rhs->SmallerY() ) return 1;
    return 0;
  }
  else if ( StartX() < rhs->StartX() )
  {
    return -1;
  }
  else //if ( StartX() > rhs->StartX() ) 
  {
    return 1;
  }
}

//_____________________________________________________________________________
double AliMUONSegment::Top() const 
{
  /// Max Y of the segment
  return TMath::Max(fStartY,fEndY); 
}

//_____________________________________________________________________________
double AliMUONSegment::Distance() const 
{
  /// Length of the segment
  return TMath::Sqrt((fStartX-fEndX)*(fStartX-fEndX) +
                     (fStartY-fEndY)*(fStartY-fEndY)); 
}

//_____________________________________________________________________________
void AliMUONSegment::Print(Option_t*) const
{
  /// Printout
  cout << AsString() << endl;
}

//_____________________________________________________________________________
const char* AliMUONSegment::AsString() const 
{
  /// Return a string representation of this object
  return Form("[ (%10.5f,%10.5f) -> (%10.5f,%10.5f) %s ] (d=%e)",fStartX,fStartY,fEndX,fEndY,
              IsLeftEdge() ? "L" : ( IsRightEdge() ? "R" : ( IsHorizontal() ? "H" : "" )),
              Distance() ); 
}  

//_____________________________________________________________________________
void 
AliMUONSegment::Set(Double_t xstart, Double_t ystart, Double_t xend, Double_t yend)
{
  /// Set start and end point, and (re)compute internal values
  fStartX = xstart;
  fEndX = xend;
  fStartY = ystart;
  fEndY = yend;
  fSmallerY = TMath::Min(fStartY,fEndY); 
  fIsHorizontal = AreEqual(fStartY,fEndY); 
  fIsVertical = AreEqual(fStartX,fEndX); 
  fIsLeftEdge = fIsVertical && ( fStartY > fEndY );
  fIsRightEdge = fIsVertical && ( fStartY < fEndY );
  fIsAPoint = ( Distance() < fgkPrecision );
}

