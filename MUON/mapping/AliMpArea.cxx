// $Id$
// Category: basic
//
// Class AliMpArea
// ----------------
// Class that defines a rectangle area positioned in plane..
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>

#include "AliMpArea.h"

ClassImp(AliMpArea)


//_____________________________________________________________________________
AliMpArea::AliMpArea(const TVector2& position, const TVector2& dimensions)
  : TObject(),
    fPosition(position),
    fDimensions(dimensions),
    fValidity(true) {
//
  // Check dimensions
  if (fDimensions.X() <= 0. || fDimensions.Y() <=0.) {
    fDimensions = TVector2();
    fValidity = false;
  }  
}

//_____________________________________________________________________________
AliMpArea::AliMpArea()
  : TObject(),
    fPosition(TVector2()),
    fDimensions(TVector2()), 
    fValidity(false) {
//
}

//_____________________________________________________________________________
AliMpArea::AliMpArea(const AliMpArea& rhs):
  TObject(rhs),
  fPosition(rhs.fPosition),
  fDimensions(rhs.fDimensions) {
//
}

//_____________________________________________________________________________
AliMpArea::~AliMpArea() {
//
}

//
// operators
//

//______________________________________________________________________________
AliMpArea& AliMpArea::operator = (const AliMpArea& right)
{
// Assignement operator

  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  TObject::operator=(right);

  fPosition = right.fPosition;
  fDimensions = right.fDimensions;
  fValidity = right.fValidity;

  return *this;
} 

//
// public methods
//

//_____________________________________________________________________________
Double_t AliMpArea::LeftBorder() const
{
// Returns the position of the left edge.
// --

  return fPosition.X() - fDimensions.X();
}

//_____________________________________________________________________________
Double_t AliMpArea::RightBorder() const
{
// Returns the position of right edge.
// --

  return fPosition.X() + fDimensions.X();
}

//_____________________________________________________________________________
Double_t AliMpArea::UpBorder() const
{
// Returns the position of the up edge.
// --

  return fPosition.Y() + fDimensions.Y();
}

//_____________________________________________________________________________
Double_t AliMpArea::DownBorder() const
{
// Returns the position of the down edge.
// --

  return fPosition.Y() - fDimensions.Y();
}

//_____________________________________________________________________________
TVector2 AliMpArea::LeftDownCorner() const
{
// Returns position of the left down corner.
// --

  return TVector2(LeftBorder(), DownBorder());
}  

//_____________________________________________________________________________
TVector2 AliMpArea::LeftUpCorner() const
{
// Returns position of the left up corner.
// --

  return TVector2(LeftBorder(), UpBorder());
}  

//_____________________________________________________________________________
TVector2 AliMpArea::RightDownCorner() const
{
// Returns position of the right down corner.
// --

  return TVector2(RightBorder(), DownBorder());
}  


//_____________________________________________________________________________
TVector2 AliMpArea::RightUpCorner() const
{
// Returns position of the right up corner.
// --

  return TVector2(RightBorder(), UpBorder());
}  

//_____________________________________________________________________________
ostream& operator<< (ostream &stream,const AliMpArea& area)
{
  stream << "Area: position: (" 
         << area.Position().X() << ", " << area.Position().Y() << ") " 
	 << " dimensions: (" 
         << area.Dimensions().X() << ", " << area.Dimensions().Y() << ") " 
	 << endl;
  return stream;
}

