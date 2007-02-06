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
// $MpId: AliMpArea.cxx,v 1.8 2006/05/24 13:58:29 ivana Exp $
// Category: basic
//
// Class AliMpArea
// ----------------
// Class that defines a rectangle area positioned in plane..
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpArea.h"
#include "AliMpConstants.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpArea)
/// \endcond

//_____________________________________________________________________________
AliMpArea::AliMpArea(const TVector2& position, const TVector2& dimensions)
  : TObject(),
    fPosition(position),
    fDimensions(dimensions),
    fValidity(true) 
{
/// Standard constructor

  // Check dimensions
  if (  fDimensions.X() < - AliMpConstants::LengthTolerance() || 
        fDimensions.Y() < - AliMpConstants::LengthTolerance() || 
      ( fDimensions.X() < AliMpConstants::LengthTolerance() && 
        fDimensions.Y() < AliMpConstants::LengthTolerance() ) )
  {
    fDimensions = TVector2();
    fValidity = false;
  }  
}

//_____________________________________________________________________________
AliMpArea::AliMpArea()
  : TObject(),
    fPosition(TVector2()),
    fDimensions(TVector2()), 
    fValidity(false) 
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpArea::AliMpArea(const AliMpArea& rhs):
  TObject(rhs),
  fPosition(rhs.fPosition),
  fDimensions(rhs.fDimensions), 
  fValidity(rhs.fValidity) 
{
/// Copy constructor
}

//_____________________________________________________________________________
AliMpArea::~AliMpArea() 
{
/// Destructor
}

//
// operators
//

//______________________________________________________________________________
AliMpArea& AliMpArea::operator = (const AliMpArea& right)
{
/// Assignment operator

  // check assignment to self
  if (this == &right) return *this;

  // base class assignment
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
/// Return the position of the left edge.

  return fPosition.X() - fDimensions.X();
}

//_____________________________________________________________________________
Double_t AliMpArea::RightBorder() const
{
/// Return the position of right edge.

  return fPosition.X() + fDimensions.X();
}

//_____________________________________________________________________________
Double_t AliMpArea::UpBorder() const
{
/// Return the position of the up edge.

  return fPosition.Y() + fDimensions.Y();
}

//_____________________________________________________________________________
Double_t AliMpArea::DownBorder() const
{
/// Return the position of the down edge.

  return fPosition.Y() - fDimensions.Y();
}

//_____________________________________________________________________________
TVector2 AliMpArea::LeftDownCorner() const
{
/// Return position of the left down corner.

  return TVector2(LeftBorder(), DownBorder());
}  

//_____________________________________________________________________________
TVector2 AliMpArea::LeftUpCorner() const
{
/// Return position of the left up corner.

  return TVector2(LeftBorder(), UpBorder());
}  

//_____________________________________________________________________________
TVector2 AliMpArea::RightDownCorner() const
{
/// Return position of the right down corner.

  return TVector2(RightBorder(), DownBorder());
}  


//_____________________________________________________________________________
TVector2 AliMpArea::RightUpCorner() const
{
/// Return position of the right up corner.

  return TVector2(RightBorder(), UpBorder());
}  

//_____________________________________________________________________________
void
AliMpArea::Print(Option_t*) const
{
  cout << (*this) << endl;
}

//_____________________________________________________________________________
ostream& operator<< (ostream &stream,const AliMpArea& area)
{
/// Output streaming

  stream << "Area: position: (" 
         << area.Position().X() << ", " << area.Position().Y() << ") " 
	 << " dimensions: (" 
         << area.Dimensions().X() << ", " << area.Dimensions().Y() << ") " 
  << " valid: " << (area.IsValid()==true ? "YES":"NO")
	 << endl;
  return stream;
}

