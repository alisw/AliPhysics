////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/Region.hpp"
#include "AliRoot/Point.hpp"
#include <TMath.h>
#include "AliHLTMUONUtils.h"

ClassImp(AliHLTMUONRegion)


AliHLTMUONRegion::AliHLTMUONRegion() :
	TObject(), fLeft(0), fRight(0), fBottom(0), fTop(0)
{
	fLeft = fRight = fBottom = fTop = 0.0;
}


AliHLTMUONRegion::AliHLTMUONRegion(Float_t left, Float_t right, Float_t bottom, Float_t top)
	 : TObject(), fLeft(0), fRight(0), fBottom(0), fTop(0)
{
// Construct a region of interest from the boundary box defined by the borders
// 'left', 'right', 'top' and 'bottom'.

	if (left > right)
	{
		fLeft = fRight = fBottom = fTop = 0.0;
		Error("AliHLTMUONRegion", "parameter left (%f) is larger than right (%f).", left, right);
	}
	else if (bottom > top)
	{
		fLeft = fRight = fBottom = fTop = 0.0;
		Error("AliHLTMUONRegion", "parameter bottom (%f) is larger than top (%f).", bottom, top);
	}
	else
	{
		fLeft = left;
		fRight = right;
		fBottom = bottom;
		fTop = top;
	}
}


void AliHLTMUONRegion::Left(Float_t value)
{
// Set the left border of the ROI.

	if (value > fRight)
		Error("Left", "Trying to assign fLeft (%f) larger than fRight (%f).", value, fRight);
	else
		fLeft = value;
}


void AliHLTMUONRegion::Right(Float_t value)
{
// Set the right border of the ROI.

	if (value < fLeft)
		Error("Right", "Trying to assign fRight (%f) smaller than fLeft (%f).", value, fLeft);
	else
		fRight = value;
}


void AliHLTMUONRegion::Bottom(Float_t value)
{
// Set the bottom border of the ROI.

	if (value > fTop)
		Error("Bottom", "Trying to assign fBottom (%f) larger than fTop (%f).", value, fTop);
	else
		fBottom = value;
}


void AliHLTMUONRegion::Top(Float_t value)
{
// Set the top border of the ROI.

	if (value < fBottom)
		Error("Top", "Trying to assign fTop (%f) smaller than fBottom (%f).", value, fBottom);
	else
		fTop = value;
}


Bool_t AliHLTMUONRegion::Contains(const AliHLTMUONPoint& p) const
{
// Checks if the point is within this region. If it is then kTRUE is returned
// otherwise kFALSE is returned.

	return 
	  fLeft <= p.X()
	  && p.X() <= fRight 
	  && fBottom <= p.Y()
	  && p.Y() <= fTop;
}


ostream& operator << (ostream& os, const AliHLTMUONRegion& r)
{
	os << "[(" << r.fLeft << ", " << r.fRight << "), (" << r.fLeft << ", " << r.fRight << ")]";
	return os;
}

