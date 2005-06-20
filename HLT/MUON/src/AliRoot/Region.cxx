////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/Region.hpp"
#include "AliRoot/Point.hpp"
#include <TMath.h>
#include "Utils.hpp"

ClassImp(AliMUONHLT::Region);

namespace AliMUONHLT
{


Region::Region() : TObject()
{
	fLeft = fRight = fBottom = fTop = 0.0;
};


Region::Region(const Float_t left, const Float_t right, const Float_t bottom, const Float_t top)
	 : TObject()
{
	if (left > right)
	{
		fLeft = fRight = fBottom = fTop = 0.0;
		Error("Region", "parameter left (%f) is larger than right (%f).", left, right);
	}
	else if (bottom > top)
	{
		fLeft = fRight = fBottom = fTop = 0.0;
		Error("Region", "parameter bottom (%f) is larger than top (%f).", bottom, top);
	}
	else
	{
		fLeft = left;
		fRight = right;
		fBottom = bottom;
		fTop = top;
	};
};


void Region::Left(const Float_t value)
{
	if (value > fRight)
		Error("Left", "Trying to assign fLeft (%f) larger than fRight (%f).", value, fRight);
	else
		fLeft = value;
};


void Region::Right(const Float_t value)
{
	if (value < fLeft)
		Error("Right", "Trying to assign fRight (%f) smaller than fLeft (%f).", value, fLeft);
	else
		fRight = value;
};


void Region::Bottom(const Float_t value)
{
	if (value > fTop)
		Error("Bottom", "Trying to assign fBottom (%f) larger than fTop (%f).", value, fTop);
	else
		fBottom = value;
};


void Region::Top(const Float_t value)
{
	if (value < fBottom)
		Error("Top", "Trying to assign fTop (%f) smaller than fBottom (%f).", value, fBottom);
	else
		fTop = value;
};


Bool_t Region::Contains(const Point& p) const
{
	return fLeft <= p.fX and p.fX <= fRight and fBottom <= p.fY and p.fY <= fTop;
};


ostream& operator << (ostream& os, const Region& r)
{
	os << "[(" << r.fLeft << ", " << r.fRight << "), (" << r.fLeft << ", " << r.fRight << ")]";
	return os;
};


}; // AliMUONHLT
