////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/Point.hpp"
#include "TMath.h"

ClassImp(AliMUONHLT::Point);

namespace AliMUONHLT
{

Point::Point() : TObject()
{
	fX = fY = 0.0;
};

Point::Point(const Float_t x, const Float_t y) : TObject()
{
	fX = x;
	fY = y;
};

ostream& operator << (ostream& os, const Point& p)
{
	os << "[" << p.fX << ", " << p.fY << "]";
	return os;
};

}; // AliMUONHLT

