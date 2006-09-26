////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/Point.hpp"
#include "TMath.h"

ClassImp(AliHLTMUONPoint)


AliHLTMUONPoint::AliHLTMUONPoint() : TObject(), fX(0), fY(0)
{
	fX = fY = 0.0;
}

AliHLTMUONPoint::AliHLTMUONPoint(Float_t x, Float_t y) : TObject(), fX(x), fY(y)
{
	fX = x;
	fY = y;
}

ostream& operator << (ostream& os, const AliHLTMUONPoint& p)
{
	os << "[" << p.fX << ", " << p.fY << "]";
	return os;
}
