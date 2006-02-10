////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Debug/print.hpp"

#include <iostream>
using std::endl;
using std::cout;
using std::ostream;

namespace dHLT
{

ostream& operator << (ostream& os, const EventID& id)
{
	os << "<" << id.bunch << ":" << id.timestamp << ">";
	return os;
}


ostream& operator << (ostream& os, const Point& p)
{
	os << "[" << p.x << ", " << p.y << "]";
	return os;
}


ostream& operator << (ostream& os, const ParticleSign s)
{
	switch (s)
	{
	case Minus:       os << "Minus";   break;
	case Plus:        os << "Plus";    break;
	case UnknownSign: os << "Unknown"; break;
	default:          os << "FAULT!";
	}
	return os;
}


ostream& operator << (ostream& os, const TriggerRecord& rec)
{
	char* signstr;
	switch (rec.sign)
	{
	case Minus:       signstr = "Minus  "; break;
	case Plus:        signstr = "Plus   "; break;
	case UnknownSign: signstr = "Unknown"; break;
	default:          signstr = "FAULT!!";
	}
	
	os << "{ sign = " << signstr << ", pt = " 
	   << rec.pt << ", impact on station 1: "
	   << rec.station1impact << " , station 2: "
	   << rec.station2impact << " }";
	return os;
}


ostream& operator << (ostream& os, const ChamberID chamber)
{
	Int ch = (Int)chamber;
	if (0 <= ch && ch <= 9)
		os << ch + 1;
	else
		os << "FAULT!!";
	return os;
}


} // dHLT
