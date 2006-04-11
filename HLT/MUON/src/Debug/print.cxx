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


ostream& operator << (ostream& os, const AliHLTMUONCoreEventID& id)
{
	os << "<" << id.fBunch << ":" << id.fTimeStamp << ">";
	return os;
}


ostream& operator << (ostream& os, const AliHLTMUONCorePoint& p)
{
	os << "[" << p.fX << ", " << p.fY << "]";
	return os;
}


ostream& operator << (ostream& os, const AliHLTMUONCoreParticleSign s)
{
	switch (s)
	{
	case kSignMinus:   os << "Minus";   break;
	case kSignPlus:    os << "Plus";    break;
	case kUnknownSign: os << "Unknown"; break;
	default:           os << "FAULT!";
	}
	return os;
}


ostream& operator << (ostream& os, const AliHLTMUONCoreTriggerRecord& rec)
{
	char* signstr;
	switch (rec.fSign)
	{
	case kSignMinus:   signstr = "Minus  "; break;
	case kSignPlus:    signstr = "Plus   "; break;
	case kUnknownSign: signstr = "Unknown"; break;
	default:           signstr = "FAULT!!";
	}
	
	os << "{ sign = " << signstr << ", pt = " 
	   << rec.fPt << ", impact on station 1: "
	   << rec.fStation1impact << " , station 2: "
	   << rec.fStation2impact << " }";
	return os;
}


ostream& operator << (ostream& os, const AliHLTMUONCoreChamberID chamber)
{
	Int ch = (Int)chamber;
	if (0 <= ch && ch <= 9)
		os << ch + 1;
	else
		os << "FAULT!!";
	return os;
}

