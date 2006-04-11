////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCOREEVENTID_H
#define ALIHLTMUONCOREEVENTID_H

#include "BasicTypes.hpp"
#include "Utils.hpp"


// Must correct this definition!!!!
struct AliHLTMUONCoreEventID
{
	UInt fBunch;
	UInt fTimeStamp;


	AliHLTMUONCoreEventID(UInt bunch = 0, UInt timestamp = 0)
	{
		fBunch = bunch;
		fTimeStamp = timestamp;
	};


	bool operator == (const AliHLTMUONCoreEventID& rhs)
	{
		return fBunch == rhs.fBunch and fTimeStamp == rhs.fTimeStamp;
	};

	bool operator != (const AliHLTMUONCoreEventID& rhs)
	{
		return not (*this == rhs);
	};
};


#endif // ALIHLTMUONCOREEVENTID_H
