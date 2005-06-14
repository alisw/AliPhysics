////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_EVENT_ID_HPP
#define dHLT_EVENT_ID_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"

namespace dHLT
{


// Must correct this definition!!!!
struct EventID
{
	UInt bunch;
	UInt timestamp;


	EventID(UInt bunch = 0, UInt timestamp = 0)
	{
		this->bunch = bunch;
		this->timestamp = timestamp;
	};


	bool operator == (const EventID& rhs)
	{
		return bunch == rhs.bunch and timestamp == rhs.timestamp;
	};

	bool operator != (const EventID& rhs)
	{
		return not (*this == rhs);
	};
};


} // dHLT

#endif // dHLT_EVENT_ID_HPP
