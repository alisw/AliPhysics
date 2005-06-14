////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_EVENT_HPP
#define dHLT_SYSTEM_EVENT_HPP

#include "BasicTypes.hpp"
#include "SystemTypes.hpp"

namespace dHLT
{
namespace System
{


class EventSignal
{
public:

	EventSignal();
	~EventSignal();

	void Lock();
	bool TryLock();
	void Unlock();

private:

	SystemEvent mutex;

};


} // System
} // dHLT

#endif // dHLT_SYSTEM_EVENT_HPP


