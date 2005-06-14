////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_MUTEX_CONDITION_HPP
#define dHLT_SYSTEM_MUTEX_CONDITION_HPP

#include "BasicTypes.hpp"
#include "System/SystemTypes.hpp"
#include "System/Mutex.hpp"

namespace dHLT
{
namespace System
{


class MutexCondition
{
public:

	MutexCondition();
	~MutexCondition();

	void Wait(Mutex& mutex);

	// timeout in milliseconds.
	bool Wait(Mutex& mutex, const UInt timeout);

	void Signal(bool broadcast = true);

private:

	SystemMutexCondition condition;

};


} // System
} // dHLT

#endif // dHLT_SYSTEM_MUTEX_CONDITION_HPP


