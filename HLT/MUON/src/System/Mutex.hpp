////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_MUTEX_HPP
#define dHLT_SYSTEM_MUTEX_HPP

#include "BasicTypes.hpp"
#include "System/SystemTypes.hpp"

namespace dHLT
{
namespace System
{


class MutexCondition;


class Mutex
{
public:

	Mutex();
	~Mutex();

	void Lock();
	bool TryLock();
	void Unlock();

private:

	friend class MutexCondition;
	SystemMutex mutex;

};


} // System
} // dHLT

#endif // dHLT_SYSTEM_MUTEX_HPP


