////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include <windows.h>

#include "System/MutexCondition.hpp"

namespace dHLT
{
namespace System
{


MutexCondition::MutexCondition()
{
	condition = CreateEvent(NULL, FALSE, FALSE, NULL);
	if (condition == NULL) throw System::Error();
};


MutexCondition::~MutexCondition()
{
	CloseHandle(condition);
};


void MutexCondition::Wait(Mutex& mutex)
{
	mutex.Unlock();
	DWORD result = WaitForSingleObject(condition, INFINITE);
	mutex.Lock();
	if (result == WAIT_FAILED) throw System::Error();
};


bool MutexCondition::Wait(Mutex& mutex, const UInt timeout)
{
	mutex.Unlock();
	DWORD result = WaitForSingleObject(condition, timeout);
	mutex.Lock();
	if (result == WAIT_FAILED) throw System::Error();
	if (result == WAIT_TIMEOUT)
		return false;
	else
		return true;
};


void MutexCondition::Signal(bool broadcast)
{
	if ( PulseEvent(condition) != TRUE )
		throw System::Error();
};


} // System
} // dHLT

