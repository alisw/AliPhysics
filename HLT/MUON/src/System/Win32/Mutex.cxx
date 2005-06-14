////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include <windows.h>
#include <winbase.h>

#include "System/Mutex.hpp"

namespace dHLT
{
namespace System
{


Mutex::Mutex()
{
	InitializeCriticalSection(&mutex);
};


Mutex::~Mutex()
{
	DeleteCriticalSection(&mutex);
};


void Mutex::Lock()
{
	EnterCriticalSection(&mutex);
};


bool Mutex::TryLock()
{
	//return TryEnterCriticalSection(&mutex);
	EnterCriticalSection(&mutex);
	return true;
};


void Mutex::Unlock()
{
	LeaveCriticalSection(&mutex);
};


} // System
} // dHLT

