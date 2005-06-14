////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/MutexCondition.hpp"
#include "System/SystemError.hpp"

#include <sys/time.h>
#include <pthread.h>
#include <asm/errno.h>

namespace dHLT
{
namespace System
{


MutexCondition::MutexCondition()
{
	int result = pthread_cond_init(&condition, NULL);
	if (result != 0) throw System::Error(result);
};


MutexCondition::~MutexCondition()
{
	int result = pthread_cond_destroy(&condition);
	if (result != 0) throw System::Error(result);
};


void MutexCondition::Wait(Mutex& mutex)
{
	int result = pthread_cond_wait(&condition, &mutex.mutex);
	if (result != 0) throw System::Error(result);
};


bool MutexCondition::Wait(Mutex& mutex, const UInt timeout)
{
	struct timespec abstime;
	struct timeval now;

	if ( gettimeofday(&now, NULL) != 0 ) throw System::Error();
	abstime.tv_sec = timeout / 1000 + now.tv_sec;
	abstime.tv_nsec = (timeout % 1000) * 1000000 + now.tv_usec * 1000;

	int result = pthread_cond_timedwait(&condition, &mutex.mutex, &abstime);
	if (result == ETIMEDOUT)
		return false;
	else if (result != 0)
		throw System::Error(result);
	else
		return true;
};


void MutexCondition::Signal(bool broadcast)
{
	int result;
	if (broadcast)
		result = pthread_cond_broadcast(&condition);
	else
		result = pthread_cond_signal(&condition);
	if (result != 0) throw System::Error(result);
};


} // System
} // dHLT

