////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/Mutex.hpp"
#include "System/SystemError.hpp"

#include <pthread.h>
#include <asm/errno.h>

namespace dHLT
{
namespace System
{


Mutex::Mutex()
{
	pthread_mutex_init(&mutex, NULL);
};


Mutex::~Mutex()
{
	int result = pthread_mutex_destroy(&mutex);
	if (result != 0) throw System::Error(result);
};


void Mutex::Lock()
{
	pthread_mutex_lock(&mutex);
};


bool Mutex::TryLock()
{
	if ( pthread_mutex_trylock(&mutex) == EBUSY )
		return false;
	else
		return true;
};


void Mutex::Unlock()
{
	pthread_mutex_unlock(&mutex);
};


} // System
} // dHLT

