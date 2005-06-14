////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/Thread.hpp"
#include "System/SystemError.hpp"

#include <pthread.h>
#include <unistd.h>

#include <iostream>
using namespace std;

namespace
{

void* ThreadStartRoutine(void* thispointer)
{
	using dHLT::System::Thread;
	Thread* thread = static_cast<Thread*>( thispointer );
	return (void*) thread->Execute();
};

} // end of namespace


namespace dHLT
{
namespace System
{


Thread::Thread()
{
	// Nothing to be done.
};


void Thread::Start()
{
	int result = pthread_create(&thread, NULL, &ThreadStartRoutine, (void*)this);
	if (result != 0) throw System::Error(result);
};


Thread::~Thread()
{
	pthread_cancel(thread);
};


Int Thread::Priority() const
{
	struct sched_param param;
	int policy;
	int result = pthread_getschedparam(thread, &policy, &param);
	if (result != 0) throw System::Error(result);
	return param.sched_priority;
};


void Thread::Priority(const Int value)
{
	struct sched_param param;
	param.sched_priority = value;
	int result = pthread_setschedparam(thread, SCHED_OTHER, &param);
	if (result != 0) throw System::Error(result);
};


Int Thread::WaitForThread() const
{
	void* returnvalue;
	int result = pthread_join(thread, &returnvalue);
	if (result != 0) throw System::Error(result);
	return (Int) returnvalue;
};


void Thread::Sleep(const UInt milliseconds) const
{
	struct timespec ts;
	ts.tv_sec = milliseconds / 1000;
	ts.tv_nsec = (milliseconds % 1000) * 1000000;
	int result = nanosleep(&ts, NULL);
	if (result != 0) throw System::Error();
};


void Thread::SwitchContext() const
{
	if (sched_yield() != 0) throw System::Error();
};


void Thread::Exit(const Int exitcode)
{
	pthread_exit( (void*)exitcode );
};


} // System
} // dHLT

