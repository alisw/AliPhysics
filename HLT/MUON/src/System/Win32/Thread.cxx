////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include <windows.h>

#include "System/Thread.hpp"
#include "System/SystemError.hpp"

namespace
{

DWORD WINAPI ThreadStartRoutine(LPVOID thispointer)
{
	using dHLT::System::Thread;
	Thread* thread = static_cast<Thread*>( thispointer );
	return thread->Execute();
};

} // anonymous namespace


namespace dHLT
{
namespace System
{


Thread::Thread()
{
	thread = CreateThread(NULL, 0, &ThreadStartRoutine, (LPVOID)this, CREATE_SUSPENDED, NULL);
	if (thread == NULL) throw System::Error();
};


void Thread::Start()
{
	if ( ResumeThread(thread) == -1 ) throw System::Error();
};


Thread::~Thread()
{
	CloseHandle(thread);
};


Int Thread::Priority() const
{
	return GetThreadPriority(thread);
};


void Thread::Priority(const Int value)
{
	if ( SetThreadPriority(thread, value) != TRUE )
		throw System::Error();
};


Int Thread::WaitForThread() const
{
	DWORD result = WaitForSingleObject(thread, INFINITE);
	if (result == WAIT_FAILED) throw System::Error();
	DWORD exitcode;
	if ( GetExitCodeThread(thread, &exitcode) != TRUE )
		throw System::Error();
	return exitcode;
};


void Thread::Sleep(const UInt milliseconds) const
{
	Sleep(milliseconds);
};


void Thread::SwitchContext() const
{
	Sleep(0);
};


void Thread::Exit(const Int exitcode)
{
	ExitThread(exitcode);
};


} // System
} // dHLT

