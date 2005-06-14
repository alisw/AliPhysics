////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_THREAD_HPP
#define dHLT_SYSTEM_THREAD_HPP

#include "BasicTypes.hpp"
#include "System/SystemTypes.hpp"

namespace dHLT
{
namespace System
{


class Thread
{
public:

	Thread();
	virtual ~Thread();
	
	// Start() should only be called once.
	void Start();

	Int Priority() const;

	void Priority(const Int value);

	Int WaitForThread() const;

	virtual Int Execute() = 0;

protected:

	void Sleep(const UInt milliseconds) const;

	void SwitchContext() const;

	void Exit(const Int exitcode);

private:

	SystemThread thread;
};


} // System
} // dHLT

#endif // dHLT_SYSTEM_THREAD_HPP


