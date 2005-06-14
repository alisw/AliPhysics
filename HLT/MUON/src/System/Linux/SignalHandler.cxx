////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/SignalHandler.hpp"
#include "System/SystemError.hpp"
#include "Utils.hpp"

#include <signal.h>
#include <vector>

namespace
{
	// The currenly installed signal handlers.
	std::vector<dHLT::System::SignalHandler*> handlers;
	
	
	// Global system signal handler routine to be installed by SignalHandler class.
	void sighandler(int signalnumber, siginfo_t*, void*)
	{
		std::vector<dHLT::System::SignalHandler*>::iterator h = handlers.begin();
		while (h != handlers.end())
		{
			if ( (*h)->SignalNumber() == signalnumber )
				(*h)->HandleSignal();
			h++;
		};
	};
	
}; // end of namespace


namespace dHLT
{
namespace System
{


SignalHandler::SignalHandler(int signalnumber)
{
	// Add this hanlder to the handlers list and remember our signal number.
	handlers.push_back(this);
	signum = signalnumber;
	
	// Replace the system signal handler.
	struct sigaction sa;
	sa.sa_sigaction = &sighandler;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_SIGINFO;
	int result = sigaction(signalnumber, &sa, NULL);
	if (result != 0) throw System::Error();
};


SignalHandler::~SignalHandler()
{
	// Remove handler.
	std::vector<SignalHandler*>::iterator i = handlers.begin();
	while (i != handlers.end())
	{
		if ( *i == this )
		{
			handlers.erase(i);
			break;
		};
		i++;
	};
};


} // System
} // dHLT
