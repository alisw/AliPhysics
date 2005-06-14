////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "DDL/L2SignalReceiver.hpp"
#include "Error.hpp"

#include <string.h>
#include <iostream>
#include <asm/errno.h>

namespace dHLT
{
namespace DDL
{


L2SignalReceiver::L2SignalReceiver(const UShort port)
	: receiver(port)
{
	receiver.eventhandler = this;
};


L2SignalReceiver::~L2SignalReceiver()
{
	receiver.Terminate();
};


void L2SignalReceiver::Run()
{
	// Enter an event handling loop:
	while (not terminate)
	{
		try
		{
			// Try handle events or timeout every 50 milliseconds
			// to check if we were signaled to terminate or not. 
			receiver.HandleEvents(50);
		}
		catch (System::Error& e)
		{
			if (e.ErrorCode() != EINTR) throw;
		};
	};
};


void L2SignalReceiver::OnMessage(
		const char* message, const UInt length, const System::Address& from
	)
{
	UInt* x = (UInt*) message;
	char* msg = "L2OK";
	// Check that the size is 8 bytes and the first 4 bytes contain the chars L2OK.
	if ( (length == 4 + sizeof(EventID)) and (*((UInt*)msg) == *x) )
	{
		GotL2( *((EventID*)(message+4)) );
	};
};


}; // DDL
}; // dHLT
