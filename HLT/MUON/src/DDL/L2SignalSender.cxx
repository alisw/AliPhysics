////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "DDL/L2SignalSender.hpp"
#include "Error.hpp"
#include "Utils.hpp"
#include <asm/errno.h>

#include <string.h>
#include <iostream>

namespace dHLT
{
namespace DDL
{


L2SignalSender::L2SignalSender(const UShort port)
	: sender(port)
{
	DebugMsg(9, "L2SignalSender::L2SignalSender(" << port << ")");
	sender.eventhandler = this;
};


L2SignalSender::~L2SignalSender()
{
	sender.Terminate();
};


void L2SignalSender::SignalL2(const EventID eventid)
{
	char buffer[4 + sizeof(EventID)];
	UInt length = 4 + sizeof(EventID);
	char* message = &buffer[0];
	strcpy(message, "L2OK");
	memcpy(message + 4, &eventid, sizeof(EventID));
	sender.Send(message, length);
};


void L2SignalSender::Run()
{
	// Enter an event handling loop:
	while (not terminate)
	{
		try
		{
			// Try handle events or timeout every 50 milliseconds
			// to check if we were signaled to terminate or not. 
			sender.HandleEvents(50);
		}
		catch (System::Error& e)
		{
			if (e.ErrorCode() != EINTR) throw;
		};
	};
};


}; // DDL
}; // dHLT
