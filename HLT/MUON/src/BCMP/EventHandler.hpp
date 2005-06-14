////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BCMP_EVENT_HANDLER_HPP
#define dHLT_BCMP_EVENT_HANDLER_HPP

#include "BCMP/EventQueue.hpp"

#include <iostream>
using std::endl;
using std::cerr;

namespace dHLT
{
namespace BCMP
{


class EventHandler
{
public:
	
	/* Called whenever a new sender or receiver has connected.
	 */
	virtual void OnConnect(const System::Address& address) {};
	
	/* Called whenever a sender or receiver has disconnected.
	 */
	virtual void OnDisconnect(const System::Address& address) {};
	
	/* Called whenever a connection timed out and was lost.
	 */
	virtual void OnConnectionLost(const System::Address& address) {};
	
	/* Called when a message has been received.
	 */
	virtual void OnMessage(
			const char* message, const UInt length,
			const System::Address& from
		) {};
	
	/* Called when a protocol error message has been received.
	   The default behaviour is to write the message to standard output.
	 */
	virtual void OnError(
			const char* message, const Int errorcode,
			const System::Address& from
		)
	{
		cerr << "Error: [" << errorcode << "] from " << from << endl;
		cerr << "       " << message << endl;
	};
};


} // BCMP
} // dHLT

#endif // dHLT_BCMP_EVENT_HANDLER_HPP
