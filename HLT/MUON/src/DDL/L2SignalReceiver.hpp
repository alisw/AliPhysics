////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DDL_L2_SIGNAL_RECEIVER_HPP
#define dHLT_DDL_L2_SIGNAL_RECEIVER_HPP

#include "BCMP/Receiver.hpp"
#include "EventID.hpp"

namespace dHLT
{
namespace DDL
{


class L2SignalReceiver : public BCMP::EventHandler
{
public:

	L2SignalReceiver(const UShort port = 4900);
	virtual ~L2SignalReceiver();

	void Terminate() { terminate = true; };
	
	void Run();
	
	virtual void GotL2(const EventID eventid) = 0;
	
	virtual void OnMessage(
			const char* message, const UInt length,
			const System::Address& from
		);
	
private:

	bool terminate;
	BCMP::Receiver receiver;
};


}; // DDL
}; // dHLT

#endif // dHLT_DDL_L2_SIGNAL_RECEIVER_HPP

