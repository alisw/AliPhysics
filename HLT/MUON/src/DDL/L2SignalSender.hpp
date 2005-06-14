////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DDL_L2_SIGNAL_SENDER_HPP
#define dHLT_DDL_L2_SIGNAL_SENDER_HPP

#include "BCMP/Sender.hpp"
#include "EventID.hpp"

namespace dHLT
{
namespace DDL
{


class L2SignalSender : public BCMP::EventHandler
{
public:

	L2SignalSender(const UShort port = 4900);
	virtual ~L2SignalSender();

	void Terminate() { terminate = true; };
	
	void SignalL2(const EventID eventid);
	
	void Run();
	
private:

	bool terminate;
	BCMP::Sender sender;
};


}; // DDL
}; // dHLT

#endif // dHLT_DDL_L2_SIGNAL_SENDER_HPP

