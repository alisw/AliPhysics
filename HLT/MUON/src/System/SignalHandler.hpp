////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_SIGNAL_HANDLER_HPP
#define dHLT_SYSTEM_SIGNAL_HANDLER_HPP

namespace dHLT
{
namespace System
{


/* Define a macro for quick implementation of signal handlers.
   Example:
   
      HandleSignals(SIGINT,
          cout << "got SIGINT" << endl;
      );
 */
#define HandleSignals(signum, code) \
	class HandleSignals_##signum : public dHLT::System::SignalHandler \
	{ \
	public: \
		HandleSignals_##signum() : dHLT::System::SignalHandler(signum) {}; \
		virtual void HandleSignal() { code }; \
	} handlesignals_##signum;


class SignalHandler
{
public:

	SignalHandler(int signalnumber);
	virtual ~SignalHandler();
	
	virtual void HandleSignal() = 0;
	
	int SignalNumber() { return signum; };
	
private:

	int signum;
};


} // System
} // dHLT

#endif // dHLT_SYSTEM_SIGNAL_HANDLER_HPP
