////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <asm/errno.h>

using std::endl;
using std::cout;
using std::cerr;

#include "System/SystemError.hpp"
#include "System/Socket.hpp"
#include "BCMP/Receiver.hpp"
#include "DDL/FileList.hpp"
#include "System/SignalHandler.hpp"
#include "Error.hpp"
#include "Utils.hpp"

using namespace dHLT;
using DDL::FileList;


bool terminate = false;

HandleSignals(SIGINT,
	terminate = true;
);


class Handler : public BCMP::EventHandler
{
public:
	
	virtual void OnConnect(const System::Address& address)
	{
		cout << "Connected to: " << address << endl;
	};
	
	virtual void OnDisconnect(const System::Address& address)
	{
		cout << "Disonnected to: " << address << endl;
	};
	
	virtual void OnConnectionLost(const System::Address& address)
	{
		cout << "Connection to: " << address << " lost." << endl;
	};
	
	virtual void OnMessage(
			const char* message, const UInt length,
			const System::Address& from
		)
	{
		cout << "Got message of " << length << " bytes from: " << from << endl;
		char* str = new char[length+1];
		memcpy(str, message, length);
		str[length] = '\0';
		cout << "    : " << str << endl;
		delete [] str;
	};
};

Handler handler;



void PrintUsage()
{
	cerr << "Usage: > TriggerDDLFilePublisher [-r] <filelist> [<port>]" << endl;
	cerr << "        -r : Optional flag indicating to recursively go through sub directories" << endl;
	cerr << "             and add files to the file list for publishing." << endl;
	cerr << "<filelist> : Required file name containing a file list (one file/directory per line)" << endl;
	cerr << "             of files to publish." << endl;
	cerr << "    <port> : Optional port number to listen to for L2 signals." << endl;
};


/* Parse command line and fill the port number and filelist.
   If the port number or filelist is untouched then the default
   values are used.
   Returns true if the command line is parsed properly.
 */
bool ParseCommandLine(const int argc, const char** argv, UShort& port, FileList& filelist)
{
	// TODO: parse command line.
	// command line format:
	// > TriggerDDLFilePublisher [-r] <filelist> [<port>]
	//         -r : Optional flag indicating to recursively go through sub directories
	//              and add files to the file list for publishing.
	// <filelist> : The file name containing a file list (one file/directory per line)
	//              of files to publish.
	//     <port> : The port number to listen to for L2 signals.
	
	PrintUsage();
	return false;
};


int main(const int argc, const char** argv)
{
	try
	{
		UShort port = 4900;
		FileList filelist;

		bool valid_arguments = ParseCommandLine(argc, argv, port, filelist);
		if (not valid_arguments) return 1;
		
		BCMP::Receiver receiver(&handler, port);
		cout << "Listening on: " << receiver.LocalAddress() << endl;
		
		// Enter an event handling loop:
		while (not terminate)
		{
			try
			{
				// Try handle events or timeout every 50 milliseconds
				// to check if we were signaled to terminate or not. 
				receiver.HandleEvents(10);
			}
			catch (System::Error& e)
			{
				if (e.ErrorCode() != EINTR) throw;
			};
		};
	}
	catch (const Error& e)
	{
		cerr << "Error [" << e.ErrorCode() << "]: " << e << endl;
	};
	
	cout << "done." << endl;

	return 0;
};
