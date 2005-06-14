////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
// Author: Gareth de Vaux
// Email:  devaux@lhc.phy.uct.ac.za | dhlt@lordcow.org
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <signal.h>
#include <asm/errno.h>

using std::endl;
using std::cout;
using std::cerr;

#include "System/SystemError.hpp"
#include "System/Socket.hpp"
#include "DDL/L2SignalReceiver.hpp"
#include "DDL/FileList.hpp"
#include "System/SignalHandler.hpp"
#include "Debug/print.hpp"
#include "Error.hpp"
#include "Utils.hpp"

using namespace dHLT;
using DDL::FileList;


class L2SignalReceiver : public DDL::L2SignalReceiver
{
public:

	L2SignalReceiver(FileList& list, const UShort port = 4900)
		: DDL::L2SignalReceiver(port)
	{
		filelist = &list;
	};

	virtual void GotL2(const EventID eventid)
	{
		cout << "Got signal: " << eventid << endl;
	};
	
private:

	FileList* filelist;
};


L2SignalReceiver* signalreceiver = NULL;

// If we get a SIGINT ( keyboard Ctrl+C ) signal or SIGHUP then terminate the program.
HandleSignals(SIGINT,
	if (signalreceiver != NULL) signalreceiver->Terminate();
);
HandleSignals(SIGHUP,
	if (signalreceiver != NULL) signalreceiver->Terminate();
);


void PrintUsage()
{
	cerr << "Usage: > triggerDDL [-r] <filelist> [<port>]" << endl;
	cerr << "        -r : Optional flag indicating to recursively go through sub directories" << endl;
	cerr << "             and add files to the file list for publishing." << endl;
	cerr << "<filelist> : Required file name containing a file list (one file/directory per line)" << endl;
	cerr << "             of files to publish." << endl;
	cerr << "    <port> : Optional port number to listen to for L2 signals." << endl;
};


bool CheckIfEmpty(const FileList& filelist)
{
	if (filelist.Empty())
	{
		cerr << "Error: No files added to file list. Is the path empty?" << endl;
		return false;
	}
	else
		return true;
};


/* Parse command line and fill the port number and filelist.
   If the port number or filelist is untouched then the default
   values are used.
   Returns true if the command line is parsed properly.
 */
bool ParseCommandLine(const int argc, const char** argv, UShort& port, FileList& filelist)
{
	for (int i = 1; i < argc; i++)
	{
		if (	strcmp(argv[i], "-h") == 0 or strcmp(argv[i], "-help") == 0
			or strcmp(argv[i], "--help") == 0)
		{
			PrintUsage();
			return false;
		};
	};

	if (argc == 2)
	{
		filelist = argv[1];
		return CheckIfEmpty(filelist);
	}

	if (argc == 3)
	{
		if (strcmp(argv[1],"-r") == 0)
		{
			filelist.Add(argv[2], true);
			return CheckIfEmpty(filelist);
		}
		else
		{
			int x = atoi(argv[2]);
			if ( not (0 < x and x < 65536) )
			{
				cerr << "Error: Port is out of range, valid range is [1..65535]." << endl;
				return false;
			}
			else
			{
				port = x;
				filelist.Add(argv[1], false);
				return CheckIfEmpty(filelist);
			}
		}
	}

	if (argc == 4)
	{
		if (strcmp(argv[1],"-r") != 0)
		{
			cerr << "Error: Unknown flag '" << argv[1] << "', expected '-r'." << endl;
			return false;
		};
		
		int x = atoi(argv[3]);
		if ( not (0 < x and x < 65536) )
		{
			cerr << "Error: Port is out of range, valid range is [1..65535]." << endl;
			return false;
		}
		else
		{
			port = x;
			filelist.Add(argv[2], true);
			return CheckIfEmpty(filelist);
		}
	}

	PrintUsage();
	return false;
};


int main(const int argc, const char** argv)
{
	int returncode = 0;
	try
	{
		UShort port = 4900;
		FileList filelist;

		bool valid_arguments = ParseCommandLine(argc, argv, port, filelist);
		if (valid_arguments)
		{
			L2SignalReceiver receiver(filelist, port);
			signalreceiver = &receiver;
			receiver.Run();
		}
		else
			returncode = 1;
	}
	catch (const System::Error& e)
	{
		cerr << "Error [" << e.ErrorCode();
		if (e.ErrorCode() == EACCES)
			cerr << "]: Could not allocate socket. ";
		else
			cerr << "]: ";
		cerr << e << endl;
		returncode = 2;
	}
	catch (const Error& e)
	{
		cerr << "Error [" << e.ErrorCode() << "]: " << e << endl;
		returncode = 2;
	}
	catch (...)
	{
		cerr << "Error: Unknown exception!" << endl;
		returncode = 2;
	};

	return returncode;
};

