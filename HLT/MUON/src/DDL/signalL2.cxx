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
#include <netdb.h>
#include <unistd.h>
#include <signal.h>

using std::endl;
using std::cout;
using std::cerr;

#include "System/SystemError.hpp"
#include "System/Socket.hpp"
#include "System/Thread.hpp"
#include "Error.hpp"
#include "BCMP/Sender.hpp"
#include "DDL/L2SignalSender.hpp"
#include "System/SignalHandler.hpp"
#include "EventID.hpp"
#include "Debug/print.hpp"

using namespace dHLT;
using System::Address;


class L2SignalSender : public DDL::L2SignalSender
{
public:

	L2SignalSender(const UShort port = 4900) : DDL::L2SignalSender(port)
	{
	};
	
	virtual void OnConnect(const System::Address& address)
	{
		cout << "Connected to: " << address << endl;
		EventID eventid(2, 1);
		cout << "Signaling event: " << eventid << endl;
		SignalL2(eventid);
		Terminate();
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
	};
};


L2SignalSender* signalsender = NULL;

// If we get a SIGINT ( keyboard Ctrl+C ) signal or SIGHUP then terminate the program.
HandleSignals(SIGINT,
	if (signalsender != NULL) signalsender->Terminate();
);
HandleSignals(SIGHUP,
	if (signalsender != NULL) signalsender->Terminate();
);


void PrintUsage()
{
	cerr << "Usage: > signalL2 [-p <port>] [-n <num>] [-delay <seconds>] [-gaus | -flat | -exp] [-sigma <sigma>]" << endl;
	cerr << " -p <port> : Optional flag to specify the port number to send to." << endl;
	cerr << "  -n <num> : Optional flag specifying the number of times to send a signal." << endl;
	cerr << "             <num> must be a positive integer." << endl;
	cerr << "-delay <seconds> : Optional delay period between signals in seconds." << endl;
	cerr << "             This can be any positive floating point number. The default value is 1 sec." << endl;
	cerr << "     -gaus : Optional parameter specifying to select the delay periods from a Gaussian" << endl;
	cerr << "             distribution. The <sigma> parameter is used as the Gaussian width, and" << endl;
	cerr << "             <seconds> as the mean value." << endl;
	cerr << "     -flat : Optional parameter specifying to select the delay periods from a flat" << endl;
	cerr << "             distribution. The distribution is a box starting at <seconds> - <sigma> / 2" << endl;
	cerr << "             and ending at <seconds> + <sigma> / 2." << endl;
	cerr << "      -exp : Optional parameter specifying to select the delay periods from an" << endl;
	cerr << "             exponencial distribution. The exponent is chosen such that delays shorter" << endl;
	cerr << "             than <seconds> are equaly probable to delays larger than <seconds>" << endl;
	cerr << "-sigma <sigma> : Optional parameter used by the distribution modifiers. Default value is 1." << endl;
	cerr << "             It can be any positive floating point value." << endl;
};


enum DistributionType
{
	Delta,
	Gaussian,
	Flat,
	Exponencial
};

 
bool Odd(UInt i, const char** argv, DistributionType& distribtype)
{
	if (strcmp(argv[i],"-gaus") == 0)
	{
		distribtype = Gaussian;
		return true;
	}

	if (strcmp(argv[i],"-flat") == 0)
	{
		distribtype = Flat;
		return true;
	}

	if (strcmp(argv[i],"-exp") == 0)
	{
		distribtype = Exponencial;
		return true;
	}

	return false;
}


bool validatePort = true;

bool Even(UInt i, const char** argv, UShort& port, UInt& count, Float& delay, Float& sigma)
{
	if (strcmp(argv[i],"-p") == 0)
	{
		int x = atoi(argv[i+1]);

		if ( not (0 < x and x < 65536) )
		{                                
			validatePort = false;
			return false;
		}

		else
		{
			port = x;
			return true;
		}
	}

	if (strcmp(argv[i],"-n") == 0)
	{
		count = atoi(argv[i+1]);
		return true;
	}

	if (strcmp(argv[i],"-delay") == 0)
	{
		delay = atof(argv[i+1]);
		return true;
	}

	if (strcmp(argv[i],"-sigma") == 0)
	{
		sigma = atof(argv[i+1]);
		return true;
	}

	return false;
}


#define Even1 Even(1, argv, port, count, delay, sigma)
#define Even2 Even(2, argv, port, count, delay, sigma)
#define Even3 Even(3, argv, port, count, delay, sigma)
#define Even4 Even(4, argv, port, count, delay, sigma)
#define Even5 Even(5, argv, port, count, delay, sigma)
#define Even6 Even(6, argv, port, count, delay, sigma)
#define Even7 Even(7, argv, port, count, delay, sigma)
#define Even8 Even(8, argv, port, count, delay, sigma)

#define Odd1 Odd(1, argv, distribtype)
#define Odd3 Odd(3, argv, distribtype)
#define Odd5 Odd(5, argv, distribtype)
#define Odd7 Odd(7, argv, distribtype)
#define Odd9 Odd(9, argv, distribtype)

#define Compare1and3 strcmp(argv[1],argv[3])!=0
#define Compare1and4 strcmp(argv[1],argv[4])!=0
#define Compare1and5 strcmp(argv[1],argv[5])!=0
#define Compare1and6 strcmp(argv[1],argv[6])!=0
#define Compare1and7 strcmp(argv[1],argv[7])!=0
#define Compare1and8 strcmp(argv[1],argv[8])!=0
#define Compare2and4 strcmp(argv[2],argv[4])!=0
#define Compare2and6 strcmp(argv[2],argv[6])!=0
#define Compare2and8 strcmp(argv[2],argv[8])!=0
#define Compare3and5 strcmp(argv[3],argv[5])!=0
#define Compare3and6 strcmp(argv[3],argv[6])!=0
#define Compare3and7 strcmp(argv[3],argv[7])!=0
#define Compare3and8 strcmp(argv[3],argv[8])!=0
#define Compare4and6 strcmp(argv[4],argv[6])!=0
#define Compare4and8 strcmp(argv[4],argv[8])!=0
#define Compare5and7 strcmp(argv[5],argv[7])!=0
#define Compare5and8 strcmp(argv[5],argv[8])!=0
#define Compare6and8 strcmp(argv[6],argv[8])!=0


/* Parse command line and fill the parameters.
   Returns true if the command line is parsed properly.
 */
bool ParseCommandLine(const int argc, const char** argv,
		UShort& port, UInt& count, DistributionType& distribtype, Float& delay, Float& sigma)
{
	for (int i=1; i<argc; i++)
	{
		if (strcmp(argv[i], "-h") == 0 or strcmp(argv[i], "-help") == 0
			or strcmp(argv[i], "--help") == 0)
		{
			PrintUsage();
			return false;
		}
	}

	switch(argc)
	{
		case 2:
			if (Odd1)
				return true;
			break;

		case 3:
			if (Even1)
				return true;
			break;

		case 4:
			if ((Odd1 and Even2) or (Even1 and Odd3))
				return true;
			break;

		case 5:
			if (Even1 and Even3 and Compare1and3)
				return true;
			break;

		case 6:
			if (	(Odd1 and Even2 and Even4 and Compare2and4)
				or (Even1 and Odd3 and Even4 and Compare1and4)
				or (Even1 and Even3 and Odd5 and Compare1and3)
			   )
				return true;
			break;

		case 7:
			if (Even1 and Even3 and Even5 and Compare1and3 and Compare1and5 and Compare3and5)
				return true;
			break;

		case 8:
			if (	(Odd1 and Even2 and Even4 and Even6 and Compare2and4 and Compare2and6 and Compare4and6)
				or (Even1 and Odd3 and Even4 and Even6 and Compare1and4 and Compare1and6 and Compare4and6)
				or (Even1 and Even3 and Odd5 and Even6 and Compare1and3 and Compare1and6 and Compare3and6)
				or (Even1 and Even3 and Even5 and Odd7 and Compare1and3 and Compare1and5 and Compare3and5)
			   )
				return true;
			break;

		case 9:
			if (	Even1 and Even3 and Even5 and Even7 and Compare1and3 and Compare1and5
				and Compare1and7 and Compare3and5 and Compare3and7 and Compare5and7
			   )
				return true;
			break;

		case 10:
			if (	(Odd1 and Even2 and Even4 and Even6 and Even8 and Compare2and4 and
				 Compare2and6 and Compare2and8 and Compare4and6 and Compare4and8 and Compare6and8)
				or (Even1 and Odd3 and Even4 and Even6 and Even8 and Compare1and4 and
				    Compare1and6 and Compare1and8 and Compare4and6 and Compare4and8 and Compare6and8)
				or (Even1 and Even3 and Odd5 and Even6 and Even8 and Compare1and3 and
				    Compare1and6 and Compare1and8 and Compare3and6 and Compare3and8 and Compare6and8)
				or (Even1 and Even3 and Even5 and Odd7 and Even8 and Compare1and3 and
				    Compare1and5 and Compare1and8 and Compare3and5 and Compare3and8 and Compare5and8)
				or (Even1 and Even3 and Even5 and Even7 and Odd9 and Compare1and3 and
				    Compare1and5 and Compare1and7 and Compare3and5 and Compare3and7 and Compare5and7))
				return true;
	}

	if (validatePort)
		PrintUsage();
	else
		cerr << "Error: Port is out of range, valid range is [1..65535]." << endl;

	return false;
}



int main(const int argc, const char** argv)
{
	int returncode = 0;
	try
	{
		UShort port = 4900;
		UInt sigcount = 1;
		DistributionType distribtype = Delta;
		Float delay = 1.0;  // 1 second
		Float sigma = 1.0;
		
		bool valid_arguments = ParseCommandLine(argc, argv, port, sigcount, distribtype, delay, sigma);
		if (valid_arguments)
		{
			L2SignalSender sender(port);
			signalsender = &sender;
			cout << "Waiting for connections..." << endl;
			sender.Run();
		}
		else
			returncode = 1;
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
