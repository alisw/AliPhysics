/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$

/**
 * \file   hltdimserver.cxx
 * \author Artur Szostak <artursz@iafrica.com>
 * \date   16 Feb 2011
 * \brief  Standalone DIM server to publish HLT information like the reconstructed vertex.
 */

#include "TSystem.h"
#include "TString.h"
#include "TH1.h"
#include "AliHLTDimServer.h"
#include "AliHLTHOMERReader.h"
#include "AliHLTLogging.h"
#include "AliHLTMessage.h"
#include "AliHLTOUT.h"

#include <cstring>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <iostream>
using std::cerr;
using std::endl;


#define DEFAULT_DIM_DNS             "alidcsdimdns.cern.ch" 
#define DEFAULT_SERVER_NAME         "ALICE_HLT"
#define DEFAULT_VERTEX_SERVICE_NAME "VERTEX"


#define CMDLINE_ERROR 1
#define SYSTEM_ERROR 2
#define FATAL_ERROR 3


/**
 * The signal handler class is used to trap quit and interupt signals from
 * the system so as to give the server a chance to shutdown cleanly.
 */
class AliHLTSignalHandler : public TSignalHandler, public AliHLTLogging
{
public:
	
	AliHLTSignalHandler(ESignals sig) : TSignalHandler(sig, kFALSE), AliHLTLogging() {}
	
	virtual Bool_t Notify()
	{
		fgTerminationSignaled = true;
		return TSignalHandler::Notify();
	}
	
	static bool TerminationSignaled() { return fgTerminationSignaled; }
	
private:
	
	static bool fgTerminationSignaled;
};


bool AliHLTSignalHandler::fgTerminationSignaled = false;


/**
 * The AliHLTDCSPublisherServer class implements the HLT specific DIM
 * server for publishing HLT information. In particular the HLT vertex
 * information.
 */
class AliHLTDCSPublisherServer : public AliHLTLogging
{
public:
	AliHLTDCSPublisherServer();
	~AliHLTDCSPublisherServer();
	
	/**
	 * Parses the command line arguments given and initialises relevant
	 * internal variables.
	 * \param argc  Number of arguments as given in main().
	 * \param argv  Array of arguments as given in main().
	 * \return  A status flag suitable for returning from main(), containing
	 *          one of EXIT_SUCCESS, CMDLINE_ERROR.
	 */
	int Init(int argc, char** argv);
	
	/**
	 * Implements the main server loop. Connects to the HOMER port specified
	 * and pushes the data from the chain via DIM to DCS.
	 * \return  A status flag suitable for returning from main(), containing
	 *      one of EXIT_SUCCESS, FATAL_ERROR or SYSTEM_ERROR.
	 */
	int Run();
	
private:
	
	// Do not allow copying of this class.
	AliHLTDCSPublisherServer(const AliHLTDCSPublisherServer&);
	AliHLTDCSPublisherServer& operator = (const AliHLTDCSPublisherServer&);
	
	
	struct AliLuminosityRegion
	{
		float fX;  /// X coordinate of beam spot centroid [micrometers].
		float fY;  /// Y coordinate of beam spot centroid [micrometers].
		float fZ;  /// Z coordinate of beam spot centroid [mm].
		float fSizeX;  /// Beam spot size in X direction (RMS) [micrometers].
		float fSizeY;  /// Beam spot size in Y direction (RMS) [micrometers].
		float fSizeZ;  /// Beam spot size in Z direction (RMS) [mm].
		float fDxdz;  /// Beam spot tilt angle dx/dz [microradians].
		float fDydz;  /// Beam spot tilt angle dy/dz [microradians].
	};
	
	/**
	 * Tries to establish the HOMER connection if not already made.
	 */
	void MakeHOMERConnection();
	
	/**
	 * Fills the luminosity region structure.
	 */
	void SetLuminosityRegion(AliLuminosityRegion& lumiRegion);
	
	/**
	 * Prints the command line usage of this program to the HLT logging facility.
	 * \param asError  If true then the output is generated as an error message
	 *      and as an informational otherwise.
	 */
	void PrintUsage(bool asError = true);
	
	TString fHomerSource;  /// Hostname of the server running the TCPDumpSubscriber for the HOMER connection.
	unsigned short fHomerPort;   /// The TCP port number on which to connect for the HOMER connection.
	TString fDimDns;  /// The name of the DIM DNS to register with.
	TString fServerName;  /// The server name to use for this DIM server.
	TString fVertexServiceName;  /// The name of the vertex service to register.
	AliHLTHOMERReader* fHomerReader;  /// HOMER reader to access HLT data.
	UInt_t fPublishPeriod;  /// The period between publishing data points to DIM in milliseconds.
	
	virtual const char* Class_Name() { return "AliHLTDCSPublisherServer"; }
};


AliHLTDCSPublisherServer::AliHLTDCSPublisherServer() :
	AliHLTLogging(),
	fHomerSource(),
	fHomerPort(0),
	fDimDns(DEFAULT_DIM_DNS),
	fServerName(DEFAULT_SERVER_NAME),
	fVertexServiceName(DEFAULT_VERTEX_SERVICE_NAME),
	fHomerReader(NULL),
	fPublishPeriod(1000)
{
	// Default constructor.
}


AliHLTDCSPublisherServer::~AliHLTDCSPublisherServer()
{
	// Default destructor.
}


int AliHLTDCSPublisherServer::Init(int argc, char** argv)
{
	// Parses the command line and initialises the internal variables.
	
	const char* homerSource = NULL;
	unsigned short homerPort = 0;
	const char* dimdns = NULL;
	const char* serverName = NULL;
	const char* vertexServiceName = NULL;
	const char* publishPeriod = NULL;

	// Parse the command line.
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-help") == 0 or strcmp(argv[i], "-h") == 0)
		{
			PrintUsage(false);
			return EXIT_SUCCESS;
		}
		else if (strcmp(argv[i], "-source") == 0 or strcmp(argv[i], "-s") == 0)
		{
			if (homerSource != NULL)
			{
				HLTWarning("Already used -source|-s with \"%s\" and port = %d before."
					   " Will override it with the last value specified with -source|-s.",
					   homerSource, int(homerPort)
				);
			}
			if (++i >= argc)
			{
				HLTError("Missing the hostname for -source|-s option.");
				PrintUsage();
				return CMDLINE_ERROR;
			}
			homerSource = argv[i];
			if (++i >= argc)
			{
				HLTError("Missing the TCP port number for -source|-s option.");
				PrintUsage();
				return CMDLINE_ERROR;
			}
			char* cpErr = NULL;
			long num = strtol(argv[i], &cpErr, 0);
			if (cpErr == NULL or *cpErr != '\0' or num < 0 or num > 65535)
			{
				HLTError("Cannot convert \"%s\" to a valid port integer number in the range [0..65535].", argv[i]);
				return CMDLINE_ERROR;
			}
			homerPort = (unsigned short)num;
		}
		else if (strcmp(argv[i], "-dimdns") == 0)
		{
			if (dimdns != NULL)
			{
				HLTWarning("Already used -dimdns with the value \"%s\" before."
					   " Will override it with the last value specified with -dimdns.",
					   dimdns
				);
			}
			if (++i >= argc)
			{
				HLTError("Missing the DIM DNS hostname to register with.");
				PrintUsage();
				return CMDLINE_ERROR;
			}
			dimdns = argv[i];
		}
		else if (strcmp(argv[i], "-servername") == 0)
		{
			if (serverName != NULL)
			{
				HLTWarning("Already used -servername with the value \"%s\" before."
					   " Will override it with the last value specified with -servername.",
					   serverName
				);
			}
			if (++i >= argc)
			{
				HLTError("Missing the name of the server to use for this DIM server.");
				PrintUsage();
				return CMDLINE_ERROR;
			}
			serverName = argv[i];
		}
		else if (strcmp(argv[i], "-vertexservicename") == 0)
		{
			if (vertexServiceName != NULL)
			{
				HLTWarning("Already used -vertexservicename with the value \"%s\" before."
					   " Will override it with the last value specified with -vertexservicename.",
					   vertexServiceName
				);
			}
			if (++i >= argc)
			{
				HLTError("Missing the name of the service to publish with DIM.");
				PrintUsage();
				return CMDLINE_ERROR;
			}
			vertexServiceName = argv[i];
		}
		else if (strcmp(argv[i], "-publishperiod") == 0)
		{
			if (publishPeriod != NULL)
			{
				HLTWarning("Already used -publishperiod with the value \"%s\" before."
					   " Will override it with the last value specified with -publishperiod.",
					   publishPeriod
				);
			}
			if (++i >= argc)
			{
				HLTError("Missing the publishing period value in milliseconds.");
				PrintUsage();
				return CMDLINE_ERROR;
			}
			publishPeriod = argv[i];
			char* cpErr = NULL;
			unsigned long num = strtoul(argv[i], &cpErr, 0);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert \"%s\" to a valid unsigned integer.", argv[i]);
				return CMDLINE_ERROR;
			}
			fPublishPeriod = UInt_t(num);
		}
		else
		{
			HLTError("Unknown command line option \"%s\".", argv[i]);
			PrintUsage();
			return CMDLINE_ERROR;
		}
	}
	
	if (homerSource == NULL)
	{
		HLTError("No HOMER source hostname or port was specified with the -source|-s option.");
		PrintUsage();
		return CMDLINE_ERROR;
	}
	
	fHomerSource = homerSource;
	fHomerPort = homerPort;
	if (dimdns != NULL) fDimDns = dimdns;
	if (serverName != NULL) fServerName = serverName;
	if (vertexServiceName != NULL) fVertexServiceName = vertexServiceName;
	return EXIT_SUCCESS;
}


int AliHLTDCSPublisherServer::Run()
{
	// Runs the main loop of the DCS publisher server.
	
	HLTImportant("Starting HLT DCS publishing server");
	
	// Create and register the DIM server instance.
	AliHLTDimServer dimServer(fServerName.Data());
	if (dimServer.Init(fDimDns.Data()) != 0)
	{
		HLTError("Could not register with DIM DNS \"%s\"", fDimDns.Data());
		return SYSTEM_ERROR;
	}
	HLTInfo("Registered DIM server with DNS \"%s\"", fDimDns.Data());
	
	// Create and register the service.
	AliLuminosityRegion lumiRegion;
	memset(&lumiRegion, 0x0, sizeof(lumiRegion));
	AliHLTDimServer::AliHLTDimService* service = new AliHLTDimServer::AliHLTDimService("F:8", &lumiRegion, sizeof(lumiRegion), fVertexServiceName.Data());
	if (dimServer.RegisterService(service) != 0)
	{
		HLTError("Could not create new DIM service \"%s\"", fVertexServiceName.Data());
		return SYSTEM_ERROR;
	}
	HLTInfo("Created \"%s\" DIM service.", fVertexServiceName.Data());
	
	// Finally start the DIM server loop.
	if (dimServer.Start() != 0)
	{
		HLTError("Could not start the DIM server \"%s\".", fServerName.Data());
		return SYSTEM_ERROR;
	}
	HLTInfo("Started DIM server \"%s\".", fServerName.Data());
	
	while (not AliHLTSignalHandler::TerminationSignaled())
	{
		MakeHOMERConnection();
		if (fHomerReader != NULL)
		{
			int result = fHomerReader->ReadNextEvent(10000000); // timeout of 10 seconds.
			if (result != ETIMEDOUT and result != 0)
			{
				// Do not show as error message since it is expected to loose
				// HOMER connections as runs stop and start.
				HLTDebug("Error occurred trying to read from HOMER connection to %s:%d",
					fHomerSource.Data(), int(fHomerPort)
				);
				continue;
			}
			if (result != ETIMEDOUT)
			{
				SetLuminosityRegion(lumiRegion);
				HLTDebug("Update luminosity region: Beam spot (x [um], y [um], z [mm]) = (%f, %f, %f);"
					" size (x [um], y [um], z [mm]) = (%f, %f, %f); tilt (dx/dz [urad], dy/dz [urad]) = (%f, %f).",
					lumiRegion.fX, lumiRegion.fY, lumiRegion.fZ,
					lumiRegion.fSizeX, lumiRegion.fSizeY, lumiRegion.fSizeZ,
					lumiRegion.fDxdz, lumiRegion.fDydz
				);
				service->Update();
				gSystem->Sleep(fPublishPeriod);
			}
			else
			{
				memset(&lumiRegion, 0x0, sizeof(lumiRegion));
				service->Update();
			}
		}
		else
		{
			memset(&lumiRegion, 0x0, sizeof(lumiRegion));
			service->Update();
			gSystem->Sleep(10000);  // sleep for 10 seconds.
		}
	}
	
	HLTImportant("Stopping HLT DCS publishing server");
	
	// Stop the DIM server loop.
	if (dimServer.Stop() != 0)
	{
		HLTError("Could not stop the DIM server properly.");
		return SYSTEM_ERROR;
	}
	HLTInfo("Stopped DIM server \"%s\".", fServerName.Data());
	
	return EXIT_SUCCESS;
}


void AliHLTDCSPublisherServer::MakeHOMERConnection()
{
	// Makes a connection to the HLT data via HOMER.
	
	if (fHomerReader != NULL)
	{
		// Delete the old connection object if there was an error and try again.
		if (fHomerReader->GetConnectionStatus() != 0 and fHomerReader->GetConnectionStatus() != ETIMEDOUT)
		{
			HLTDebug("Previous connection to HOMER source %s:%d has an error (code = %d).",
				 fHomerSource.Data(), int(fHomerPort), fHomerReader->GetConnectionStatus()
			);
			delete fHomerReader;
			fHomerReader = NULL;
		}
	}
	
	if (fHomerReader == NULL)
	{
		fHomerReader = new AliHLTHOMERReader(fHomerSource.Data(), fHomerPort);
		if (fHomerReader->GetConnectionStatus() != 0)
		{
			HLTDebug("Could not establish connection to HOMER source %s:%d, error code = %d.",
				 fHomerSource.Data(), int(fHomerPort), fHomerReader->GetConnectionStatus()
			);
			delete fHomerReader;
			fHomerReader = NULL;
		}
		else
		{
			HLTInfo("Established new HOMER connection to %s:%d.", fHomerSource.Data(), int(fHomerPort));
		}
	}
}


void AliHLTDCSPublisherServer::SetLuminosityRegion(AliLuminosityRegion& lumiRegion)
{
	// Sets the values for the luminosity region.
	
	assert(fHomerReader != NULL);
	
	const char* histName[6] = {"primVertexX", "primVertexY", "primVertexZ", "spdVertexX", "spdVertexY", "spdVertexZ"};
	TH1* hist[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
	
	// First find the histograms in the list of data blocks read from HOMER.
	for (unsigned long i = 0; i < fHomerReader->GetBlockCnt(); ++i)
	{
		homer_uint64 type = fHomerReader->GetBlockDataType(i);
		homer_uint32 origin = fHomerReader->GetBlockDataOrigin(i);
		AliHLTComponentDataType dt;
		AliHLTComponent::SetDataType(dt, AliHLTOUT::ByteSwap64(type), AliHLTOUT::ByteSwap32(origin));
		HLTDebug("Received data block of type %s, specification = 0x%X",
			 AliHLTComponent::DataType2Text(dt).c_str(), fHomerReader->GetBlockDataSpec(i)
		);
		
		if (dt != (kAliHLTDataTypeHistogram | kAliHLTDataOriginAny)) continue;
		
		TObject* obj = AliHLTMessage::Extract(fHomerReader->GetBlockData(i), fHomerReader->GetBlockDataLength(i));
		if (obj == NULL)
		{
			HLTWarning("Problem extracting object from data for event %ul (0x%Xl), from data block %u, type = %s",
				   fHomerReader->GetEventID(), fHomerReader->GetEventID(), i,
				   AliHLTComponent::DataType2Text(dt).c_str()
			);
			continue;
		}
		TH1* h = dynamic_cast<TH1*>(obj);
		if (h == NULL)
		{
			delete obj;
			continue;
		}
		HLTDebug("Found histogram object \"%s\".", h->GetName());
		bool assigned = false;
		for (int j = 0; j < 6; ++j)
		{
			if (TString(h->GetName()) == histName[j])
			{
				HLTDebug("Marking histogram \"%s\" as found in hist[%d].", h->GetName(), j);
				// Note: This object will be marked for now and deleted later
				// at the end of this method.
				hist[j] = h;
				assigned = true;
				break;
			}
		}
		if (not assigned) delete obj;
	}
	
	// Now get the Mean and RMS values as the centroid and beam size values,
	// Scale the values from centimetres to the required units of microns and mm.
	memset(&lumiRegion, 0x0, sizeof(lumiRegion));
	if (hist[3] != NULL)
	{
		lumiRegion.fX = hist[3]->GetMean() * 1e4;
		lumiRegion.fSizeX = hist[3]->GetRMS() * 1e4;
	}
	if (hist[4] != NULL)
	{
		lumiRegion.fY = hist[4]->GetMean() * 1e4;
		lumiRegion.fSizeY = hist[4]->GetRMS() * 1e4;
	}
	if (hist[5] != NULL)
	{
		lumiRegion.fZ = hist[5]->GetMean() * 10;
		lumiRegion.fSizeZ = hist[5]->GetRMS() * 10;
	}
	
	//TODO: for now just fill zeros into the tilt angles since they are not yet calculated.
	lumiRegion.fDxdz = 0;
	lumiRegion.fDydz = 0;
	
	for (int j = 0; j < 6; ++j)
	{
		if (hist[j] != NULL) delete hist[j];
	}
}


void AliHLTDCSPublisherServer::PrintUsage(bool asError)
{
	// Prints the command line usage.
	
	std::stringstream os;
	os << "Usage: hltdimserver -source|-s <hostname> <port> [-help|-h] [-dimdns <DIM-DNS-name>]" << endl;
	os << "         [-servername <server-name>] [-vertexservicename <service-name>]" << endl;
	os << "Options:" << endl;
	os << " -source | -s <hostname> <port>" << endl;
	os << "       Specifies the hostname and port to connect to via HOMER to extract data from" << endl;
	os << "       the online chain from a TCPDumpSubscriber." << endl;
	os << " -help | -h" << endl;
	os << "       Displays this message." << endl;
	os << " -dimdns <DIM-DNS-name>" << endl;
	os << "       The name of the DIM DNS to register the this DIM server with." << endl;
	os << " -servername <server-name>" << endl;
	os << "       The name to register with for this DIM server." << endl;
	os << " -vertexservicename <service-name>" << endl;
	os << "       The name to register with of the vertex information service." << endl;
	os << " -publishperiod <time>" << endl;
	os << "       The amount of time in milliseconds to wait between publishing data points to" << endl;
	os << "       DIM. This controls the publishing rate." << endl;
	
	if (asError)
	{
		HLTError(os.str().c_str());
	}
	else
	{
		HLTInfo(os.str().c_str());
	}
}


int main(int argc, char** argv)
{
	int returnCode = EXIT_SUCCESS;
	try
	{
		gSystem->AddSignalHandler(new AliHLTSignalHandler(kSigPipe));
		gSystem->AddSignalHandler(new AliHLTSignalHandler(kSigQuit));
		gSystem->AddSignalHandler(new AliHLTSignalHandler(kSigInterrupt));
		gSystem->AddSignalHandler(new AliHLTSignalHandler(kSigTermination));
		AliHLTDCSPublisherServer server;
		returnCode = server.Init(argc, argv);
		if (returnCode == EXIT_SUCCESS) returnCode = server.Run();
	}
	catch (...)
	{
		cerr << "FATAL ERROR: An unknown exception occurred!" << endl;
	}
	return returnCode;
}
