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
#include "AliHLTLogging.h"
#include "AliHLTMessage.h"
#include "AliHLTOUT.h"
#include "AliZMQhelpers.h"
#include "zmq.h"
#include "AliHLTLumiRegComponent.h"

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

using namespace AliZMQhelpers;


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
	 * Implements the main server loop. Connects to the ZMQ port specified
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
	 * Initialises the luminosity buffer with values indicating an invalid data point.
	 */
	void InitLumiRegionAsInvalid(AliLuminosityRegion& lumiRegion);
	
	/**
	 * Fills the luminosity region structure.
	 */
	int SetLuminosityRegion(AliLuminosityRegion& lumiRegion);
	
	/**
	 * Prints the command line usage of this program to the HLT logging facility.
	 * \param asError  If true then the output is generated as an error message
	 *      and as an informational otherwise.
	 */
	void PrintUsage(bool asError = true);
	
  void* fZMQcontext;
  void* fZMQin;
  int fZMQtimeout;
	TString fZMQconfigIN;  /// ZMQ source coordinates
	TString fDimDns;  /// The name of the DIM DNS to register with.
	TString fServerName;  /// The server name to use for this DIM server.
	TString fVertexServiceName;  /// The name of the vertex service to register.
	UInt_t fPublishPeriod;  /// The period between publishing data points to DIM in milliseconds.
  
  //hold histograms
  TH1F* fPrimary[3]; 
  TH1F* fPrimaryDefMult[3]; 
	
	virtual const char* Class_Name() { return "AliHLTDCSPublisherServer"; }
};


AliHLTDCSPublisherServer::AliHLTDCSPublisherServer() :
	AliHLTLogging(),
  fZMQcontext(NULL),
  fZMQin(NULL),
  fZMQtimeout(1000),
	fZMQconfigIN(),
	fDimDns(DEFAULT_DIM_DNS),
	fServerName(DEFAULT_SERVER_NAME),
	fVertexServiceName(DEFAULT_VERTEX_SERVICE_NAME),
	fPublishPeriod(1000),
  fPrimary(),
  fPrimaryDefMult()
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
	
	const char* zmqSource = NULL;
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
			if (zmqSource != NULL)
			{
				HLTWarning("Already used -source|-s with \"%s\""
					   " Will override it with the last value specified with -source|-s.",
					   zmqSource
				);
			}
			if (++i >= argc)
			{
				HLTError("no value specified for -source|-s option.");
				PrintUsage();
				return CMDLINE_ERROR;
			}
			zmqSource = argv[i];
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
	
	if (zmqSource == NULL)
	{
		HLTError("No ZMQ source was specified with the -source|-s option.");
		PrintUsage();
		return CMDLINE_ERROR;
	}
	
	fZMQconfigIN = zmqSource;
	if (dimdns != NULL) fDimDns = dimdns;
	if (serverName != NULL) fServerName = serverName;
	if (vertexServiceName != NULL) fVertexServiceName = vertexServiceName;

  fZMQcontext = zmq_ctx_new();
  if (alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data()) < 0)
    return -1;
  
  for(int i=0;i<3;i++){fPrimary[i]=NULL;}
  for(int i=0;i<3;i++){fPrimaryDefMult[i]=NULL;}

	return EXIT_SUCCESS;
}


int AliHLTDCSPublisherServer::Run()
{
	// Runs the main loop of the DCS publisher server.
	
	HLTImportant("Starting HLT DCS publishing server");
	
	AliLuminosityRegion lumiRegion;
	InitLumiRegionAsInvalid(lumiRegion);

	//Create and register the DIM server instance.
	AliHLTDimServer dimServer(fServerName.Data());
	if (dimServer.Init(fDimDns.Data()) != 0)
	{
		HLTError("Could not register with DIM DNS \"%s\"", fDimDns.Data());
		return SYSTEM_ERROR;
	}
	HLTInfo("Registered DIM server with DNS \"%s\"", fDimDns.Data());
	
	// Create and register the service.
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
    if (SetLuminosityRegion(lumiRegion)==0){
      HLTInfo("Update luminosity region: Beam spot (x [um], y [um], z [mm]) = (%f, %f, %f);"
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
			InitLumiRegionAsInvalid(lumiRegion);
      HLTInfo("Update luminosity region: Beam spot (x [um], y [um], z [mm]) = (%f, %f, %f);"
          " size (x [um], y [um], z [mm]) = (%f, %f, %f); tilt (dx/dz [urad], dy/dz [urad]) = (%f, %f).",
          lumiRegion.fX, lumiRegion.fY, lumiRegion.fZ,
          lumiRegion.fSizeX, lumiRegion.fSizeY, lumiRegion.fSizeZ,
          lumiRegion.fDxdz, lumiRegion.fDydz
          );
			service->Update();
      gSystem->Sleep(fPublishPeriod);
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

void AliHLTDCSPublisherServer::InitLumiRegionAsInvalid(AliLuminosityRegion& lumiRegion)
{
	// Initialises the luminosity regions with values indicating an invalid data point.
	
  //reset the histograms
  for (int i=0;i<3;i++){
      delete fPrimary[i]; fPrimary[i]=NULL;
  }
  for (int i=0;i<3;i++){
      delete fPrimaryDefMult[i]; fPrimaryDefMult[i]=NULL;
  }

	memset(&lumiRegion, 0x0, sizeof(lumiRegion));
	lumiRegion.fX = -1000.;
	lumiRegion.fY = -1000.;
	lumiRegion.fZ = -1000.;
	lumiRegion.fSizeX = -1000.;
	lumiRegion.fSizeY = -1000.;
	lumiRegion.fSizeZ = -1000.;
	lumiRegion.fDxdz = 1000.;
	lumiRegion.fDydz = -1000.;
}


int AliHLTDCSPublisherServer::SetLuminosityRegion(AliLuminosityRegion& lumiRegion)
{
	// Sets the values for the luminosity region.
	int rc = 0;
  std::string primaryName[3] = {"hXTRKVtx", "hYTRKVtx", "hZTRKVtx" };
  std::string primaryDefMultName[3] = {"hXTRKDefMult", "hYTRKDefMult", "hZTRKDefMult"};

  //invalidate previous values
	InitLumiRegionAsInvalid(lumiRegion);


  //send a request if we are using REQ
  int sockettype = alizmq_socket_type(fZMQin);
  if (sockettype == ZMQ_REQ)
  {
    //here we send an empty request which will get us our data
    //and a command to reset the merger after
    HLTInfo("sending request");
    alizmq_msg_send("*","",fZMQin,ZMQ_SNDMORE);
    alizmq_msg_send("CONFIG","reset",fZMQin,0);
  }

  //wait for the data
  zmq_pollitem_t sockets[] = { { fZMQin, 0, ZMQ_POLLIN, 0 } };
  rc = zmq_poll(sockets, 1, (sockettype == ZMQ_REQ) ? fZMQtimeout : -1);
  if (rc<0)
  {
    //eternal interrupt - exit
    return 0;
  }
  if (!(sockets[0].revents & ZMQ_POLLIN))
  {
    //server died, reinit socket
    HLTInfo("server not responding...");
    alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data());
    return 0;
  }

  aliZMQmsg message;
  alizmq_msg_recv(&message, fZMQin, 0);

  //sort incoming histograms into proper slots
  int numberOfUsefulHistograms=0;
  for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
  {
    //Printf("have some data!");
    TObject* object;
    alizmq_msg_iter_data(i, object);
    if (!object) {
      //Printf("no tobject!");
      continue;
    }
    //Printf("have an TObject");
    TH1F* hist = dynamic_cast<TH1F*>(object);
    if (!hist) {delete object; continue;}
    std::string name = hist->GetName();
    //Printf("is a TH1F named: %s",name);
    bool useful=false;
    for (int i=0;i<3;i++){
      if (name.find(primaryName[i])!=std::string::npos) {
        delete fPrimary[i]; fPrimary[i]=NULL;
        //Printf("matched name: %s, i: %i",name,i);
        fPrimary[i]=hist;
        useful=true;
        numberOfUsefulHistograms++;
      }
    }
    for (int i=0;i<3;i++){
      if (name.find(primaryDefMultName[i])!=std::string::npos){
        delete fPrimaryDefMult[i]; fPrimaryDefMult[i]=NULL;
        //Printf("matched name: %s, i: %i",name,i);
        fPrimaryDefMult[i]=hist;
        useful=true;
        numberOfUsefulHistograms++;
      }
    }
    if (!useful) delete object;
  }
  
  alizmq_msg_close(&message);

  if (numberOfUsefulHistograms!=6) {
    Printf("received %i histograms!, should be 6",numberOfUsefulHistograms);
    return 1;
  }

  Printf("%p %p %p %p %p %p", fPrimary[0],fPrimary[1],fPrimary[2],fPrimaryDefMult[0],fPrimaryDefMult[1],fPrimaryDefMult[2]);
  if (!(fPrimary[0] && fPrimary[1] && fPrimary[2] && fPrimaryDefMult[0] && fPrimaryDefMult[1] && fPrimaryDefMult[2]))
    return 0;
  Printf("starting fit");

  //do the fits
  Float_t meanVtx[3] ={0.,0.,0.};
  Float_t sigmaVtx[3] ={0.,0.,0.};
  /*Int_t fitVtxResults =*/ AliHLTLumiRegComponent::FitPositions(fPrimary, meanVtx, sigmaVtx);
  Float_t meanLR[3] ={0.,0.,0.};
  Float_t sigmaLR[3] ={0.,0.,0.};
  Int_t lumiRegResults = AliHLTLumiRegComponent::LuminousRegionExtraction(fPrimaryDefMult, meanLR, sigmaLR);

  lumiRegion.fX = meanVtx[0];
  lumiRegion.fSizeX = (lumiRegResults==1 || lumiRegResults==3) ? sigmaVtx[0] : sigmaLR[0];
  lumiRegion.fY = meanVtx[1];
  lumiRegion.fSizeY = (lumiRegResults==2 || lumiRegResults==3) ? sigmaVtx[1] : sigmaLR[1];
  lumiRegion.fZ = meanVtx[2];
  lumiRegion.fSizeZ = sigmaVtx[2];
	//TODO: for now just fill -999 into the tilt angles since they are not yet calculated.
	lumiRegion.fDxdz = -999.;
	lumiRegion.fDydz = -999.;
  return 0;
}

void AliHLTDCSPublisherServer::PrintUsage(bool asError)
{
	// Prints the command line usage.
	
	std::stringstream os;
	os << "Usage: hltdimserver -source|-s <ZMQconfig> [-help|-h] [-dimdns <DIM-DNS-name>]" << endl;
	os << "         [-servername <server-name>] [-vertexservicename <service-name>]" << endl;
	os << "Options:" << endl;
	os << " -source | -s <ZMQmode><@|-><transport>://<endpoint>" << endl;
	os << "       e.g. PULL@tcp://*:60201,>tcp://192.168.56.1:3456" << endl;
	os << "            REQ@inproc://someSocket,>ipc:///tmp/somesocket" << endl;
	os << "            SUB>ipc:///tmp/somesocket" << endl;
	os << "            > connects, @ binds, multiple endpoints are allowed" << endl;
	os << "            modes are any of: PULL,REQ,SUB (other ones dont make much sense here)" << endl;
	os << "            see zeromq docs for details << endl;" << endl;
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
