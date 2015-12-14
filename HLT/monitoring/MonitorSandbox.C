/**************************************************************************
 * Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
 * \file MonitorSandbox.C
 * \date 8 Oct 2011
 * \author Artur Szostak <artursz@iafrica.com>
 * \brief Macro providing a very rudimentary implementation for running a monitoring chain asynchronously to the HLT main chain.
 * 
 * This macro provides a small light weight implementation of a monitoring chain.
 * It connects the the HLT system using the HOMER interface. This allows extracting
 * online data in an asynchronous manner without affecting the main chain in any
 * way. The extracted data is then passed through a small configuration defined
 * inside the Configure() routine. The output of which is then published again via
 * a very light weight implementation of the HOMER interface. The results then
 * can be connected to by DQM using the HOMER interface to display the histogram
 * results.
 * 
 * The only modifications needed to be done by end users should be in Configure(),
 * where the HLT mini monitoring chain is defined. There are two default components
 * already created, the "source" and "sink" components. These are the end points
 * for the monitoring chain. The "source" components simply publishes the data
 * blocks as received via the HOMER interface into the chain. While the "sink"
 * component extracts any output data blocks it gets, makes a copy in memory,
 * which are ultimately published to any clients connected to this process.
 * 
 * An example of running the macro is given below:
 * \code
 *   $ aliroot -b -q -l MonitorSandbox.C'"localhost",49000,49001,0.1,"local://$ALICE_ROOT/OCDB")'
 * \endcode
 * where the parameters are, the host name to which to connect with HOMER to fetch
 * the data, the port number to which to connect HOMER to, the server side port
 * on which other clients should connect to the monitoring chain i.e. DQM, the
 * HOMER polling period in seconds and finally the CDB path to use.
 * For production it should be "local:///opt/HCDB"
 */

#include "AliHLTHOMERReader.h"
#include "AliHLTHOMERWriter.h"
#include "AliHLTOUT.h"
#include "AliHLTMisc.h"
#include "AliHLTSystem.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTPluginBase.h"
#include "AliHLTConfiguration.h"
#include "AliHLTHOMERLibManager.h"
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTDataSource.h"
#include "AliHLTProcessor.h"
#include "AliHLTDataSink.h"
#include "AliHLTMessage.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTCTPData.h"
#include "AliCDBManager.h"
#include "AliDAQ.h"
#include "TServerSocket.h"
#include "TObjArray.h"
#include "TSocket.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

#ifndef __CINT__
#include <sys/ioctl.h>
#include <arpa/inet.h>
#endif


enum
{
	kGlobalHistoChainConfig = 1,
	kCompressionHistoChainConfig = 2,
	kCorrelationMinBiasHistoChainConfig = 3,
	kCorrelationSemiCentralHistoChainConfig = 4,
	kCorrelationCentralHistoChainConfig = 5,
	kCorrelationAllTriggerHistoChainConfig = 6,
	kEmcalHistoChainConfig = 7,
	kMuonHistoChainConfig = 8,
	kCorrelationsHistoChainConfig = 9
};


/**
 * The configure macro can be edited by the end used to add more components to the
 * monitoring chain.
 * The 3 things to edit are:
 * 1) Add the additional library where the component is defined if not already done.
 * 2) Add the components configuration with AliHLTConfiguration. The parent of the
 *    component should ultimately be called "source", which is the MonitorHomerPublisher
 *    which injects the data coming from the main HLT chain via HOMER.
 * 3) Add the name of the new component (first parameter used in AliHLTConfiguration)
 *    to the listOfChains string so that the monitor component's output will be
 *    collected and published to DQM. Remember its a space separated string of names.
 * \param configId The exact configuration to use for this instance.
 */
void Configure(int configId)
{
	AliHLTSystem* system = AliHLTPluginBase::GetInstance();
	TString listOfChains = "";
	
	// Add additional libraries that need to be loaded here:
	system->LoadComponentLibraries("libAliHLTUtil.so");
	system->LoadComponentLibraries("libAliHLTGlobal.so");
	system->LoadComponentLibraries("libAliHLTTPC.so");
	system->LoadComponentLibraries("libAliHLTEMCAL.so");
	system->LoadComponentLibraries("libAliHLTMUON.so");
	
	AliHLTConfiguration reader("source", "MonitorHomerPublisher", "", "");
	AliHLTLogging log;
	
	switch (configId)
	{
	case kGlobalHistoChainConfig:
		{
		  const char* stdCuts="Track_TPCclus>0&&Track_DCAr<7&&Track_DCAr>-7&&Track_pt>0.3&&Track_eta<0.9&&Track_eta>-0.9";
		  TString config;
		  config += "-interval 20 -maxentries 250000 -maxmemory 4000000";
		  config += " -histogram TrackPt(100,0,10) -size 1000 -expression Track_pt -cut Track_TPCclus>0";
		  config += " -histogram TrackPhi(180,0,360) -size 1000 -expression Track_phi -cut Track_TPCclus>0";
		  config += " -histogram TrackMultiplicity(60,1,600) -size 1000 -expression trackcount";
		  config += " -histogram TrackMultiplicityTrend -size 1000 -expression trackcount:event";
		  config += Form(" -histogram TrackMultiplicityPrimary(60,1,600) -size 1000 -expression Sum$(%s)",stdCuts);
		  config += Form(" -histogram TrackMultiplicityPrimaryTrend -size 1000 -expression Sum$(%s):event",stdCuts);
		  config += Form(" -histogram TrackMultiplicityBackground(60,1,600) -size 1000 -expression Sum$(!(%s))",stdCuts);
		  config += Form(" -histogram TrackMultiplicityBackgroundTrend -size 1000 -expression Sum$(!(%s)):event",stdCuts);
		  config += Form(" -histogram TrackEta(100,-2,2) -size 1000 -expression Track_eta -cut %s",stdCuts);
		  config += Form(" -histogram TrackTPCclus(200,0,200) -size 1000 -expression Track_TPCclus -cut %s",stdCuts);
		  config += " -histogram TrackITSclus(7,-0.5,6.5) -size 1000 -expression Track_ITSclus";
		  config += Form(" -histogram TrackTheta(90,0,180) -size 1000 -expression Track_theta -cut %s",stdCuts);
		  config += " -histogram TrackDCAr(200,-50,50) -size 1000 -expression Track_DCAr -cut Track_TPCclus>0";
		  config += " -histogram TrackCharge(7,-1.75,1.75) -size 1000 -expression Track_charge -cut Track_TPCclus>0";
		  config += " -histogram VertexXY(100,-2.5,2.5,100,-2.5,2.5) -size 1000 -expression vertexY:vertexX -cut nContributors>3 -opt colz";
		  config += " -histogram VertexX(100,-2.5,2.5) -size 1000 -expression vertexX -cut nContributors>3";
		  config += " -histogram VertexY(100,-2.5,2.5) -size 1000 -expression vertexY -cut nContributors>3";
		  config += " -histogram VertexZ(600,-30,30) -size 1000 -expression vertexZ -cut nContributors>3";
		  config += " -histogram VertexTrendX -size 1000 -expression vertexX:event -cut nContributors>3";
		  config += " -histogram VertexTrendY -size 1000 -expression vertexY:event -cut nContributors>3";
		  AliHLTConfiguration tpcGlobalHisto("tpcGlobalHisto", "GlobalHisto", "source", config.Data());
		  listOfChains += " tpcGlobalHisto";
		}
		break;
	case kCompressionHistoChainConfig:
		{
		  AliHLTConfiguration blockfilter("blockfilter", "BlockFilter", "source", "-datatype 'HWCLUST1' 'TPC ' -datatype 'REMCLSCM' 'TPC ' -datatype 'COMPDESC' 'TPC '");
		AliHLTConfiguration compressionHisto("compressionHisto", "TPCDataCompressorMonitor", "blockfilter", "-pushback-period=20");
		listOfChains += " compressionHisto";
		}
		break;
	case kCorrelationMinBiasHistoChainConfig:
		{
		AliHLTConfiguration correlationHisto1("correlationHisto1", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CPBI1 -addTrigger CPBI2");
		AliHLTConfiguration renameHisto1("renameHisto1", "ObjectRenamer", "correlationHisto1", "-suffix _minbias");
		listOfChains += " renameHisto1";
		}
		break;
	case kCorrelationSemiCentralHistoChainConfig:
		{
		AliHLTConfiguration correlationHisto2("correlationHisto2", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CVLN");
		AliHLTConfiguration renameHisto2("renameHisto2", "ObjectRenamer", "correlationHisto2", "-suffix _semicentral");
		listOfChains += " renameHisto2";
		}
		break;
	case kCorrelationCentralHistoChainConfig:
		{
		AliHLTConfiguration correlationHisto3("correlationHisto3", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CVHN");
		AliHLTConfiguration renameHisto3("renameHisto3", "ObjectRenamer", "correlationHisto3", "-suffix _central");
		listOfChains += " renameHisto3";
		}
		break;
	case kCorrelationAllTriggerHistoChainConfig:
		{
		AliHLTConfiguration correlationHisto4("correlationHisto4", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CPBI1 -addTrigger CPBI2 -addTrigger CVLN -addTrigger CVHN");
		listOfChains += " correlationHisto4";
		}
		break;
	case kEmcalHistoChainConfig:
		{
		AliHLTConfiguration emcalHisto("emcalHisto", "EmcalClusterMonitor", "source", "-pushback-period=20");
		listOfChains += " emcalHisto";
		}
		break;
	case kMuonHistoChainConfig:
		{
		AliHLTConfiguration muonHisto("muonHisto", "MUONClusterHistogrammer", "source", "-pushback-period=20");
		listOfChains += " muonHisto";
		}
		break;
	case kCorrelationsHistoChainConfig:
		{
		AliHLTConfiguration correlationHisto1("correlationHisto1", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CPBI1 -addTrigger CPBI2");
		AliHLTConfiguration renameHisto1("renameHisto1", "ObjectRenamer", "correlationHisto1", "-suffix _minbias");
		listOfChains += " renameHisto1";
		AliHLTConfiguration correlationHisto2("correlationHisto2", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CVLN");
		AliHLTConfiguration renameHisto2("renameHisto2", "ObjectRenamer", "correlationHisto2", "-suffix _semicentral");
		listOfChains += " renameHisto2";
		AliHLTConfiguration correlationHisto3("correlationHisto3", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CVHN");
		AliHLTConfiguration renameHisto3("renameHisto3", "ObjectRenamer", "correlationHisto3", "-suffix _central");
		listOfChains += " renameHisto3";
		AliHLTConfiguration correlationHisto4("correlationHisto4", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CPBI1 -addTrigger CPBI2 -addTrigger CVLN -addTrigger CVHN");
		listOfChains += " correlationHisto4";
		}
		break;
		
	default:
		log.LoggingVarargs(kHLTLogError, "", FUNCTIONNAME(), __FILE__, __LINE__,
					"Could not configure the monitoring chain, unknown configId = %d.",
					configId
				);
		return;
	}
	
        //receive blocks via ZMQ
        AliHLTConfiguration zmqSource("ZMQsource","ZMQsource",listOfChains.Data(),"in=REQ+tcp://ecs1.internal:60321 ZMQneverBlock=0 ZMQrequestTimeout=1000000");
        listOfChains += " ZMQsource";
	
        //caching relay
        AliHLTConfiguration monitorRelay("monitorRelay", "MonitoringRelay", listOfChains.Data(), "-check-object"); 
        
        //send blocks via ZMQ
        AliHLTConfiguration zmqSink("ZMQsink","ZMQsink", "monitorRelay", "out=PUB>tcp://ecs1:60212");

	AliHLTConfiguration writer("sink", "MonitorMemoryWriter", "monitorRelay ZMQsink", "");
	
	system->BuildTaskList("sink");
}

////////////////////////////////////////////////////////////////////////////////
//      The following code need not be changed if you are not an expert.      //
////////////////////////////////////////////////////////////////////////////////

/// The server socket to which DQM clients connect.
TServerSocket* gServerSocket = NULL;

/// The HOMER library manager through which instances of AliHLTHOMERReader and
/// AliHLTHOMERWriter are created.
AliHLTHOMERLibManager* gHomerLib = NULL;

/// The list of active client TSocket objects.
TObjArray gClientList;

/// The current active HOMER reader instance.
AliHLTHOMERReader* gHomerReader = NULL;


void DeleteGlobalObjects();


/**
 * Creates the global objects for the server socket and HOMER library manager.
 * \return true if all the global objects could be created and false in error.
 */
bool CreateGlobalObjects(int serverPort)
{
	DeleteGlobalObjects();  // delete anything that already exists.
	
	// Create the server socket.
	gServerSocket = new TServerSocket(serverPort, kTRUE);
	if (gServerSocket == NULL)
	{
		cerr << "ERROR: Could not create server socket on port " << serverPort << "." << endl;
		return false;
	}
	gServerSocket->SetOption(kNoBlock, 1);
	if (not gServerSocket->IsValid())
	{
		cerr << "ERROR: Could not start server socket on port " << serverPort << "." << endl;
		return false;
	}
	
	gHomerLib = new AliHLTHOMERLibManager;
	if (gHomerLib == NULL)
	{
		cerr << "ERROR: Could not create HOMER library manager." << endl;
		return false;
	}
	
	gClientList.SetOwner(kTRUE);
	
	return true;
}


/**
 * Deletes the global objects created in CreateGlobalObjects().
 */
void DeleteGlobalObjects()
{
	gClientList.Clear();
	if (gServerSocket != NULL)
	{
		delete gServerSocket;
		gServerSocket = NULL;
	}
	if (gHomerLib != NULL)
	{
		delete gHomerLib;
		gHomerLib = NULL;
	}
}


/**
 * Tries to open a new connection to the HOMER interface on the given
 * port and address if not already connected.
 * \param hostName The name of the host machine to connect to.
 * \param port  The port number on which the HOMER interface is listening.
 */
void OpenHomerConnection(const char* hostName, int port)
{
	if (gHomerReader != NULL) return;
	
	gHomerReader = gHomerLib->OpenReader(hostName, port);
	if (gHomerReader == NULL or gHomerReader->GetConnectionStatus() != 0)
	{
		//cerr << "ERROR: Could not open HOMER connection to " << hostName << ":" << port << endl;
		return;
	}
	else
	{
		cout << "Opened HOMER connection to " << hostName << ":" << port << endl;
	}
}


/**
 * This class is used to feed data blocks from a HOMER reader into a HLT chain
 * running under AliHLTSystem.
 * A global instance of the current AliHLTHOMERReader to used should be set with
 * the HomerReader() method. The class does not fetch the next event automatically.
 * It is expected that the run is stepped through one event at a time and the
 * outer loop must call AliHLTHOMERReader::ReadNextEvent() directly.
 * 
 * \note This class is not thread safe.
 */
class AliHLTMonitorHomerPublisher : public AliHLTDataSource
{
public:
	AliHLTMonitorHomerPublisher() : AliHLTDataSource() {}
	virtual ~AliHLTMonitorHomerPublisher() {}

	virtual const char* GetComponentID() { return "MonitorHomerPublisher"; }
	virtual AliHLTComponentDataType GetOutputDataType() { return kAliHLTAnyDataType; }
	
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
	{
		constBase = 128*1024*1024;
		inputMultiplier = 0;
	}
	
	virtual AliHLTComponent* Spawn() { return new AliHLTMonitorHomerPublisher; }
	
	static const AliHLTHOMERReader* HomerReader() { return fgHomerReader; }
	static void HomerReader(const AliHLTHOMERReader* reader) { fgHomerReader = reader; }

protected:
	
	virtual int DoInit(int /*argc*/, const char** /*argv*/) { return 0; }
	virtual int DoDeinit() { return 0; }

	virtual int GetEvent(
			const AliHLTComponentEventData& evtData,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr, 
			AliHLTUInt32_t& size,
			AliHLTComponentBlockDataList& outputBlocks
		);

	using AliHLTDataSource::GetEvent;
	
private:
	// Do not allow copying of this class.
	AliHLTMonitorHomerPublisher(const AliHLTMonitorHomerPublisher&);
	AliHLTMonitorHomerPublisher& operator=(const AliHLTMonitorHomerPublisher&);
	
	static const AliHLTHOMERReader* fgHomerReader;

	ClassDef(AliHLTMonitorHomerPublisher, 0)
};


ClassImp(AliHLTMonitorHomerPublisher)


const AliHLTHOMERReader* AliHLTMonitorHomerPublisher::fgHomerReader = NULL;


/**
 * This method just copies the data blocks from fgHomerReader and injects them into
 * the HLT chain. The fgHomerReader HOMER reader should already have the ReadNextEvent()
 * method called successfully.
 */
int AliHLTMonitorHomerPublisher::GetEvent(
		const AliHLTComponentEventData& /*evtData*/,
		AliHLTComponentTriggerData& /*trigData*/,
		AliHLTUInt8_t* outputPtr,
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks
	)
{
	if (fgHomerReader == NULL) return 0;
	AliHLTUInt32_t capacity = size;
	AliHLTUInt8_t* currentOutput = outputPtr;
	size = 0;
	for (unsigned long i = 0; i < fgHomerReader->GetBlockCnt(); i++)
	{
		if (size + fgHomerReader->GetBlockDataLength(i) > capacity)
		{
			return -ENOMEM;
		}
		homer_uint64 type = AliHLTOUT::ByteSwap64(fgHomerReader->GetBlockDataType(i));
		homer_uint32 origin = AliHLTOUT::ByteSwap32(fgHomerReader->GetBlockDataOrigin(i));
		AliHLTComponentDataType dt;
		memcpy(&dt.fID, &type, sizeof(type));
		memcpy(&dt.fOrigin, &origin, sizeof(origin));
		
		HLTInfo("Injecting new data block %s", AliHLTComponent::DataType2Text(dt).c_str());
		
		memcpy(currentOutput, fgHomerReader->GetBlockData(i), fgHomerReader->GetBlockDataLength(i));
		
		AliHLTComponentBlockData bd;
		FillBlockData(bd);
		bd.fPtr = outputPtr;
		bd.fOffset = currentOutput - outputPtr;
		bd.fSize = fgHomerReader->GetBlockDataLength(i);
		bd.fDataType = dt;
		bd.fSpecification = fgHomerReader->GetBlockDataSpec(i);
		outputBlocks.push_back(bd);
		
		size += fgHomerReader->GetBlockDataLength(i);
		currentOutput += fgHomerReader->GetBlockDataLength(i);
	}
	return 0;
}


class AliHLTObjectRenamer : public AliHLTProcessor
{
public:
	AliHLTObjectRenamer() : AliHLTProcessor(), fSuffix("") {}
	virtual ~AliHLTObjectRenamer() {}

	virtual const char* GetComponentID() { return "ObjectRenamer"; }
	virtual void GetInputDataTypes( vector<AliHLTComponentDataType>& list) { list.push_back(kAliHLTAnyDataType); } 
	virtual AliHLTComponentDataType GetOutputDataType() { return kAliHLTAnyDataType; }
	
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
	{
		constBase = 1024*1024;
		inputMultiplier = 1;
	}
	
	virtual AliHLTComponent* Spawn() { return new AliHLTObjectRenamer; }

protected:
	
	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit() { return 0; }

	virtual int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

	using AliHLTProcessor::DoEvent;
	
private:
	// Do not allow copying of this class.
	AliHLTObjectRenamer(const AliHLTObjectRenamer&);
	AliHLTObjectRenamer& operator=(const AliHLTObjectRenamer&);
	
	TString fSuffix;   /// Suffix to appent to object names.

	ClassDef(AliHLTObjectRenamer, 0)
};


ClassImp(AliHLTObjectRenamer)


int AliHLTObjectRenamer::DoInit(int argc, const char** argv)
{
	fSuffix = "";
	bool suffixSet = false;
	for (int i = 0; i < argc; i++)
	{
		if (strcmp( argv[i], "-suffix" ) == 0)
		{
			if (suffixSet)
			{
				HLTWarning("The suffix has already been specified."
					" Will replace previous value given by -suffix."
				);
			}
			if ( argc <= i+1 )
			{
				HLTError("The suffix value was not specified.");
				return -EINVAL;
			}
			fSuffix = argv[i+1];
			i++;
			continue;
		}
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}
	
	return 0;
}


template <typename ListType>
TObject* CloneList(const TObject* obj, const TString& suffix)
{
	if (obj->IsA() == ListType::Class())
	{
		const ListType* oldArray = static_cast<const ListType*>(obj);
		ListType* newArray = new ListType;
		newArray->SetOwner(kTRUE);
		TIter next(oldArray);
		const TObject* tmpObj = NULL;
		while ((tmpObj = next()) != NULL)
		{
			TString name = tmpObj->GetName();
			name += suffix;
			TObject* newObj = tmpObj->Clone(name);
			newArray->Add(newObj);
		}
		return newArray;
	}
	return NULL;
}


int AliHLTObjectRenamer::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
	const TObject* obj = NULL;
	for (obj = GetFirstInputObject(); obj != NULL; obj = GetNextInputObject())
	{
		TObject* newObj = NULL; CloneList<TObjArray>(obj, fSuffix);
		if (newObj == NULL) newObj = CloneList<TList>(obj, fSuffix);
		if (newObj == NULL)
		{
			TString name = obj->GetName();
			name += fSuffix;
			newObj = obj->Clone(name);
		}
		PushBack(newObj, GetDataType(), GetSpecification());
		delete newObj;
	}
	return 0;
}


/**
 * This class is used to collect the data blocks generated by monitoring components
 * into temporary memory buffers which can then be fetched through GetOutputBlocks()
 * and handled in any manner appropriate. In this case they are compiled into
 * a HOMER format image to be sent to any connected DQM clients.
 * The data blocks are constantly accumulated in the fgOutputBlocks array, even
 * between events. So the outer loop must step through events one at a time using
 * the AliHLTSystem::Run() method and call ClearOutputBlocks() on each iteration
 * before calling the Run() method.
 * 
 * \note This class is not thread safe.
 */
class AliHLTMonitorMemoryWriter : public AliHLTDataSink
{
public:
	AliHLTMonitorMemoryWriter() : AliHLTDataSink() {}
	virtual ~AliHLTMonitorMemoryWriter() { ClearOutputBlocks(); }

	virtual const char* GetComponentID() { return "MonitorMemoryWriter"; }
	virtual void GetInputDataTypes( vector<AliHLTComponentDataType>& list) { list.push_back(kAliHLTAnyDataType); } 
	virtual AliHLTComponent* Spawn() { return new AliHLTMonitorMemoryWriter; }
	
	static const std::vector<AliHLTComponentBlockData>& GetOutputBlocks() { return fgOutputBlocks; }
	static void ClearOutputBlocks();
	static AliHLTUInt64_t CurrentEventId() { return fgEventId; }

protected:
	
	virtual int DoInit(int /*argc*/, const char** /*argv*/) { return 0; }
	virtual int DoDeinit() { return 0; }

	virtual int DumpEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks, 
			AliHLTComponentTriggerData& trigData
		);

	using AliHLTDataSink::DumpEvent;
	
private:
	// Do not allow copying of this class.
	AliHLTMonitorMemoryWriter(const AliHLTMonitorMemoryWriter&);
	AliHLTMonitorMemoryWriter& operator=(const AliHLTMonitorMemoryWriter&);
	
	static std::vector<AliHLTComponentBlockData> fgOutputBlocks;
	static AliHLTUInt64_t fgEventId;

	ClassDef(AliHLTMonitorMemoryWriter, 0)
};


ClassImp(AliHLTMonitorMemoryWriter)


std::vector<AliHLTComponentBlockData> AliHLTMonitorMemoryWriter::fgOutputBlocks;
AliHLTUInt64_t AliHLTMonitorMemoryWriter::fgEventId = 0x0;

/**
 * This method simply collects all input objects into the fgOutputBlocks buffer.
 * \return 0 on success and non-zero on failure.
 */
int AliHLTMonitorMemoryWriter::DumpEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks, 
		AliHLTComponentTriggerData& /*trigData*/
	)
{
	fgEventId = evtData.fEventID;
	AliHLTMonitorMemoryWriter::ClearOutputBlocks();
	for (unsigned int i = 0; i < evtData.fBlockCnt; ++i)
	{
		AliHLTComponentBlockData bd;
		FillBlockData(bd);
		bd.fPtr = new char[blocks[i].fSize];
		bd.fOffset = 0;
		bd.fSize = blocks[i].fSize;
		bd.fDataType = blocks[i].fDataType;
		bd.fSpecification = blocks[i].fSpecification;
		
		char* sourceBuffer = reinterpret_cast<char*>(blocks[i].fPtr);
		memcpy(bd.fPtr, sourceBuffer + blocks[i].fOffset, bd.fSize);
		
		fgOutputBlocks.push_back(bd);
	}
	return 0;
}


/**
 * Deletes all memory allocated for the data blocks collected in the fgOutputBlocks
 * buffers in the DumpEvent() method. 
 */
void AliHLTMonitorMemoryWriter::ClearOutputBlocks()
{
	for (size_t i = 0; i < fgOutputBlocks.size(); ++i)
	{
		char* buffer = reinterpret_cast<char*>( fgOutputBlocks[i].fPtr );
		delete [] buffer;
	}
	fgOutputBlocks.clear();
}


/**
 * Read the next event from the HOMER reader interface.
 * \param timeout  The maximum amount of time to wait for the next event in micro-seconds.
 */
void ReadEventFromHomer(int timeout = 1000000)
{
	AliHLTMonitorHomerPublisher::HomerReader(NULL);
	//outblocks are cleared in the AliHLTMonitorMemoryWriter::DumpEvent
	//method, otherwise they would be empty most of the time.
	//Getting the blocks would only succeed by chance.
	//AliHLTMonitorMemoryWriter::ClearOutputBlocks();
	if (gHomerReader == NULL) return;
	
	int result = gHomerReader->ReadNextEvent(timeout);
	if (result != 0)
	{
		if (result != ETIMEDOUT)
		{
		  //cout << "Closed HOMER connection." << endl;
			gHomerLib->DeleteReader(gHomerReader);
			gHomerReader = NULL;
		}
		return;
	}
	
	AliHLTLogging log;
	log.LoggingVarargs(kHLTLogInfo, "", FUNCTIONNAME(), __FILE__, __LINE__,
			"Received event 0x%llX with %d data blocks.",
			gHomerReader->GetEventID(), gHomerReader->GetBlockCnt()
		);
	AliHLTMonitorHomerPublisher::HomerReader(gHomerReader);
}


/**
 * Fetches needed run information from the ECS file on protal-ecs1.
 * \param runNumber  Filled with the run number sent by ECS.
 * \param ctpTriggerClasses  Filled with the CTP trigger class string sent by ECS.
 * \param detectorList  Detectors taking part in the run.
 */
void FetchRunInformation(int& runNumber, TString& ctpTriggerClasses, TString& detectorList)
{
  //const char* cmd1 = "ssh ecs0 'cat $ECS_PROXY_RUNDIR/log/ECSProxy*.log | grep \"ECS-Proxy received parameter\" | sort | tail -n11 | sed -e \"s/.*ECS-Proxy received parameter //\"'";
  //const char* cmd1 = "ssh ecs0 'tac $(ls -1t $ECS_PROXY_RUNDIR/log/ECSProxy*.log | head -1) | sed -n \"1,/.* ENGAGE.*/p\" | sed \"s/^.*received parameter //g\"'";
  const char* cmd1 = "ssh ecs0 'tac $(ls -1t /opt/HLT/dist/ecs/log/ECSProxy*.log | head -1) | sed -n \"1,/.*\\<ENGAGE\\>.*/p\"|grep \"received parameter\"|sed \"s/^.*received parameter //g\"|sed \"s/ from ECS.//\"'";
  printf("executing: %s\n", cmd1);
	FILE* pipe = gSystem->OpenPipe(cmd1, "r");
	if (pipe == NULL)
	{
		runNumber = -1;
		ctpTriggerClasses = "";
		detectorList = "";
		return;
	}
	char buf[1024*1024];
	size_t pos = fread(buf, 1, sizeof(buf), pipe);
	buf[pos] = '\0';
	gSystem->ClosePipe(pipe);
	TString str = buf;
	TString token;
	Ssiz_t from = 0;
	while (str.Tokenize(token, from, "\n"))
	{
	  token.ReplaceAll(" from ECS.","");
		if (token.Contains("RUN_NUMBER"))
		{
			token.ReplaceAll("RUN_NUMBER :  ", "");
			runNumber = token.Atoi();
		}
		else if (token.Contains("CTP_TRIGGER_CLASS"))
		{
		  token.ReplaceAll("CTP_TRIGGER_CLASS :  ", "CTP_TRIGGER_CLASS=");
		  ctpTriggerClasses = token;
		}
		else if (token.Contains("DETECTOR_LIST"))
		{
		  token.ReplaceAll("DETECTOR_LIST :  ", "");
		  detectorList = token;
		}
	}
}


/**
 * Fetches the CTP trigger classes from the global trigger object.
 * \returns the CTP trigger classes and 0x0 if none was found.
 */
AliHLTTriggerMask_t GetTriggerClassesForCurrentEvent()
{
  AliHLTTriggerMask_t trigMask = 0;
  if (gHomerReader == NULL) return trigMask;
  AliHLTLogging log;
  long long unsigned int eventID = gHomerReader->GetEventID();

  for (unsigned long i = 0; i < gHomerReader->GetBlockCnt(); i++)
    {
      homer_uint64 type = AliHLTOUT::ByteSwap64(gHomerReader->GetBlockDataType(i));
      homer_uint32 origin = AliHLTOUT::ByteSwap32(gHomerReader->GetBlockDataOrigin(i));
      AliHLTComponentDataType dt;
      memcpy(&dt.fID, &type, sizeof(type));
      memcpy(&dt.fOrigin, &origin, sizeof(origin));
      if (dt == kAliHLTDataTypeGlobalTrigger)
	{
	  AliHLTHOMERBlockDesc* block = new AliHLTHOMERBlockDesc();
	  block->SetBlock(
			  const_cast<void*>(gHomerReader->GetBlockData(i)), 
			  gHomerReader->GetBlockDataLength(i),
			  origin,
			  type,
			  (ULong_t) gHomerReader->GetBlockDataSpec(i)
			  );

	  TObject* obj = block->GetTObject();
	  AliHLTGlobalTriggerDecision* decision;
	  if(obj && (decision=dynamic_cast<AliHLTGlobalTriggerDecision*>(obj) ) ) 
	    {
	      const TObject* ctpobj = decision->InputObjects().FindObject("AliHLTCTPData");
	      const AliHLTCTPData* ctpdata = (ctpobj != NULL ? dynamic_cast<const AliHLTCTPData*>(ctpobj) : NULL);
	      if (ctpdata != NULL)
		{
		trigMask = ctpdata->Triggers();
		delete block;
		break;
		}
	      else
		{
		  log.LoggingVarargs(kHLTLogError, "", FUNCTIONNAME(), __FILE__, __LINE__,
				     "Could not find CTP data object in global trigger for event %lld (0x%016llX).",
				     eventID, eventID
				     );
		}
	    }
	  delete block;
	}
    }
  log.LoggingVarargs(kHLTLogInfo, "", FUNCTIONNAME(), __FILE__, __LINE__,
		     "Using CTP trigger class bits 0x%016llx for event %lld (0x%016llX).",
		     trigMask, eventID, eventID
		     );

  return trigMask;
}


/**
 * Processes the event data read in ReadEventFromHomer() by passing it through a
 * HLT chain using AliHLTSystem. The chain is configured in Configure() and only
 * single stepped in this routine.
 * Initial run parameters are taken from the RunManager files generated on portal-ecs1
 * every second. Thus if a new run is started the stop chain command is sent instead
 * of a normal processing step and the next time this routine is called the
 * AliHLTSystem::Run() will generate the Start-Of-Run sequence as if a new run
 * is started.
 * \param checkPeriod  The time in seconds between checks of the run information from ECS.
 */
int ProcessEvent(double checkPeriod = 1)
{
	AliHLTSystem* system = AliHLTPluginBase::GetInstance();
	
	AliHLTLogging log;
	static int previousRunNumber = -1;
	static AliHLTEventID_t lastEventID = -1;
	static AliHLTUInt32_t participatingDetectors = 0;
	
	static double lastCheckOfRunInfo = 0;
	double now = TTimeStamp().AsDouble();
	if (now > lastCheckOfRunInfo + checkPeriod)
	{
		lastCheckOfRunInfo = now;
		int runNumber = -1;
		TString ctpTriggerClasses;
		TString detectorList;
		FetchRunInformation(runNumber, ctpTriggerClasses, detectorList);
		std::cout << runNumber << std::endl
			  << ctpTriggerClasses << std::endl
			  << detectorList << std::endl;
		if (previousRunNumber != runNumber and previousRunNumber != -1){
		  log.LoggingVarargs(kHLTLogError, "", FUNCTIONNAME(), __FILE__, __LINE__,
				     "Run changed. Aborting processing. Restart from scratch..."
				     );
		  return -1;
		}
		if (previousRunNumber != runNumber and runNumber != -1)
		{
			// Send stop command.
			system->Run(0, 1);

			// Initialise new run.
			log.LoggingVarargs(kHLTLogInfo, "", FUNCTIONNAME(), __FILE__, __LINE__,
					"New run %d.",
					runNumber
				);
			log.LoggingVarargs(kHLTLogInfo, "", FUNCTIONNAME(), __FILE__, __LINE__,
					"Setting new ECS options: '%s'",
					ctpTriggerClasses.Data()
				);
			log.LoggingVarargs(kHLTLogInfo, "", FUNCTIONNAME(), __FILE__, __LINE__,
					"Participating detectors: '%s'",
					detectorList.Data()
				);
			TString opt = "ECS=";
			opt += ctpTriggerClasses;
			system->ScanOptions(opt.Data());
			previousRunNumber = runNumber;
			AliCDBManager::Instance()->SetRun(runNumber);
			AliHLTMisc::Instance().InitMagneticField();		
			detectorList.ReplaceAll(",", " ");
			participatingDetectors = AliDAQ::DetectorPatternOffline(detectorList) xor AliDAQ::DetectorPatternOffline("DAQ_TEST");
		}
	}
	
	if (AliCDBManager::Instance()->GetRun() == -1)
	{
		log.LoggingVarargs(kHLTLogError, "", FUNCTIONNAME(), __FILE__, __LINE__,
				"Could not fetch or set the run number from ECS parameters. Aborting processing."
			);
		return -1;
	}
	
	if (gHomerReader == NULL) return -1;
	
	int nofEvents = 1;
	int stop = 0;
	AliHLTTriggerMask_t trgMask = GetTriggerClassesForCurrentEvent();
	AliHLTUInt32_t timestamp = AliHLTUInt32_t( TTimeStamp().AsDouble() );
	AliHLTUInt32_t eventType = AliHLTUInt32_t( gHomerReader->GetEventType() );
	
	if (lastEventID != gHomerReader->GetEventID())
	{
		system->Run(nofEvents, stop, trgMask, timestamp, eventType, participatingDetectors);
	}
	lastEventID = gHomerReader->GetEventID();
	return 0;
}


/**
 * Prepares the HOMER image buffer for writing to the DQM clients using AliHLTHOMERWriter.
 */
void WriteOutputBuffer(char* buffer, size_t& size)
{
	AliHLTHOMERWriter* writer = gHomerLib->OpenWriter();
	
	homer_uint64 eventType = 0x0;
	
	// Loop through all the data blocks collected in the monitor memory writer
	// and add them to the HOMER writer image being created.
	const std::vector<AliHLTComponentBlockData>& outputBlocks = AliHLTMonitorMemoryWriter::GetOutputBlocks();
	for (size_t i = 0; i < outputBlocks.size(); ++i)
	{
		// If we find the data block marking the data type then
		// record the type of the event but discard that data block.
		// We dont want to forward it to the DQM clients.
		if (outputBlocks[i].fDataType == kAliHLTDataTypeEvent)
		{
			eventType = outputBlocks[i].fSpecification;
			continue;
		}
		
		homer_uint64 homerHeader[kCount_64b_Words];
		HOMERBlockDescriptor homerDescriptor(homerHeader);
		memset( homerHeader, 0, sizeof(homer_uint64)*kCount_64b_Words );
		homerDescriptor.Initialize();
		homer_uint64 id=0;
		homer_uint64 origin=0;
		memcpy(&id, outputBlocks[i].fDataType.fID, sizeof(homer_uint64));
		memcpy(((AliHLTUInt8_t*)&origin)+sizeof(homer_uint32), outputBlocks[i].fDataType.fOrigin, sizeof(homer_uint32));
		homerDescriptor.SetType(AliHLTOUT::ByteSwap64(id));
		homerDescriptor.SetSubType1(AliHLTOUT::ByteSwap64(origin));
		homerDescriptor.SetSubType2(outputBlocks[i].fSpecification);
		homerDescriptor.SetBlockSize(outputBlocks[i].fSize);
		homerDescriptor.SetByteOrder(kHOMERNativeByteOrder);
		writer->AddBlock(homerHeader, outputBlocks[i].fPtr);
	}
	
	// Now generate the HOMER image buffer with the right run number and event type.
	size_t neededSize = writer->GetTotalMemorySize();
	if (neededSize <= size)
	{
		homer_uint64 eventNr = AliHLTMonitorMemoryWriter::CurrentEventId();
		homer_uint64 statusFlags = 0x0;  // No status flags.
		homer_uint64 nodeID = 0x0;   // No node ID.
		writer->Copy(buffer, eventType, eventNr, statusFlags, nodeID);
		size = neededSize;
	}
	else
	{
		cerr << "ERROR: Output buffer too small. Size is " << size
			<< " bytes, but need " << neededSize << " bytes." << endl;
		size = 0;
	}
	
	gHomerLib->DeleteWriter(writer);
}


/**
 * Accept new incoming DQM client connections on the server socket.
 */
void AcceptNewConnections()
{
  //loop!
  while(true) {
	TSocket* client = gServerSocket->Accept();
	if (client == reinterpret_cast<TSocket*>(-1)) return;  // no pending client.
	if (client == NULL)
	{
		cerr << "ERROR: Could not accept new client connection." << endl;
		return;
	}
	client->SetOption(kNoBlock, 1);
	gClientList.Add(client);
  }
}


// Structure used to keep track of transmission progress to DQM clients.
struct AliHLTConnectionInfo
{
	size_t fWrittenSize;   // Number of bytes of the HOMER image buffer written.
	size_t fReceivedSize;  // Number of bytes of commands received from the client.
	char fBuffer[256];     // The buffer to which incoming client commands are received.
};


/**
 * Writes the prepared HOMER image buffer to the listening DQM clients over TCP connections.
 * The list of clients is stored in gClientList as a list of TSocket objects.
 * The output interface replicates a very light weight version of the HOMER interface.
 * \param buffer  The HOMER image buffer to write.
 * \param size    The size in bytes of the buffer.
 */
void WriteToClients(const char* buffer, size_t size)
{
	double startTime = TTimeStamp().AsDouble();
	double timeoutPeriod = 2.;
	
	std::vector<AliHLTConnectionInfo> connectionInfo;
	for (int i = 0; i < gClientList.GetEntriesFast(); ++i)
	{
		AliHLTConnectionInfo newconn;
		memset(&newconn, 0x0, sizeof(newconn));
		connectionInfo.push_back(newconn);
	}
	
	bool allDone = false;
	while (not allDone)
	{
		int totalClients = 0;
		int clientsDone = 0;
		for (int i = 0; i < gClientList.GetEntriesFast(); ++i)
		{
			TSocket* client = static_cast<TSocket*>( gClientList.UncheckedAt(i) );
			if (client == NULL) continue;  // might be NULL because the connection was lost.
			++totalClients;
			AliHLTConnectionInfo& info = connectionInfo[i];
			if (info.fReceivedSize >= 8 and info.fWrittenSize >= size) ++clientsDone;
			
			if (info.fReceivedSize < 8)
			{
				// Read the client socket for command.
				int bytesRead = client->RecvRaw(&info.fBuffer[info.fReceivedSize], sizeof(info.fBuffer) - info.fReceivedSize, kDontBlock);
				if (bytesRead >= 0)
				{
					info.fReceivedSize += bytesRead;
				}
				else if (bytesRead != -4) // was it a non-blocking operation or an error
				{
					gClientList[i] = NULL;
					delete client;
				}
			}
			
			bool commandOk = false;
			if (info.fReceivedSize >= 8)
			{
				// Check that the command was correct.
				char tmpbuf[1024];
				memset(tmpbuf, 0x0, sizeof(tmpbuf));
				memcpy(tmpbuf, &info.fBuffer[0], info.fReceivedSize);
				TString msg(tmpbuf);
				if (msg == "GET ONE\n" or msg == "MOD BIN\nGET ONE\n")
				{
					commandOk = true;
				}
				else
				{
					gClientList[i] = NULL;
					delete client;
				}
			}
			
			if (commandOk and info.fWrittenSize < size)
			{
				// Write the output to the client.
				if (info.fWrittenSize == 0)
				{
					UInt_t sizeOfStream = htonl(size);
					client->SendRaw(&sizeOfStream, sizeof(UInt_t));
				}
				int bytesWritten = client->SendRaw(buffer + info.fWrittenSize, size - info.fWrittenSize);
				if (bytesWritten >= 0)
				{
					info.fWrittenSize += bytesWritten;
				}
				else if (bytesWritten != -4) // was it a non-blocking operation or an error
				{
					gClientList[i] = NULL;
					delete client;
				}
			}
		}
		
		if (totalClients == clientsDone) allDone = true;
		//if (TTimeStamp().AsDouble() > startTime + timeoutPeriod) break;
		//gSystem->Sleep(10);
	}
	
	gClientList.Compress();
}


/**
 * The entry point routine for the MonitorSandbox program.
 * \param hostName  The name of the machine to connect to where the HLT online
 *     data can be fetched. This is normally a TCPDumpSubscriber.
 * \param port  The port number on which HOMER is listening to which we connect.
 * \param serverPort  The server port number to use to which DQM clients connect.
 * \param pollPeriod  The time in seconds between polls of the HOMER interface when
 *     fetching the next event.
 * \param cdbpath The path to the CDB to use.
 * \param configId  Identifier for the configuration to use.
 */
void MonitorSandbox(
		const char* hostName = "localhost",
		int port = 49000,
		int serverPort = 49001,
		double pollPeriod = 0.1,  // seconds
		const char* cdbpath = "local://$ALICE_ROOT/OCDB",
		int configId = 1,
		double runInfoPollPeriod = 10  // seconds
	)
{
	AliCDBManager::Instance()->SetDefaultStorage(cdbpath);
	
	AliHLTSystem* system = AliHLTPluginBase::GetInstance();
	AliHLTComponentHandler* componentHandler = system->GetComponentHandler();
	componentHandler->AddComponent(new AliHLTMonitorHomerPublisher);
	componentHandler->AddComponent(new AliHLTMonitorMemoryWriter);
	componentHandler->AddComponent(new AliHLTObjectRenamer);
	
	if (not CreateGlobalObjects(serverPort)) return;
	size_t bufferCapacity = 1024*1024*128;
	char* buffer = new char[bufferCapacity];
	
	Configure(configId);
	
	do
	{
		double startTime =  TTimeStamp().AsDouble();
		
		OpenHomerConnection(hostName, port);
		ReadEventFromHomer();
		//Stop macro if run changes or on error
		//Mitigates memory leaks
		if( ProcessEvent(runInfoPollPeriod) < 0 )
		  break;
		size_t size = bufferCapacity;
		WriteOutputBuffer(buffer, size);
		AcceptNewConnections();
		WriteToClients(buffer, size);
		
		double endTime = TTimeStamp().AsDouble();
		double timeDiff = endTime - startTime;
		if (timeDiff < 1.0 and timeDiff < pollPeriod) gSystem->Sleep(int((pollPeriod - timeDiff)*1000));
	}
	while (gServerSocket->IsValid());
	
	delete [] buffer;
	DeleteGlobalObjects();
}
