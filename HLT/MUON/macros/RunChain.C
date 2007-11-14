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

/* This macro is used to run the dHLT component chain within the AliHLTSystem
 * framework, which is a simulation of the dHLT online system's output.
 * To run this macro you must be in the same directory as the galice.root
 * files or in the directory containing the rawXX/ directories produced from
 * AliRoot simulations.
 * Also make sure the LUTs are in the same directory or specify the directory
 * with the lutDir parameter option.
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "HLT/BASE/AliHLTSystem.h"
#include "HLT/BASE/AliHLTConfiguration.h"
#include "AliLog.h"
#include "TString.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

/* @param chainType  Specifies the type of chain to run. This can be one of the
 *     following:
 *       "full" - Run the full dHLT chain. (default)
 *       "ddlreco" - Run only the reconstruction of the DDL raw data up to hits
 *                   and trigger records.
 *       "tracker" - Run the tracker only using hits and trigger records from
 *                   AliRoot simulation or offline reconstruction.
 * @param firstEvent  The event number of the first event to process. (default = 0)
 * @param lastEvent  The event number of the last event to process. If this is
 *     less than firstEvent then it is set to firstEvent automatically. (default = -1)
 * @param output  Specifies the kind of output to generate. The options can be one
 *     of the following:
 *       "bin" - Generates all possible output and dumps it to binary files
 *                  using the FileWriter component.
 *       "root" - Generates all possible output and writes it to a ROOT file
 *                  using the ROOTFileWriter component.
 *       "tracks_bin" - Generates only track data and dumps it to binary files
 *                  using the FileWriter component. This option is equivalent to
 *                  'bin' if the chain type is 'ddlreco'.
 *       "tracks_root" - Generates only track data and writes it to a ROOT file
 *                  using the ROOTFileWriter component. This option is equivalent
 *                  to 'bin' if the chain type is 'ddlreco'.
 * @param dataSource  This is the data source from which to use hits and trigger
 *     records when we are running the tracker only chain. Note this parameter
 *     only makes sense if chainType = "tracker". The possible options are:
 *       "sim" - Take data from GEANT hits. (default)
 *       "rec" - Take data from offline reconstructed hits and local trigger
 *               objects.
 * @param logLevel  Specifies the amount of logging messages produced by the HLT
 *     system. The possible options are:
 *       "normal" - Shows warnings, errors and information.
 *       "max" - Shows everything including debug messages if they were compiled in.
 *       "min" - Shows only error messages.
 * @param lutDir  This is the directory in which the LUTs can be found.
 */
void RunChain(
		const char* chainType = "full",
		Int_t firstEvent = 0,
		Int_t lastEvent = -1,
		const char* output = "bin",
		const char* dataSource = "sim",
		const char* logLevel = "normal",
		const char* lutDir = "."
	)
{
	// Make sure that the lastEvent is greater than firstEvent.
	if (lastEvent < firstEvent)
		lastEvent = firstEvent;
	int eventCount = lastEvent - firstEvent + 1;
	
	bool buildDDLFilePubs = false;
	bool buildDDLRecoComps = false;
	bool buildSimDataPubs = false;
	bool buildRecDataPubs = false;
	bool buildTrackerComp = false;
	bool maxLogging = false;
	bool minLogging = false;
	bool useRootWriter = false;
	bool makeTracksOnly = false;
	
	// Parse the chainType, output, dataSource and logLevel option strings:
	TString outOpt = output;
	if (outOpt.CompareTo("bin", TString::kIgnoreCase) == 0)
	{
		useRootWriter = false;
	}
	else if (outOpt.CompareTo("root", TString::kIgnoreCase) == 0)
	{
		useRootWriter = true;
	}
	else if (outOpt.CompareTo("tracks_bin", TString::kIgnoreCase) == 0)
	{
		useRootWriter = false;
		makeTracksOnly = true;
	}
	else if (outOpt.CompareTo("tracks_root", TString::kIgnoreCase) == 0)
	{
		useRootWriter = true;
		makeTracksOnly = true;
	}
	else
	{
		cerr << "ERROR: Unknown option for output: " << output << endl;
		return;
	}
	
	TString chainOpt = chainType;
	if (chainOpt.CompareTo("full", TString::kIgnoreCase) == 0)
	{
		buildDDLFilePubs = true;
		buildDDLRecoComps = true;
		buildTrackerComp = true;
	}
	else if (chainOpt.CompareTo("ddlreco", TString::kIgnoreCase) == 0)
	{
		buildDDLFilePubs = true;
		buildDDLRecoComps = true;
	}
	else if (chainOpt.CompareTo("tracker", TString::kIgnoreCase) == 0)
	{
		buildTrackerComp = true;
		
		TString dataOpt = dataSource;
		if (dataOpt.CompareTo("sim", TString::kIgnoreCase) == 0)
		{
			buildSimDataPubs = true;
		}
		else if (dataOpt.CompareTo("rec", TString::kIgnoreCase) == 0)
		{
			buildRecDataPubs = true;
		}
		else
		{
			cerr << "ERROR: Unknown option for dataSource: " << dataSource << endl;
			return;
		}
	}
	else
	{
		cerr << "ERROR: Unknown option for chainType: " << chainType << endl;
		return;
	}
	
	TString logOpt = logLevel;
	if (logOpt.CompareTo("normal", TString::kIgnoreCase) == 0)
	{
		// nothing to do.
	}
	else if (logOpt.CompareTo("max", TString::kIgnoreCase) == 0)
	{
		maxLogging = true;
	}
	else if (logOpt.CompareTo("min", TString::kIgnoreCase) == 0)
	{
		minLogging = true;
	}
	else
	{
		cerr << "ERROR: Unknown option for logLevel: " << logLevel << endl;
		return;
	}
	
	// If we are supposed to make tracks only but are in a ddlreco chain
	// then we clearly can only generate the DDL reconstructed data, so do that.
	if (makeTracksOnly && ! buildTrackerComp)
	{
		makeTracksOnly = false;
	}
	
	// Now we can initialise the AliHLTSystem...
	AliHLTSystem sys;
	
	// Start by setting up the logging.
	if (maxLogging)
	{
		AliLog::SetGlobalLogLevel(AliLog::kMaxType);
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	if (minLogging)
	{
		sys.SetGlobalLoggingLevel(kHLTLogError);
	}
	
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	sys.LoadComponentLibraries("libAliHLTMUON.so");

	// The DDL file publishers are only needed if we create the ddlreco or
	// full chains. The filename lists are built assuming the aliroot rawXX/
	// directory structure.
	if (buildDDLFilePubs)
	{
		string cmd13 = "-datatype 'DDLTRACK' 'MUON' -dataspec 0x001000";
		string cmd14 = "-datatype 'DDLTRACK' 'MUON' -dataspec 0x002000";
		string cmd15 = "-datatype 'DDLTRACK' 'MUON' -dataspec 0x004000";
		string cmd16 = "-datatype 'DDLTRACK' 'MUON' -dataspec 0x008000";
		string cmd17 = "-datatype 'DDLTRACK' 'MUON' -dataspec 0x010000";
		string cmd18 = "-datatype 'DDLTRACK' 'MUON' -dataspec 0x020000";
		string cmd19 = "-datatype 'DDLTRACK' 'MUON' -dataspec 0x040000";
		string cmd20 = "-datatype 'DDLTRACK' 'MUON' -dataspec 0x080000";
		string cmd21 = "-datatype 'DDLTRIGR' 'MUON' -dataspec 0x100000";
		string cmd22 = "-datatype 'DDLTRIGR' 'MUON' -dataspec 0x200000";
		for (int i = firstEvent; i < lastEvent+1; i++)
		{
			if (i != 0)
			{
				cmd13 += " -nextevent";
				cmd14 += " -nextevent";
				cmd15 += " -nextevent";
				cmd16 += " -nextevent";
				cmd17 += " -nextevent";
				cmd18 += " -nextevent";
				cmd19 += " -nextevent";
				cmd20 += " -nextevent";
				cmd21 += " -nextevent";
				cmd22 += " -nextevent";
			}
			char buf[16];
			sprintf(buf, "%d", i);
			cmd13 += " -datafile raw"; cmd13 += buf; cmd13 += "/MUONTRK_2572.ddl";
			cmd14 += " -datafile raw"; cmd14 += buf; cmd14 += "/MUONTRK_2573.ddl";
			cmd15 += " -datafile raw"; cmd15 += buf; cmd15 += "/MUONTRK_2574.ddl";
			cmd16 += " -datafile raw"; cmd16 += buf; cmd16 += "/MUONTRK_2575.ddl";
			cmd17 += " -datafile raw"; cmd17 += buf; cmd17 += "/MUONTRK_2576.ddl";
			cmd18 += " -datafile raw"; cmd18 += buf; cmd18 += "/MUONTRK_2577.ddl";
			cmd19 += " -datafile raw"; cmd19 += buf; cmd19 += "/MUONTRK_2578.ddl";
			cmd20 += " -datafile raw"; cmd20 += buf; cmd20 += "/MUONTRK_2579.ddl";
			cmd21 += " -datafile raw"; cmd21 += buf; cmd21 += "/MUONTRG_2816.ddl";
			cmd22 += " -datafile raw"; cmd22 += buf; cmd22 += "/MUONTRG_2817.ddl";
		}

		AliHLTConfiguration pubDDL13("pubDDL13", "FilePublisher", NULL, cmd13.c_str());
		AliHLTConfiguration pubDDL14("pubDDL14", "FilePublisher", NULL, cmd14.c_str());
		AliHLTConfiguration pubDDL15("pubDDL15", "FilePublisher", NULL, cmd15.c_str());
		AliHLTConfiguration pubDDL16("pubDDL16", "FilePublisher", NULL, cmd16.c_str());
		AliHLTConfiguration pubDDL17("pubDDL17", "FilePublisher", NULL, cmd17.c_str());
		AliHLTConfiguration pubDDL18("pubDDL18", "FilePublisher", NULL, cmd18.c_str());
		AliHLTConfiguration pubDDL19("pubDDL19", "FilePublisher", NULL, cmd19.c_str());
		AliHLTConfiguration pubDDL20("pubDDL20", "FilePublisher", NULL, cmd20.c_str());
	        AliHLTConfiguration pubDDL21("pubDDL21", "FilePublisher", NULL, cmd21.c_str());
	        AliHLTConfiguration pubDDL22("pubDDL22", "FilePublisher", NULL, cmd22.c_str());
	}
	
	// Build the DDL reconstructor components for all the DDLs 13 to 22, that
	// is for chambers 7 to 10 and trigger stations. We only need to build
	// these components if we are are building the ddlreco or full chains.
	if (buildDDLRecoComps)
	{
		AliHLTConfiguration recDDL13("recDDL13", "MUONHitReconstructor", "pubDDL13", TString("-ddl 13 -buspatchmap ") + lutDir + TString("/BusToDetElem.dat -lut ") + lutDir + TString("/Lut13.dat"));
		AliHLTConfiguration recDDL14("recDDL14", "MUONHitReconstructor", "pubDDL14", TString("-ddl 14 -buspatchmap ") + lutDir + TString("/BusToDetElem.dat -lut ") + lutDir + TString("/Lut14.dat"));
		AliHLTConfiguration recDDL15("recDDL15", "MUONHitReconstructor", "pubDDL15", TString("-ddl 15 -buspatchmap ") + lutDir + TString("/BusToDetElem.dat -lut ") + lutDir + TString("/Lut15.dat"));
		AliHLTConfiguration recDDL16("recDDL16", "MUONHitReconstructor", "pubDDL16", TString("-ddl 16 -buspatchmap ") + lutDir + TString("/BusToDetElem.dat -lut ") + lutDir + TString("/Lut16.dat"));
		AliHLTConfiguration recDDL17("recDDL17", "MUONHitReconstructor", "pubDDL17", TString("-ddl 17 -buspatchmap ") + lutDir + TString("/BusToDetElem.dat -lut ") + lutDir + TString("/Lut17.dat"));
		AliHLTConfiguration recDDL18("recDDL18", "MUONHitReconstructor", "pubDDL18", TString("-ddl 18 -buspatchmap ") + lutDir + TString("/BusToDetElem.dat -lut ") + lutDir + TString("/Lut18.dat"));
		AliHLTConfiguration recDDL19("recDDL19", "MUONHitReconstructor", "pubDDL19", TString("-ddl 19 -buspatchmap ") + lutDir + TString("/BusToDetElem.dat -lut ") + lutDir + TString("/Lut19.dat"));
		AliHLTConfiguration recDDL20("recDDL20", "MUONHitReconstructor", "pubDDL20", TString("-ddl 20 -buspatchmap ") + lutDir + TString("/BusToDetElem.dat -lut ") + lutDir + TString("/Lut20.dat"));	
		AliHLTConfiguration recDDL21("recDDL21", "MUONTriggerReconstructor", "pubDDL21", TString("-ddl 21 -lut ") + lutDir + TString("/Lut21.dat -suppress_partial_triggers"));
		AliHLTConfiguration recDDL22("recDDL22", "MUONTriggerReconstructor", "pubDDL22", TString("-ddl 22 -lut ") + lutDir + TString("/Lut22.dat -suppress_partial_triggers"));
	}

	TString startEventStr = "-firstevent ";
	startEventStr += firstEvent;
        
	// Build the data source components to take data from simulated hits if
	// we are building the tracker only chain with the 'sim' data source.
	if (buildSimDataPubs)
	{
		AliHLTConfiguration recDDL13("recDDL13", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane left  -chamber 7");
		AliHLTConfiguration recDDL14("recDDL14", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane right -chamber 7");
		AliHLTConfiguration recDDL15("recDDL15", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane left  -chamber 8");
		AliHLTConfiguration recDDL16("recDDL16", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane right -chamber 8");
		AliHLTConfiguration recDDL17("recDDL17", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane left  -chamber 9");
		AliHLTConfiguration recDDL18("recDDL18", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane right -chamber 9");
		AliHLTConfiguration recDDL19("recDDL19", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane left  -chamber 10");
		AliHLTConfiguration recDDL20("recDDL20", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane right -chamber 10");
		AliHLTConfiguration recDDL21("recDDL21", "MUONTriggerRecordsSource", NULL, startEventStr + " -hitdata -plane left");
		AliHLTConfiguration recDDL22("recDDL22", "MUONTriggerRecordsSource", NULL, startEventStr + " -hitdata -plane right");
	}
	
        // Build the data source components to take data from offline reconstructed
        // objects if we are building the tracker only chain with the 'rec' data source.
	if (buildRecDataPubs)
	{
		AliHLTConfiguration recDDL13("recDDL13", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane left  -chamber 7");
		AliHLTConfiguration recDDL14("recDDL14", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane right -chamber 7");
		AliHLTConfiguration recDDL15("recDDL15", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane left  -chamber 8");
		AliHLTConfiguration recDDL16("recDDL16", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane right -chamber 8");
		AliHLTConfiguration recDDL17("recDDL17", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane left  -chamber 9");
		AliHLTConfiguration recDDL18("recDDL18", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane right -chamber 9");
		AliHLTConfiguration recDDL19("recDDL19", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane left  -chamber 10");
		AliHLTConfiguration recDDL20("recDDL20", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane right -chamber 10");
		AliHLTConfiguration recDDL21("recDDL21", "MUONTriggerRecordsSource", NULL, startEventStr + " -recdata -plane left");
		AliHLTConfiguration recDDL22("recDDL22", "MUONTriggerRecordsSource", NULL, startEventStr + " -recdata -plane right");
	}
	
	// Build the tracker component if we are building the tracker only or
	// full chains.
	if (buildTrackerComp)
	{
		AliHLTConfiguration tracker("tracker", "MUONMansoTrackerFSM", "recDDL13 recDDL14 recDDL15 recDDL16 recDDL17 recDDL18 recDDL19 recDDL20 recDDL21 recDDL22", "");
	}

	// Build the data sink to subscribe only to what has been created and
	// to the data source we actaully want.
	TString sources = "";
	if (makeTracksOnly)
	{
		sources += "tracker ";
	}
	else
	{
		if (buildTrackerComp)
			sources += "tracker ";
		sources += "recDDL13 recDDL14 recDDL15 recDDL16 recDDL17 recDDL18 recDDL19 recDDL20 recDDL21 recDDL22";
	}
	if (useRootWriter)
	{
		AliHLTConfiguration convert("convert", "MUONRootifier", sources, "");
		AliHLTConfiguration sink("sink", "ROOTFileWriter", "convert", "-datafile output.root -specfmt");
	}
	else
	{
		AliHLTConfiguration sink("sink", "FileWriter", sources, "-datafile output.dat -specfmt");
	}
	
	// Build and run the chain's tasks.
	sys.BuildTaskList("sink");
	sys.Run(eventCount);
}
