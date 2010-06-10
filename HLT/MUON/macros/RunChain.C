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
 * \ingroup macros
 * \file RunChain.C
 * \brief Macro for running the dHLT chain in a standalone mode.
 *
 * This macro is used to run the dHLT component chain within the AliHLTSystem
 * framework, which is a simulation of the dHLT online system's output.
 * To run this macro you must be in the same directory as the galice.root
 * files or in the directory containing the rawXX/ directories produced from
 * AliRoot simulations.
 * Also make sure you specify CDB as for the lutDir (lookup table directory)
 * or that the appropriate LUTs are in the same directory as your working
 * directory, or you can specify the directory with the lutDir parameter option.
 *
 * The simplest way to run this macro with defaults is to copy "rootlogon.C" from
 * $ALICE_ROOT/HLT/MUON/macros into your current working directory, then from
 * the shell command prompt run the following command:
 * \code
 *   > aliroot -b -q -l $ALICE_ROOT/HLT/MUON/macros/RunChain.C+
 * \endcode
 *
 * \author Artur Szostak <artursz@iafrica.com>
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliRawReader.h"
#include "AliHLTOfflineInterface.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliHLTSystem.h"
#include "AliHLTConfiguration.h"
#include "AliLog.h"
#include "TString.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

/**
 * Used to run the dimuon HLT (dHLT) chain in various configuration in a standalone
 * configuration. This is normally used for debugging, testing and can also be used
 * as an example of how to build chains for dHLT by hand.
 *
 * @param chainType  Specifies the type of chain to run. This can be one of the
 *     following:
 *       "full" - Run the full dHLT chain with manso tracker. (default)
 *       "full-tracker" - Run the full dHLT chain with the full tracker.
 *       "ddlreco" - Run only the reconstruction of the DDL raw data up to hits
 *                   and trigger records.
 *       "tracker" - Run the Manso tracker only using hits and trigger records from
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
 *     records when we are running the tracker only chain and the DDL raw data
 *     reading method when running the 'full' or 'ddlreco' chain. If the
 *     chainType = "full" or "ddlreco", then the possible options are:
 *       "file"      - Reads the raw data directly from the DDL files. (default)
 *       "rawreader" - Reads the raw data using the AliRawReader interface from
 *               the current working directory.
 *     If the chainType = "tracker", then the possible options are:
 *       "sim" - Take data from GEANT hits. (default)
 *       "rec" - Take data from offline reconstructed hits and local trigger
 *               objects.
 * @param logLevel  Specifies the amount of logging messages produced by the HLT
 *     system. The possible options are:
 *       "normal" - Shows warnings, errors and information.
 *       "debug" - Shows all messages up to the debug level only for HLT and the
 *                 same as normal for all other modules.
 *       "max" - Shows everything including debug messages if they were compiled in.
 *       "min" - Shows only error messages.
 * @param lutDir  This is the directory in which the LUTs can be found.
 *      If it is set to "CDB" (case sensitive) then the LUTs will be loaded from
 *      CDB instead. The default behaviour is to read from the local CDB store.
 * @param checkData  A flag for indicating if the event data should be checked
 *      for consistency with the AliHLTMUONDataCheckerComponent.
 * @param rawDataPath  The path of the raw data (i.e. path to the rawXX directories)
 *      or can be the file name if using the "rawreader" option for dataSource.
 * @param runNumber  Specifies the run number to use. If it is set to -1 then the
 *      run number is not set if the CDB manager already has a run number set,
 *      otherwise a default run number of 0 is used. The default value is -1.
 * @param cdbPath  This gives the CDB path to use. If it is set to NULL then
 *      the CDB path is not set if the CDB manager already has a default storage
 *      CDB path set, otherwise a default value of "local://$ALICE_ROOT/OCDB" is used.
 *      The default value is NULL.
 * @param tryrecover  If this is true then the "-tryrecover" flag is set in the
 *      raw data reconstruction components. This is useful if when running RunChain
 *      log messages appear indicating that there was a problem decoding the raw data.
 *      The "-tryrecover" flag will turn on recovery logic in the raw data decoders
 *      to try and overcome errors in the data.
 */
void RunChain(
		const char* chainType = "full",
		Int_t firstEvent = 0,
		Int_t lastEvent = -1,
		const char* output = "bin",
		const char* dataSource = "file",
		const char* logLevel = "normal",
		const char* lutDir = "CDB",
		bool checkData = false,
		const char* rawDataPath = "./",
		Int_t runNumber = -1,
		const char* cdbPath = NULL,
		bool tryrecover = false
	)
{
	// Setup the CDB default storage and run number if nothing was set.
	AliCDBManager* cdbManager = AliCDBManager::Instance();
	if (cdbManager == NULL)
	{
		cerr << "ERROR: Global CDB manager object does not exist." << endl;
		return;
	}
	if (runNumber != -1)
	{
		cdbManager->SetRun(runNumber);
	}
	else if (cdbManager->GetRun() == -1)
	{
		cdbManager->SetRun(0);
	}
	if (cdbPath != NULL)
	{
		cdbManager->SetDefaultStorage(cdbPath);
	}
	else if (cdbManager->GetDefaultStorage() == NULL)
	{
		cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	}
	
	if (cdbManager->GetDefaultStorage() == NULL)
	{
		cerr << "ERROR: There is no value for the default CDB storage, cannot continue." << endl;
		return;
	}

	// Make sure that the lastEvent is greater than firstEvent.
	if (lastEvent < firstEvent)
		lastEvent = firstEvent;
	int eventCount = lastEvent - firstEvent + 1;
	
	bool buildDDLFilePubs = false;
	bool buildRawReaderPubs = false;
	bool buildDDLRecoComps = false;
	bool buildSimDataPubs = false;
	bool buildRecDataPubs = false;
	bool buildTrackerComp = false;
	bool buildFullTrackerComp = false;
	bool maxLogging = false;
	bool debugLogging = false;
	bool minLogging = false;
	bool useRootWriter = false;
	bool makeTracksOnly = false;
	bool buildDecisionComp = false;
	bool buildESDComp = false;
	
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
		cerr << "ERROR: Unknown option for output: '" << output << "'" << endl;
		return;
	}
	
	TString chainOpt = chainType;
	if (chainOpt.CompareTo("full", TString::kIgnoreCase) == 0)
	{
		buildDDLRecoComps = true;
		buildTrackerComp = true;
		buildDecisionComp = true;
		buildESDComp = true;
		
		TString dataOpt = dataSource;
		if (dataOpt.CompareTo("file", TString::kIgnoreCase) == 0)
		{
			buildDDLFilePubs = true;
		}
		else if (dataOpt.CompareTo("rawreader", TString::kIgnoreCase) == 0)
		{
			buildRawReaderPubs = true;
		}
		else
		{
			cerr << "ERROR: Unknown option for dataSource: '" << dataSource
				<< "'. Valid options are: 'file' or 'rawreader'" << endl;
			return;
		}
	}
	else if (chainOpt.CompareTo("full-tracker", TString::kIgnoreCase) == 0)
	{
		buildDDLRecoComps = true;
		buildFullTrackerComp = true;
		buildDecisionComp = true;
		buildESDComp = true;
		
		TString dataOpt = dataSource;
		if (dataOpt.CompareTo("file", TString::kIgnoreCase) == 0)
		{
			buildDDLFilePubs = true;
		}
		else if (dataOpt.CompareTo("rawreader", TString::kIgnoreCase) == 0)
		{
			buildRawReaderPubs = true;
		}
		else
		{
			cerr << "ERROR: Unknown option for dataSource: '" << dataSource
				<< "'. Valid options are: 'file' or 'rawreader'" << endl;
			return;
		}
	}
	else if (chainOpt.CompareTo("ddlreco", TString::kIgnoreCase) == 0)
	{
		buildDDLRecoComps = true;
		buildDecisionComp = false;
		buildESDComp = false;
		
		TString dataOpt = dataSource;
		if (dataOpt.CompareTo("file", TString::kIgnoreCase) == 0)
		{
			buildDDLFilePubs = true;
		}
		else if (dataOpt.CompareTo("rawreader", TString::kIgnoreCase) == 0)
		{
			buildRawReaderPubs = true;
		}
		else
		{
			cerr << "ERROR: Unknown option for dataSource: '" << dataSource
				<< "'. Valid options are: 'file' or 'rawreader'" << endl;
			return;
		}
	}
	else if (chainOpt.CompareTo("tracker", TString::kIgnoreCase) == 0)
	{
		buildTrackerComp = true;
		buildDecisionComp = true;
		buildESDComp = true;
		
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
			cerr << "ERROR: Unknown option for dataSource: '" << dataSource
				<< "'. Valid options are: 'sim' or 'rec'" << endl;
			return;
		}
	}
	else
	{
		cerr << "ERROR: Unknown option for chainType: '" << chainType << "'" << endl;
		return;
	}
	
	TString logOpt = logLevel;
	if (logOpt.CompareTo("normal", TString::kIgnoreCase) == 0)
	{
		// nothing to do.
	}
	else if (logOpt.CompareTo("debug", TString::kIgnoreCase) == 0)
	{
		debugLogging = true;
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
		cerr << "ERROR: Unknown option for logLevel: '" << logLevel << "'" << endl;
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
	if (debugLogging)
	{
		AliLog::SetModuleDebugLevel("HLT", AliLog::kMaxType);
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	if (minLogging)
	{
		sys.SetGlobalLoggingLevel(AliHLTComponentLogSeverity(
			kHLTLogFatal | kHLTLogError
		));
	}
	
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	sys.LoadComponentLibraries("libAliHLTMUON.so");

	// The DDL file publishers are only needed if we create the ddlreco or
	// full chains. The filename lists are built assuming the aliroot rawXX/
	// directory structure.
	if (buildDDLFilePubs)
	{
		string cmd1 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000001";
		string cmd2 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000002";
		string cmd3 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000004";
		string cmd4 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000008";
		string cmd5 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000010";
		string cmd6 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000020";
		string cmd7 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000040";
		string cmd8 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000080";
		string cmd9 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000100";
		string cmd10 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000200";
		string cmd11 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000400";
		string cmd12 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x000800";
		string cmd13 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x001000";
		string cmd14 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x002000";
		string cmd15 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x004000";
		string cmd16 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x008000";
		string cmd17 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x010000";
		string cmd18 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x020000";
		string cmd19 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x040000";
		string cmd20 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x080000";
		string cmd21 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x100000";
		string cmd22 = "-datatype 'DDL_RAW ' 'MUON' -dataspec 0x200000";
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
			cmd1 += " -datafile "; cmd1 += rawDataPath; cmd1 += "raw"; cmd1 += buf; cmd1 += "/MUONTRK_2560.ddl";
			cmd2 += " -datafile "; cmd2 += rawDataPath; cmd2 += "raw"; cmd2 += buf; cmd2 += "/MUONTRK_2561.ddl";
			cmd3 += " -datafile "; cmd3 += rawDataPath; cmd3 += "raw"; cmd3 += buf; cmd3 += "/MUONTRK_2562.ddl";
			cmd4 += " -datafile "; cmd4 += rawDataPath; cmd4 += "raw"; cmd4 += buf; cmd4 += "/MUONTRK_2563.ddl";
			cmd5 += " -datafile "; cmd5 += rawDataPath; cmd5 += "raw"; cmd5 += buf; cmd5 += "/MUONTRK_2564.ddl";
			cmd6 += " -datafile "; cmd6 += rawDataPath; cmd6 += "raw"; cmd6 += buf; cmd6 += "/MUONTRK_2565.ddl";
			cmd7 += " -datafile "; cmd7 += rawDataPath; cmd7 += "raw"; cmd7 += buf; cmd7 += "/MUONTRK_2566.ddl";
			cmd8 += " -datafile "; cmd8 += rawDataPath; cmd8 += "raw"; cmd8 += buf; cmd8 += "/MUONTRK_2567.ddl";
			cmd9 += " -datafile "; cmd9 += rawDataPath; cmd9 += "raw"; cmd9 += buf; cmd9 += "/MUONTRK_2568.ddl";
			cmd10 += " -datafile "; cmd10 += rawDataPath; cmd10 += "raw"; cmd10 += buf; cmd10 += "/MUONTRK_2569.ddl";
			cmd11 += " -datafile "; cmd11 += rawDataPath; cmd11 += "raw"; cmd11 += buf; cmd11 += "/MUONTRK_2570.ddl";
			cmd12 += " -datafile "; cmd12 += rawDataPath; cmd12 += "raw"; cmd12 += buf; cmd12 += "/MUONTRK_2571.ddl";
			cmd13 += " -datafile "; cmd13 += rawDataPath; cmd13 += "raw"; cmd13 += buf; cmd13 += "/MUONTRK_2572.ddl";
			cmd14 += " -datafile "; cmd14 += rawDataPath; cmd14 += "raw"; cmd14 += buf; cmd14 += "/MUONTRK_2573.ddl";
			cmd15 += " -datafile "; cmd15 += rawDataPath; cmd15 += "raw"; cmd15 += buf; cmd15 += "/MUONTRK_2574.ddl";
			cmd16 += " -datafile "; cmd16 += rawDataPath; cmd16 += "raw"; cmd16 += buf; cmd16 += "/MUONTRK_2575.ddl";
			cmd17 += " -datafile "; cmd17 += rawDataPath; cmd17 += "raw"; cmd17 += buf; cmd17 += "/MUONTRK_2576.ddl";
			cmd18 += " -datafile "; cmd18 += rawDataPath; cmd18 += "raw"; cmd18 += buf; cmd18 += "/MUONTRK_2577.ddl";
			cmd19 += " -datafile "; cmd19 += rawDataPath; cmd19 += "raw"; cmd19 += buf; cmd19 += "/MUONTRK_2578.ddl";
			cmd20 += " -datafile "; cmd20 += rawDataPath; cmd20 += "raw"; cmd20 += buf; cmd20 += "/MUONTRK_2579.ddl";
			cmd21 += " -datafile "; cmd21 += rawDataPath; cmd21 += "raw"; cmd21 += buf; cmd21 += "/MUONTRG_2816.ddl";
			cmd22 += " -datafile "; cmd22 += rawDataPath; cmd22 += "raw"; cmd22 += buf; cmd22 += "/MUONTRG_2817.ddl";
		}
		
		AliHLTConfiguration pubDDL1("pubDDL1", "FilePublisher", NULL, cmd1.c_str());
		AliHLTConfiguration pubDDL2("pubDDL2", "FilePublisher", NULL, cmd2.c_str());
		AliHLTConfiguration pubDDL3("pubDDL3", "FilePublisher", NULL, cmd3.c_str());
		AliHLTConfiguration pubDDL4("pubDDL4", "FilePublisher", NULL, cmd4.c_str());
		AliHLTConfiguration pubDDL5("pubDDL5", "FilePublisher", NULL, cmd5.c_str());
		AliHLTConfiguration pubDDL6("pubDDL6", "FilePublisher", NULL, cmd6.c_str());
		AliHLTConfiguration pubDDL7("pubDDL7", "FilePublisher", NULL, cmd7.c_str());
		AliHLTConfiguration pubDDL8("pubDDL8", "FilePublisher", NULL, cmd8.c_str());
		AliHLTConfiguration pubDDL9("pubDDL9", "FilePublisher", NULL, cmd9.c_str());
		AliHLTConfiguration pubDDL10("pubDDL10", "FilePublisher", NULL, cmd10.c_str());
		AliHLTConfiguration pubDDL11("pubDDL11", "FilePublisher", NULL, cmd11.c_str());
		AliHLTConfiguration pubDDL12("pubDDL12", "FilePublisher", NULL, cmd12.c_str());
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

	// Build the DDL file publishers using AliRawReaderPublisher components.
	if (buildRawReaderPubs)
	{
		string cmd1 = "-skipempty -minid 2560 -maxid 2560 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000001";
		string cmd2 = "-skipempty -minid 2561 -maxid 2561 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000002";
		string cmd3 = "-skipempty -minid 2562 -maxid 2562 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000004";
		string cmd4 = "-skipempty -minid 2563 -maxid 2563 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000008";
		string cmd5 = "-skipempty -minid 2564 -maxid 2564 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000010";
		string cmd6 = "-skipempty -minid 2565 -maxid 2565 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000020";
		string cmd7 = "-skipempty -minid 2566 -maxid 2566 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000040";
		string cmd8 = "-skipempty -minid 2567 -maxid 2567 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000080";
		string cmd9 = "-skipempty -minid 2568 -maxid 2568 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000100";
		string cmd10 = "-skipempty -minid 2569 -maxid 2569 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000200";
		string cmd11 = "-skipempty -minid 2570 -maxid 2570 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000400";
		string cmd12 = "-skipempty -minid 2571 -maxid 2571 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000800";
		string cmd13 = "-skipempty -minid 2572 -maxid 2572 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x001000";
		string cmd14 = "-skipempty -minid 2573 -maxid 2573 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x002000";
		string cmd15 = "-skipempty -minid 2574 -maxid 2574 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x004000";
		string cmd16 = "-skipempty -minid 2575 -maxid 2575 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x008000";
		string cmd17 = "-skipempty -minid 2576 -maxid 2576 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x010000";
		string cmd18 = "-skipempty -minid 2577 -maxid 2577 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x020000";
		string cmd19 = "-skipempty -minid 2578 -maxid 2578 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x040000";
		string cmd20 = "-skipempty -minid 2579 -maxid 2579 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x080000";
		string cmd21 = "-skipempty -minid 2816 -maxid 2816 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x100000";
		string cmd22 = "-skipempty -minid 2817 -maxid 2817 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x200000";

		AliHLTConfiguration pubDDL1("pubDDL1", "AliRawReaderPublisher", NULL, cmd1.c_str());
		AliHLTConfiguration pubDDL2("pubDDL2", "AliRawReaderPublisher", NULL, cmd2.c_str());
		AliHLTConfiguration pubDDL3("pubDDL3", "AliRawReaderPublisher", NULL, cmd3.c_str());
		AliHLTConfiguration pubDDL4("pubDDL4", "AliRawReaderPublisher", NULL, cmd4.c_str());
		AliHLTConfiguration pubDDL5("pubDDL5", "AliRawReaderPublisher", NULL, cmd5.c_str());
		AliHLTConfiguration pubDDL6("pubDDL6", "AliRawReaderPublisher", NULL, cmd6.c_str());
		AliHLTConfiguration pubDDL7("pubDDL7", "AliRawReaderPublisher", NULL, cmd7.c_str());
		AliHLTConfiguration pubDDL8("pubDDL8", "AliRawReaderPublisher", NULL, cmd8.c_str());
		AliHLTConfiguration pubDDL9("pubDDL9", "AliRawReaderPublisher", NULL, cmd9.c_str());
		AliHLTConfiguration pubDDL10("pubDDL10", "AliRawReaderPublisher", NULL, cmd10.c_str());
		AliHLTConfiguration pubDDL11("pubDDL11", "AliRawReaderPublisher", NULL, cmd11.c_str());
		AliHLTConfiguration pubDDL12("pubDDL12", "AliRawReaderPublisher", NULL, cmd12.c_str());
		AliHLTConfiguration pubDDL13("pubDDL13", "AliRawReaderPublisher", NULL, cmd13.c_str());
		AliHLTConfiguration pubDDL14("pubDDL14", "AliRawReaderPublisher", NULL, cmd14.c_str());
		AliHLTConfiguration pubDDL15("pubDDL15", "AliRawReaderPublisher", NULL, cmd15.c_str());
		AliHLTConfiguration pubDDL16("pubDDL16", "AliRawReaderPublisher", NULL, cmd16.c_str());
		AliHLTConfiguration pubDDL17("pubDDL17", "AliRawReaderPublisher", NULL, cmd17.c_str());
		AliHLTConfiguration pubDDL18("pubDDL18", "AliRawReaderPublisher", NULL, cmd18.c_str());
		AliHLTConfiguration pubDDL19("pubDDL19", "AliRawReaderPublisher", NULL, cmd19.c_str());
		AliHLTConfiguration pubDDL20("pubDDL20", "AliRawReaderPublisher", NULL, cmd20.c_str());
	        AliHLTConfiguration pubDDL21("pubDDL21", "AliRawReaderPublisher", NULL, cmd21.c_str());
	        AliHLTConfiguration pubDDL22("pubDDL22", "AliRawReaderPublisher", NULL, cmd22.c_str());
	}
	
	// Build the DDL reconstructor components for all the DDLs 13 to 22, that
	// is for chambers 7 to 10 and trigger stations. We only need to build
	// these components if we are are building the ddlreco or full chains.
	if (buildDDLRecoComps)
	{
		const char* recoverFlag = tryrecover ? "-tryrecover" : "";
		for (int k = 1; k <= 22; k++)
		{
			string compId = Form("recDDL%d", k);
			string name = (k <= 20) ? "MUONHitReconstructor" : "MUONTriggerReconstructor";
			string parent = Form("pubDDL%d", k);
			string cmd;
			const char* extraInfoFlags = k < 21 ? "-makeclusters -makechannels" : "-makedebuginfo";
			if (TString(lutDir) == "CDB")
			{
				const char* path = cdbManager->GetDefaultStorage()->GetURI().Data();
				cmd = Form("-ddl %d -cdbpath %s -run %d %s %s",
					k, path, cdbManager->GetRun(), recoverFlag, extraInfoFlags
				);
			}
			else
			{
				cmd = Form("-ddl %d -lut %s/Lut%d.dat %s %s",
					k, lutDir, k, recoverFlag, extraInfoFlags
				);
			}
			if (k >= 21)
			{
				cmd += " -suppress_partial_triggers -dont_use_crateid -dont_use_localid";
			}
			AliHLTConfiguration recDDL(compId.c_str(), name.c_str(), parent.c_str(), cmd.c_str());
		}
	}

	TString startEventStr = "-firstevent ";
	startEventStr += firstEvent;
	
	// Build the data source components to take data from simulated hits if
	// we are building the tracker only chain with the 'sim' data source.
	if (buildSimDataPubs)
	{
		AliHLTConfiguration recDDL1("recDDL1", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane left  -chamber 1");
		AliHLTConfiguration recDDL2("recDDL2", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane right -chamber 1");
		AliHLTConfiguration recDDL3("recDDL3", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane left  -chamber 2");
		AliHLTConfiguration recDDL4("recDDL4", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane right -chamber 2");
		AliHLTConfiguration recDDL5("recDDL5", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane left  -chamber 3");
		AliHLTConfiguration recDDL6("recDDL6", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane right -chamber 3");
		AliHLTConfiguration recDDL7("recDDL7", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane left  -chamber 4");
		AliHLTConfiguration recDDL8("recDDL8", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane right -chamber 4");
		AliHLTConfiguration recDDL9("recDDL9", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane left  -chamber 5");
		AliHLTConfiguration recDDL10("recDDL10", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane right -chamber 5");
		AliHLTConfiguration recDDL11("recDDL11", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane left  -chamber 6");
		AliHLTConfiguration recDDL12("recDDL12", "MUONRecHitsSource", NULL, startEventStr + " -simdata -plane right -chamber 6");
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
		AliHLTConfiguration recDDL1("recDDL1", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane left  -chamber 1");
		AliHLTConfiguration recDDL2("recDDL2", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane right -chamber 1");
		AliHLTConfiguration recDDL3("recDDL3", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane left  -chamber 2");
		AliHLTConfiguration recDDL4("recDDL4", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane right -chamber 2");
		AliHLTConfiguration recDDL5("recDDL5", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane left  -chamber 3");
		AliHLTConfiguration recDDL6("recDDL6", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane right -chamber 3");
		AliHLTConfiguration recDDL7("recDDL7", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane left  -chamber 4");
		AliHLTConfiguration recDDL8("recDDL8", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane right -chamber 4");
		AliHLTConfiguration recDDL9("recDDL9", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane left  -chamber 5");
		AliHLTConfiguration recDDL10("recDDL10", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane right -chamber 5");
		AliHLTConfiguration recDDL11("recDDL11", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane left  -chamber 6");
		AliHLTConfiguration recDDL12("recDDL12", "MUONRecHitsSource", NULL, startEventStr + " -recdata -plane right -chamber 6");
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
		AliHLTConfiguration tracker(
			"tracker",
			"MUONMansoTrackerFSM",
			"recDDL13 recDDL14 recDDL15 recDDL16 recDDL17 recDDL18 recDDL19 recDDL20 recDDL21 recDDL22",
			"-makecandidates"
		);
	}
	if (buildFullTrackerComp)
	{
		AliHLTConfiguration fulltracker(
			"tracker-full",
			"MUONFullTracker",
			"recDDL1 recDDL2 recDDL3 recDDL4 recDDL5 recDDL6 recDDL7 recDDL8 recDDL9 recDDL10 recDDL11"
			" recDDL12 recDDL13 recDDL14 recDDL15 recDDL16 recDDL17 recDDL18 recDDL19 recDDL20 recDDL21 recDDL22",
			""
		);
	}
	
	// Build the dHLT trigger decision component if enabled.
	if (buildDecisionComp)
	{
		const char* decisionSource = "tracker";
		if (buildFullTrackerComp) decisionSource = "tracker-full";
		AliHLTConfiguration decision("decision", "MUONDecisionComponent", decisionSource, "");
	}
	if (buildESDComp)
	{
	        const char* ESDSource = "tracker recDDL21 recDDL22 ";
		if (buildFullTrackerComp) ESDSource = "tracker-full recDDL21 recDDL22 ";
		AliHLTConfiguration ESD("ESD", AliHLTMUONConstants::ESDMakerId(), ESDSource, " -make_minimal_esd");
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
		if (buildTrackerComp) sources += "tracker ";
		if (buildFullTrackerComp) sources += "tracker-full ";
		sources += "recDDL1 recDDL2 recDDL3 recDDL4 recDDL5 recDDL6 recDDL7 recDDL8 recDDL9 recDDL10 recDDL11"
			" recDDL12 recDDL13 recDDL14 recDDL15 recDDL16 recDDL17 recDDL18 recDDL19 recDDL20 recDDL21 recDDL22";
	}
	if (buildDecisionComp)
	{
		// Add the trigger decision component if it was enabled.
		sources += " decision";
	}
	if (buildESDComp)
	{
		sources += " ESD ";
	}
	
	// Build the data checker component if so requested.
	if (checkData)
	{
		AliHLTConfiguration checker("checker", "MUONDataChecker", sources, "-warn_on_unexpected_block");
		sources = "checker"; // Only add the checker to the sources in this case since it will forward just the bad events.
	}
	
	if (useRootWriter)
	{
		AliHLTConfiguration convert("convert", "MUONRootifier", sources, "");
		AliHLTConfiguration sink("sink", "ROOTFileWriter", "convert", "-concatenate-events -datafile output.root -specfmt");
	}
	else
	{
		AliHLTConfiguration sink("sink", "FileWriter", sources, "-datafile output.dat -specfmt");
	}
	
	// Build and run the chain's tasks.
	sys.BuildTaskList("sink");
	if (buildDDLFilePubs)
	{
		sys.Run(eventCount);
	}
	else
	{
		// Setup the raw reader.
		AliRawReader* rawReader = AliRawReader::Create(rawDataPath);
		if (rawReader == NULL)
		{
			cerr << "ERROR: Could not create raw reader." << endl;
			return;
		}
		if (! rawReader->IsRawReaderValid())
		{
			cerr << "ERROR: Raw reader is not valid." << endl;
			return;
		}
		AliHLTOfflineInterface::SetParamsToComponents(NULL, rawReader);
		rawReader->NextEvent(); // Need to call this once here or we will start at the wrong event.
		// Now step through the events.
		for (int i = 0; i < firstEvent; i++) rawReader->NextEvent();
		for (int i = firstEvent; i <= lastEvent; i++)
		{
			// The "(i == lastEvent) ? 1 : 0" part is to indicate a
			// stop run command when we hit the last event.
			sys.Run(1, (i == lastEvent) ? 1 : 0);
			rawReader->NextEvent();
		}
	}
}
