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

/* $Id: $ */

/**
 * @file   dHLTrootify.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   14 May 2008
 * @brief  Command line utility to convert dHLT's internal raw data blocks into ROOT objects.
 */

#include "TClassTable.h"
#include "TString.h"
#include "AliHLTSystem.h"
#include "AliHLTConfiguration.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliLog.h"

#include <cstring>
#include <cstdlib>
#include <cassert>
#include <new>
#include <fstream>

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;


#define CMDLINE_ERROR 1
#define PARSE_ERROR 2
#define SYSTEM_ERROR 3
#define FATAL_ERROR 4
#define HLTSYSTEM_ERROR 5

namespace
{
	// CDB path and run number to use.
	const char* gCDBPath = "local://$ALICE_ROOT/OCDB";
	Int_t gRunNumber = 0;
}

/**
 * Uses AliHLTSystem and the AliHLTMUONRootifierComponent to convert the files
 * into ROOT object format.
 * @param filenames  Array of file name strings.
 * @param filetypes  Array of file types corresponding to each filename.
 * @param numOfFiles  Number of entries in the 'filenames' and 'filetypes' arrays.
 * @param outputFile  The output file name to use.
 * @param maxLogging  If set then all debug messages are printed. (default = false)
 * @param checkData  Flag indicating if the internal raw data should be checked for
 *          consistency. Errors and warnings will be generated for consistency errors,
 *          but the rootification of the data will continue.
 * @return  Returns HLTSYSTEM_ERROR if there was a problem reported by AliHLTSystem
 *          and EXIT_SUCCESS if the ROOT file was created OK. SYSTEM_ERROR is
 *          returned if there was a problem allocating memory for HLT configuration
 *          objects.
 */
int RootifyFiles(
		const char** filenames,
		AliHLTMUONDataBlockType* filetypes,
		int numOfFiles,
		const char* outputFile,
		bool maxLogging = false,
		bool checkData = false
	)
{
	AliHLTSystem sys;
	
	if (maxLogging)
	{
		AliLog::SetGlobalLogLevel(AliLog::kMaxType);
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	sys.LoadComponentLibraries("libAliHLTMUON.so");
	
	TString sources = "";
	typedef AliHLTConfiguration* PAliHLTConfiguration;
	PAliHLTConfiguration* filePubs = NULL;
	
	try
	{
		filePubs = new PAliHLTConfiguration[numOfFiles];
		// Must make sure all the pointers are NULL because we
		// need to clean up afterwords and we might fail half
		// way through the memory allocation.
		int i;
		for (i = 0; i < numOfFiles; i++)
		{
			filePubs[i] = NULL;
		}
		
		// Now start allocating the file publishers.
		for (i = 0; i < numOfFiles; i++)
		{
			TString name = "filePublisher_";
			name += filenames[i];
			sources += name + " ";
			TString params = "-datatype '";
			params += AliHLTMUONUtils::DataBlockTypeToString(filetypes[i]);
			params += "' 'MUON' -dataspec 0x0 -datafile ";
			params += filenames[i];
			filePubs[i] = new AliHLTConfiguration(
					name.Data(), "FilePublisher", NULL, params.Data()
				);
		}
	}
	catch (const std::bad_alloc&)
	{
		cerr << "ERROR: There is not enough memory to allocate another configuration object." << endl;
		
		// Make sure to clean up what was actaully allocated.
		if (filePubs != NULL)
		{
			for (int i = 0; i < numOfFiles; i++)
			{
				if (filePubs[i] != NULL)
					delete filePubs[i];
			}
			delete [] filePubs;
		}
		
		return SYSTEM_ERROR;
	}
	
	// Setup the component for data integrity checking.
	if (checkData)
	{
		TString dcparams = "-warn_on_unexpected_block -ignorespec -cdbpath ";
		dcparams += gCDBPath;
		dcparams += " -run ";
		dcparams += gRunNumber;
		AliHLTConfiguration checker(
				"checker", AliHLTMUONConstants::DataCheckerComponentId(),
				sources, dcparams
			);
		sources = "checker";
	}
	
	// Setup the component which converts the raw internal dHLT data to ROOT objects.
	AliHLTConfiguration convert("convert", AliHLTMUONConstants::RootifierComponentId(), sources, "");
	
	// Setup the ROOT file writer.
	TString params = "-concatenate-events -datafile ";
	params += outputFile;
	params += " -specfmt";
	AliHLTConfiguration sink("sink", "ROOTFileWriter", "convert", params.Data());
	
	// Build and run the HLT tasks.
	if (sys.BuildTaskList("sink") != 0) return HLTSYSTEM_ERROR;
	if (maxLogging) sys.PrintTaskList();
	if (sys.Run() != 0) return HLTSYSTEM_ERROR;
	
	// Clean up all the dynamically allocate objects.
	for (int i = 0; i < numOfFiles; i++)
	{
		delete filePubs[i];
	}
	delete [] filePubs;

	return EXIT_SUCCESS;
}

/**
 * This method decodes the type of the file.
 */
int DecodeFileType(const char* filename, AliHLTMUONDataBlockType& type)
{
	assert( filename != NULL );
	
	// Open the file and find its size.
	fstream file;
	file.open(filename, ios::in);
	if (not file)
	{
		cerr << "ERROR: Could not open the file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	
	AliHLTMUONDataBlockHeader header;
	file.read(reinterpret_cast<char*>(&header), sizeof(header));
	if (not file)
	{
		cerr << "ERROR: Could not read from file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	
	file.close();
	if (not file)
	{
		cerr << "ERROR: Could not close the file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	
	type = AliHLTMUONDataBlockType(header.fType);
	
	// Check that the type is something valid. Otherwise set it to a
	// value of kUnknownDataBlock.
	switch (type)
	{
	case kTriggerRecordsDataBlock:
	case kTrigRecsDebugDataBlock:
	case kRecHitsDataBlock:
	case kClustersDataBlock:
	case kChannelsDataBlock:
	case kMansoTracksDataBlock:
	case kMansoCandidatesDataBlock:
	case kSinglesDecisionDataBlock:
	case kPairsDecisionDataBlock:
		break;
	default: type = kUnknownDataBlock;
		break;
	}
	
	return EXIT_SUCCESS;
}

/**
 * Prints the command line usage of this program to standard out or error.
 */
void PrintUsage(bool asError = true)
{
	std::ostream& os = asError ? cerr : cout;
	os << "Usage: dHLTrootify [-help|-h] [-outfile|-o <output_file>] [-type|-t <typename>]" << endl;
	os << "         [-debug|-d] [-check|-c] [-cdbpath|-p <url>] [-run|-r <number>]" << endl;
	os << "         <filename> [<filename> ...]" << endl;
	os << "Where <filename> is the name of a file containing a raw data block." << endl;
	os << "Options:" << endl;
	os << " -help | -h" << endl;
	os << "       Displays this message." << endl;
	os << " -outfile | -o <output_file>" << endl;
	os << "       Specifies the output ROOT file to write to, where <output_file> is the" << endl;
	os << "       name of the output file." << endl;
	os << " -type | -t <typename>" << endl;
	os << "       Forces the contents of the subsequent files specified on the command" << endl;
	os << "       line to be interpreted as a specific type of data block." << endl;
	os << "       Where <typename> can be one of:" << endl;
	os << "         trigrecs - trigger records data." << endl;
	os << "         trigrecsdebug - debugging information about trigger records." << endl;
	os << "         rechits - reconstructed hits data." << endl;
	os << "         channels - channel debugging information from hit reconstruction." << endl;
	os << "         clusters - cluster debugging information from hit reconstruction." << endl;
	os << "         mansotracks - partial tracks from Manso algorithm." << endl;
	os << "         mansocandidates - track candidates considered in the Manso algorithm." << endl;
	os << "         singlesdecision - trigger decisions for single tracks." << endl;
	os << "         pairsdecision - trigger decisions for track pairs." << endl;
	os << "         autodetect - the type of the data block will be automatically" << endl;
	os << "                      detected." << endl;
	os << " -debug | -d" << endl;
	os << "       If specified then the all debug messages are printed by the AliHLTSystem." << endl;
	os << " -check | -c" << endl;
	os << "       If specified then data integrity checks are performed on the raw data." << endl;
	os << "       Warnings and errors are printed as problems are found with the data, but" << endl;
	os << "       the data will still be converted into ROOT objects as best as possible." << endl;
	os << " -cdbpath | -p <url>" << endl;
	os << "       The path to the CDB to use when running with the -check | -k option." << endl;
	os << " -run | -r <number>" << endl;
	os << "       The run number to use when running with the -check | -k option." << endl;
}

/**
 * Parses the command line.
 * @param argc  Number of arguments as given in main().
 * @param argv  Array of arguments as given in main().
 * @param filenames  Pointer to buffer storing file name strings.
 * @param filetypes  Array that receives the type of the data block expected, i.e.
 *                   the value of the -type flag for the corresponding file.
 * @param numOfFiles  Receives the number of file name strings that were found
 *                    and added to 'filenames'.
 * @param outputFile  The output file name to use. Set to NULL if none specified.
 * @param maxLogging  Set to true if maximal logging was requested.
 * @param checkData  Set to true if data integrity checking was requested.
 * @return  A status flag suitable for returning from main(), containing either
 *          EXIT_SUCCESS or CMDLINE_ERROR or SYSTEM_ERROR.
 */
int ParseCommandLine(
		int argc,
		const char** argv,
		const char** filenames,
		AliHLTMUONDataBlockType* filetypes,
		int& numOfFiles,
		const char*& outputFile,
		bool& maxLogging,
		bool& checkData
	)
{
	numOfFiles = 0;
	outputFile = NULL;
	maxLogging = false;
	checkData = false;
	bool pathSet = false;
	bool runSet = false;
	AliHLTMUONDataBlockType currentType = kUnknownDataBlock;

	// Parse the command line.
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-help") == 0 or strcmp(argv[i], "-h") == 0)
		{
			PrintUsage(false);
			return EXIT_SUCCESS;
		}
		else if (strcmp(argv[i], "-outfile") == 0 or strcmp(argv[i], "-o") == 0)
		{
			if (outputFile != NULL)
			{
				cerr << "WARNING: Already used -outfile|-o with " << outputFile
					<< " before. Will override it with the last value specified with -outfile|-o."
					<< endl;
			}
			if (++i >= argc)
			{
				cerr << "ERROR: Missing an output filename." << endl << endl;
				PrintUsage();
				return CMDLINE_ERROR;
			}
			outputFile = argv[i];
		}
		else if (strcmp(argv[i], "-type") == 0 or strcmp(argv[i], "-t") == 0)
		{
			if (++i >= argc)
			{
				cerr << "ERROR: Missing a type specifier." << endl << endl;
				PrintUsage();
				return CMDLINE_ERROR;
			}
			// Now we need to parse the typename in the command line.
			if (strcmp(argv[i], "autodetect") == 0)
			{
				currentType = kUnknownDataBlock;
			}
			else
			{
				currentType = AliHLTMUONUtils::ParseCommandLineTypeString(argv[i]);
				if (currentType == kUnknownDataBlock)
				{
					cerr << "ERROR: Invalid type name '" << argv[i]
						<< "' specified for argument " << argv[i-1]
						<< "." << endl << endl;
					PrintUsage();
					return CMDLINE_ERROR;
				}
			}
		}
		else if (strcmp(argv[i], "-debug") == 0 or strcmp(argv[i], "-d") == 0)
		{
			maxLogging = true;
		}
		else if (strcmp(argv[i], "-check") == 0 or strcmp(argv[i], "-c") == 0)
		{
			checkData = true;
		}
		else if (strcmp(argv[i], "-cdbpath") == 0 or strcmp(argv[i], "-p") == 0)
		{
			if (pathSet)
			{
				cerr << "WARNING: Already used -cdbpath|-p with '" << gCDBPath
					<< "' before. Will override it with the last value specified with -cdbpath|-p."
					<< endl;
			}
			if (++i >= argc)
			{
				cerr << "ERROR: Missing the URL for the CDB path." << endl << endl;
				PrintUsage();
				return CMDLINE_ERROR;
			}
			gCDBPath = argv[i];
			pathSet = true;
		}
		else if (strcmp(argv[i], "-run") == 0 or strcmp(argv[i], "-r") == 0)
		{
			if (runSet)
			{
				cerr << "WARNING: Already used -run|-r with " << gRunNumber
					<< " before. Will override it with the last value specified with -run|-r."
					<< endl;
			}
			if (++i >= argc)
			{
				cerr << "ERROR: Missing the run number." << endl << endl;
				PrintUsage();
				return CMDLINE_ERROR;
			}
			
			char* cpErr = NULL;
			Int_t run = Int_t( strtol(argv[i], &cpErr, 0) );
			if (cpErr == NULL or *cpErr != '\0' or run < 0)
			{
				cerr << "ERROR: Cannot convert '" << argv[i] << "' to a valid run number."
					" Expected a positive integer value." << endl;
				return CMDLINE_ERROR;
			}

			gRunNumber = run;
			runSet = true;
		}
		else
		{
			assert( numOfFiles < argc );
			filenames[numOfFiles] = argv[i];
			
			AliHLTMUONDataBlockType typeToUse = currentType;
			// Work out what the file type is from the file's header
			// if the type is not yet known.
			if (typeToUse == kUnknownDataBlock)
			{
				int result = DecodeFileType(filenames[numOfFiles], typeToUse);
				if (result != EXIT_SUCCESS) return result;
				if (typeToUse == kUnknownDataBlock)
				{
					cerr << "ERROR: Could not decode what type of data"
						" is stored in the file:"
						<< filenames[numOfFiles] << endl;
					return PARSE_ERROR;
				}
			}
			
			filetypes[numOfFiles] = typeToUse;
			numOfFiles++;
		}
	}
	
	// Now check that we have at least one filename and all the flags we need.
	if (numOfFiles == 0)
	{
		cerr << "ERROR: Missing a file name. You must specify at least one file to process."
			<< endl << endl;
		PrintUsage();
		return CMDLINE_ERROR;
	}
	
	return EXIT_SUCCESS;
}


int main(int argc, const char** argv)
{
	// Test endianess of this machine during runtime and print warning if it is not
	// little endian.
	union
	{
		int dword;
		char byte[4];
	} endianTest;
	endianTest.dword = 0x1;
	if (endianTest.byte[0] != 0x1)
	{
		cerr << "!!!! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !!!!!" << endl;
		cerr << "!!!! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !!!!!" << endl;
		cerr << "!!                                                                     !!" << endl;
		cerr << "!! This is not a little endian machine, but dHLT raw data is normally  !!" << endl;
		cerr << "!! generated in little endian format. Unless you are looking at localy !!" << endl;
		cerr << "!! created simulated data, then this program will not show you correct !!" << endl;
		cerr << "!! output.                                                             !!" << endl;
		cerr << "!!                                                                     !!" << endl;
		cerr << "!!!! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !!!!!" << endl;
		cerr << "!!!! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !!!!!" << endl;
		cerr << endl;
	}


	int numOfFiles = 0;
	const char* outputFile = NULL;
	bool maxLogging = false;
	bool checkData = false;
	int returnCode = EXIT_SUCCESS;

	try
	{
		// There will be a maximum of 'argc' number of filenames possible.
		typedef const char* AnsiString;
		const char** filename = new AnsiString[argc];
		AliHLTMUONDataBlockType* filetype = new AliHLTMUONDataBlockType[argc];
		
		returnCode = ParseCommandLine(
				argc, argv, filename, filetype, numOfFiles,
				outputFile, maxLogging, checkData
			);
		if (outputFile == NULL)
		{
			outputFile = "output.root";
		}
		if (returnCode == EXIT_SUCCESS)
		{
			returnCode = RootifyFiles(
					filename, filetype, numOfFiles,
					outputFile, maxLogging, checkData
				);
		}
		
		delete [] filename;
		delete [] filetype;
	}
	catch (...)
	{
		cerr << "FATAL ERROR: An unknown exception occurred!" << endl << endl;
		returnCode = FATAL_ERROR;
	}
	
	return returnCode;
}

