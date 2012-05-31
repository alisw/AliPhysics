/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
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
 * @file   testCorruptorComponent.C
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   5 Aug 2010
 *
 * This macro is used to test the basic functionality of the AliHLTCorruptorComponent
 * class.
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TSystem.h"
#include "TClassTable.h"
#include "TBits.h"
#include "AliLog.h"
#include "AliHLTSystem.h"
#include "AliHLTConfiguration.h"
#include <fstream>
#include <cstdlib>
#endif

/**
 * Creates the input data for the test.
 * It is just a file of 256 bytes with all zeros.
 */
void GenerateInputData(bool debug = false)
{
	using namespace std;
	const char* filename = "corruptorInputTestFile.dat";
	fstream file(filename, ios::trunc | ios::out | ios::binary);
	if (! file)
	{
		if (debug) cerr << "ERROR: Could not create file " << filename << endl;
		return;
	}
	char buffer[256];
	memset(buffer, 0x0, sizeof(buffer));
	file.write(buffer, sizeof(buffer));
	if (! file)
	{
		if (debug) cerr << "ERROR: I/O error when writing to file " << filename << endl;
		return;
	}
	file.close();
}

/**
 * Routine to check that the filtering of the data blocks works with the corruptor
 * component. Only filtered blocks should be modified.
 */
void RunChainToCheckFiltering(bool debug = false)
{
	AliHLTSystem sys;
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	if (debug)
	{
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	else
	{
		sys.SetGlobalLoggingLevel(kHLTLogError);
	}
	
	const int numSources = 4;
	const char* fileTypeNames[numSources] = {
		"-datatype SOMEDATA HLT -dataspec 0x0",
		"-datatype MOREDATA HLT -dataspec 0x0",
		"-datatype SOMEDATA TPC -dataspec 0x0",
		"-datatype SOMEDATA HLT -dataspec 0xFF"
	};
	TString allSources = "";
	for (int i = 0; i < numSources; ++i)
	{
		TString sourceName = Form("source%d", i+1);
		if (allSources.Length() > 0) allSources += " ";
		allSources += sourceName;
		AliHLTConfiguration src(
			sourceName.Data(),
			"FilePublisher",
			"",
			Form("%s -datafile corruptorInputTestFile.dat", fileTypeNames[i])
		);
	}
	
	const int numProcessors = 4;
	const char* commandLine[numProcessors] = {
		"-datatype MOREDATA HLT -dataspec 0x0",
		"-datatype SOMEDATA TPC -dataspec 0x0",
		"-datatype SOMEDATA HLT -dataspec 0xFF",
		"-dataspec 0xFF -typeid '*******' -typeid MOREDATA -dataspec 0x0"
	};
	for (int i = 0; i < numProcessors; ++i)
	{
		TString processorName = Form("processor%d", i+1);
		AliHLTConfiguration prc(
			processorName.Data(),
			"CorruptorComponent",
			allSources.Data(),
			commandLine[i]
		);
		TString sinkName = Form("sink%d", i+1);
		AliHLTConfiguration snk(
			sinkName.Data(),
			"FileWriter",
			processorName.Data(),
			Form("-specfmt -datafile corruptorOutputTestFile_sink%d.dat", i+1)
		);
		sys.BuildTaskList(sinkName.Data());
	}
	
	sys.Run(1); // Run for 1 event.
}

/**
 * Reads a file into a buffer.
 * \param [in]  filename  The name of the file to read.
 * \param [out] buffer  The buffer which is created to store the contents.
 *                      Must be deleted with 'delete [] buffer' by caller if not NULL.
 * \param [out] size  The size of the buffer created in bytes.
 * \param [in]  debug  Indicates if debug messages should be printed.
 * \returns true if the file was read, the buffer created and filled; false otherwise.
 * \note The caller becomes the owner of the allocated buffer.
 */
bool ReadFile(const char* filename, char*& buffer, size_t& size, bool debug = false)
{
	buffer = NULL;
	size = 0;
	using namespace std;
	fstream file(filename, ios::in | ios::binary);
	if (! file)
	{
		if (debug) cerr << "ERROR: Could not open file " << filename << endl;
		return false;
	}
	file.seekg(0, std::ios::end);
	size_t filesize = file.tellg();
	file.seekg(0, std::ios::beg);
	if (! file)
	{
		if (debug) cerr << "ERROR: Could not get file size for " << filename << endl;
		return false;
	}
	buffer = new char[filesize];
	if (buffer == NULL)
	{
		if (debug) cerr << "ERROR: Cannot allocate more memory for buffers." << endl;
		return false;
	}
	size = filesize;
	file.read(buffer, size);
	if (! file)
	{
		if (debug) cerr << "ERROR: I/O error when reading from file " << filename << endl;
		return false;
	}
	return true;
}

/**
 * Checks if the output data is as expected for the filtering.
 * Only the output data blocks that were filtered should be modified.
 */
bool CheckFilteringOutput(bool debug = false)
{
	const int numBuffers = 5;
	int sinknum[numBuffers] = {1, 2, 3, 4, 4};
	const char* blocknum[numBuffers] = {
		"0x00",
		"0x00",
		"0x00",
		"0x00",
		"0x01"
	};
	const char* blocktype[numBuffers] = {
		"HLT:MOREDATA",
		"TPC:SOMEDATA",
		"HLT:SOMEDATA",
		"HLT:SOMEDATA",
		"HLT:MOREDATA"
	};
	const char* blockspec[numBuffers] = {
		"0x00000000",
		"0x00000000",
		"0x000000ff",
		"0x000000ff",
		"0x00000000"
	};
	char* buffer[numBuffers];
	size_t size[numBuffers];
	for (int i = 0; i < numBuffers; ++i)
	{
		const char* filename = Form("corruptorOutputTestFile_sink%d_0x00000000_%s_%s_%s.dat",
					    sinknum[i], blocknum[i], blocktype[i], blockspec[i]
					   );
		if (! ReadFile(filename, buffer[i], size[i], debug))
		{
			if (! debug)
			{
				cerr << "ERROR: Filtering chain test did not generate correct output files."
					<< " Run with 'debug = true' for more details." << endl;
			}
			return false;
		}
		delete [] buffer[i];
	}
	
	return true;
}

/**
 * Routine to check the single bit flip option.
 */
void RunChainToCheckSingleFlips(bool debug = false)
{
	AliHLTSystem sys;
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	if (debug)
	{
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	else
	{
		sys.SetGlobalLoggingLevel(kHLTLogError);
	}
	
	AliHLTConfiguration src(
		"source",
		"FilePublisher",
		"",
		"-datatype SOMEDATA HLT -dataspec 0x0 -datafile corruptorInputTestFile.dat"
	);
	AliHLTConfiguration prc(
		"processor",
		"CorruptorComponent",
		"source",
		"-seed 123 -datatype SOMEDATA HLT -range 128 max -alignment 4 -errorcount 5 5 -singleflips"
	);
	AliHLTConfiguration snk(
		"sink",
		"FileWriter",
		"processor",
		"-datafile corruptorOutputTestFile_single.dat"
	);
	
	sys.BuildTaskList("sink");
	sys.Run(1); // Run for 1 event.
}

/**
 * Checks if the output data is as expected for the single bit flip option.
 */
bool CheckSingleFlipOutput(bool debug = false)
{
	char* buffer = NULL;
	size_t size = 0;
	if (! ReadFile("corruptorOutputTestFile_single_0x00000000_0x00_HLT:SOMEDATA.dat", buffer, size, debug))
	{
		if (! debug)
		{
			cerr << "ERROR: Single flips chain test did not generate correct"
				<< " output files. Run with 'debug = true' for more details."
				<< endl;
		}
		return false;
	}
	if (size != 256)
	{
		cerr << "ERROR: The output file is not the expected 256 bytes in size"
			<< " when testing the -singleflips option." << endl;
		delete [] buffer;
		return false;
	}
	// Check that the first 128 bytes are zero.
	for (int i = 0; i < 128; ++i)
	{
		if (buffer[i] != 0x0)
		{
			cerr << "ERROR: The first 128 bytes were not left as zeros when testing"
				<< " the -singleflips option." << endl;
			delete [] buffer;
			return false;
		}
	}
	// Check that the correct number of bit flips happened in the last 128 bytes.
	TBits bits(128*8);
	bits.Set(128*8, buffer+128);
	UInt_t bitcount = bits.CountBits();
	if (bitcount != 5)
	{
		cerr << "ERROR: When testing the -singleflips option,"
			<< " the number of bits flipped in the output buffer was "
			<< bitcount << ", but we expect exactly 5 bit flips."
			<< endl;
		delete [] buffer;
		return false;
	}
	// Check that the bit flips are only on 4 bit aligned addresses.
	for (UInt_t j = 128*8; j < 256*8; ++j)
	{
		if (bits[j] && (j & 0x7) != 0)
		{
			cerr << "ERROR: When testing the -singleflips option, bit " << j
				<< " was flipped in the output buffer,"
				<< " but it is not aligned to a 4 bit word boundary."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	
	delete [] buffer;
	return true;
}

/**
 * Routine to check the burst error option.
 */
void RunChainToCheckBurstErrors(bool debug = false)
{
	AliHLTSystem sys;
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	if (debug)
	{
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	else
	{
		sys.SetGlobalLoggingLevel(kHLTLogError);
	}
	
	AliHLTConfiguration src(
		"source",
		"FilePublisher",
		"",
		"-datatype SOMEDATA HLT -dataspec 0x0 -datafile corruptorInputTestFile.dat"
	);
	AliHLTConfiguration prc(
		"processor",
		"CorruptorComponent",
		"source",
		"-seed 123 -datatype SOMEDATA HLT -range 0 256 -alignment 32 -errorcount 3 3 -bursterrors 16 8"
	);
	AliHLTConfiguration snk(
		"sink",
		"FileWriter",
		"processor",
		"-datafile corruptorOutputTestFile_burst.dat"
	);
	
	sys.BuildTaskList("sink");
	sys.Run(1); // Run for 1 event.
}

/**
 * Checks if the output data is as expected for the burst errors.
 */
bool CheckBurstErrorOutput(bool debug = false)
{
	char* buffer = NULL;
	size_t size = 0;
	if (! ReadFile("corruptorOutputTestFile_burst_0x00000000_0x00_HLT:SOMEDATA.dat", buffer, size, debug))
	{
		if (! debug)
		{
			cerr << "ERROR: Burst error chain test did not generate correct"
				<< " output files. Run with 'debug = true' for more details."
				<< endl;
		}
		return false;
	}
	if (size != 256)
	{
		cerr << "ERROR: The output file is not the expected 256 bytes in size"
			<< " when testing the -bursterrors option." << endl;
		delete [] buffer;
		return false;
	}
	TBits bits(256*8);
	bits.Set(256*8, buffer);
	// Check that the bit flips are only on the lower part of the 32 bit words.
	for (UInt_t j = 0; j < 256*8; ++j)
	{
		if (bits[j] && (j % 32) >= 16)
		{
			cerr << "ERROR: When testing the -bursterrors option, bit " << j
				<< " was flipped in the output buffer, but only the lower"
				<< " part of any 32 bit word is supposed to be corrupted."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	
	delete [] buffer;
	return true;
}

/**
 * Routine to check the replacement options.
 */
void RunChainToCheckReplaceErrors(bool debug = false)
{
	AliHLTSystem sys;
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	if (debug)
	{
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	else
	{
		sys.SetGlobalLoggingLevel(kHLTLogError);
	}
	
	AliHLTConfiguration src1(
		"source1",
		"FilePublisher",
		"",
		"-datatype DATAAAAA HLT -dataspec 0x0 -datafile corruptorInputTestFile.dat"
	);
	AliHLTConfiguration prc1(
		"processor1",
		"CorruptorComponent",
		"source1",
		"-datatype DATAAAAA HLT -range 128:4 128:4 -replace 0x1/1 0xF 0xF/3 -seed 1"
	);
	AliHLTConfiguration src2(
		"source2",
		"FilePublisher",
		"",
		"-datatype DATABBBB HLT -dataspec 0x0 -datafile corruptorInputTestFile.dat"
	);
	AliHLTConfiguration prc2(
		"processor2",
		"CorruptorComponent",
		"source2",
		"-datatype DATABBBB HLT -seed 123 -range 8 8 -errorcount 1 1 -replace-random 64 64"
	);
	AliHLTConfiguration snk(
		"sink",
		"FileWriter",
		"processor1 processor2",
		"-datafile corruptorOutputTestFile_replace.dat"
	);
	
	sys.BuildTaskList("sink");
	sys.Run(1); // Run for 1 event.
}

/**
 * Checks if the output data is as expected for the replacement commands.
 */
bool CheckReplaceErrorsOutput(bool debug = false)
{
	char* buffer = NULL;
	size_t size = 0;
	if (! ReadFile("corruptorOutputTestFile_replace_0x00000000_0x01_HLT:DATAAAAA.dat", buffer, size, debug))
	{
		if (! debug)
		{
			cerr << "ERROR: Replace errors chain test did not generate correct"
				<< " output files. Run with 'debug = true' for more details."
				<< endl;
		}
		return false;
	}
	if (size != 256)
	{
		cerr << "ERROR: The output file is not the expected 256 bytes in size"
			<< " when testing the -replace option." << endl;
		delete [] buffer;
		return false;
	}
	// Check that the correct bits were set. All bits zero except 128:4 to 129:4
	for (UInt_t j = 0; j < 128; ++j)
	{
		if (buffer[j] != 0)
		{
			cerr << "ERROR: When testing the -replace option we found bits set in byte"
				<< j << " but the byte should be zero."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	if (AliHLTUInt8_t(buffer[128]) != 0xF0 || AliHLTUInt8_t(buffer[129]) != 0x0F)
	{
		cerr << "ERROR: When testing the -replace option we found a value of "
			<< Form("0x%2.2X", int(AliHLTUInt8_t(buffer[128]))) << " for byte 128 and "
			<< Form("0x%2.2X", int(AliHLTUInt8_t(buffer[129]))) << " for byte 129,"
			<< " but we expected values of 0xF0 and 0x0F respectively."
			<< endl;
		delete [] buffer;
		return false;
	}
	for (UInt_t j = 130; j < 256; ++j)
	{
		if (buffer[j] != 0)
		{
			cerr << "ERROR: When testing the -replace option we found bits set in byte"
				<< j << " but the byte should be zero."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	delete [] buffer;
	buffer = NULL;
	size = 0;
	// Now test the next buffer...
	if (! ReadFile("corruptorOutputTestFile_replace_0x00000000_0x00_HLT:DATABBBB.dat", buffer, size, debug))
	{
		if (! debug)
		{
			cerr << "ERROR: Replace errors chain test did not generate correct"
				<< " output files. Run with 'debug = true' for more details."
				<< endl;
		}
		return false;
	}
	if (size != 256)
	{
		cerr << "ERROR: The output file is not the expected 256 bytes in size"
			<< " when testing the -replace-random option." << endl;
		delete [] buffer;
		return false;
	}
	// Check that only bits in the second 64 bit word were set.
	for (UInt_t j = 0; j < 8; ++j)
	{
		if (buffer[j] != 0x0)
		{
			cerr << "ERROR: When testing the -replace-random option we found bits set in byte"
				<< j << " but the byte should be zero."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	if ((reinterpret_cast<AliHLTUInt64_t*>(buffer))[1] == 0x0)
	{
		cerr << "ERROR: When testing the -replace-random option we found no bits set in bytes 8 to 15." << endl;
		delete [] buffer;
		return false;
	}
	for (UInt_t j = 16; j < 256; ++j)
	{
		if (buffer[j] != 0x0)
		{
			cerr << "ERROR: When testing the -replace-random option we found bits set in byte"
				<< j << " but the byte should be zero."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	delete [] buffer;
	return true;
}

/**
 * Routine to check the insertion options.
 */
void RunChainToCheckInsertErrors(bool debug = false)
{
	AliHLTSystem sys;
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	if (debug)
	{
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	else
	{
		sys.SetGlobalLoggingLevel(kHLTLogError);
	}
	
	AliHLTConfiguration src1(
		"source1",
		"FilePublisher",
		"",
		"-datatype DATAAAAA HLT -dataspec 0x0 -datafile corruptorInputTestFile.dat"
	);
	AliHLTConfiguration prc1(
		"processor1",
		"CorruptorComponent",
		"source1",
		"-datatype DATAAAAA HLT -range 128 128 -errorcount 1 1 -insert 0xCBA/12"
	);
	AliHLTConfiguration src2(
		"source2",
		"FilePublisher",
		"",
		"-datatype DATABBBB HLT -dataspec 0x0 -datafile corruptorInputTestFile.dat"
	);
	AliHLTConfiguration prc2(
		"processor2",
		"CorruptorComponent",
		"source2",
		"-datatype DATABBBB HLT -seed 123 -range max max -errorcount 1 1 -insert-random 64 64"
	);
	AliHLTConfiguration snk(
		"sink",
		"FileWriter",
		"processor1 processor2",
		"-datafile corruptorOutputTestFile_insert.dat"
	);
	
	sys.BuildTaskList("sink");
	sys.Run(1); // Run for 1 event.
}

/**
 * Checks if the output data is as expected for the insertion commands.
 */
bool CheckInsertErrorsOutput(bool debug = false)
{
	char* buffer = NULL;
	size_t size = 0;
	if (! ReadFile("corruptorOutputTestFile_insert_0x00000000_0x01_HLT:DATAAAAA.dat", buffer, size, debug))
	{
		if (! debug)
		{
			cerr << "ERROR: Insert errors chain test did not generate correct"
				<< " output files. Run with 'debug = true' for more details."
				<< endl;
		}
		return false;
	}
	if (size != 258)
	{
		cerr << "ERROR: The output file is not the expected 258 bytes in size"
			<< " when testing the -insert option." << endl;
		delete [] buffer;
		return false;
	}
	// Check that the correct bits were set. All bits zero except 128 to 129
	for (UInt_t j = 0; j < 128; ++j)
	{
		if (buffer[j] != 0)
		{
			cerr << "ERROR: When testing the -insert option we found bits set in byte"
				<< j << " but the byte should be zero."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	if (AliHLTUInt8_t(buffer[128]) != 0xBA || AliHLTUInt8_t(buffer[129]) != 0x0C)
	{
		cerr << "ERROR: When testing the -insert option we found a value of "
			<< Form("0x%2.2X", int(AliHLTUInt8_t(buffer[128]))) << " for byte 128 and "
			<< Form("0x%2.2X", int(AliHLTUInt8_t(buffer[129]))) << " for byte 129,"
			<< " but we expected values of 0xBA and 0x0C respectively."
			<< endl;
		delete [] buffer;
		return false;
	}
	for (UInt_t j = 130; j < 258; ++j)
	{
		if (buffer[j] != 0)
		{
			cerr << "ERROR: When testing the -insert option we found bits set in byte"
				<< j << " but the byte should be zero."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	delete [] buffer;
	buffer = NULL;
	size = 0;
	// Now test the next buffer...
	if (! ReadFile("corruptorOutputTestFile_insert_0x00000000_0x00_HLT:DATABBBB.dat", buffer, size, debug))
	{
		if (! debug)
		{
			cerr << "ERROR: Insert errors chain test did not generate correct"
				<< " output files. Run with 'debug = true' for more details."
				<< endl;
		}
		return false;
	}
	if (size != 264)
	{
		cerr << "ERROR: The output file is not the expected 264 bytes in size"
			<< " when testing the -insert-random option." << endl;
		delete [] buffer;
		return false;
	}
	// Check that only bits in the last 64 bit word were set.
	for (UInt_t j = 0; j < 256; ++j)
	{
		if (buffer[j] != 0x0)
		{
			cerr << "ERROR: When testing the -insert-random option we found bits set in byte"
				<< j << " but the byte should be zero."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	if ((reinterpret_cast<AliHLTUInt64_t*>(buffer))[32] == 0x0)
	{
		cerr << "ERROR: When testing the -insert-random option we found no bits set in bytes 256 to 263." << endl;
		delete [] buffer;
		return false;
	}
	delete [] buffer;
	return true;
}

/**
 * Routine to check the remove option.
 */
void RunChainToCheckRemoveErrors(bool debug = false)
{
	AliHLTSystem sys;
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	if (debug)
	{
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	else
	{
		sys.SetGlobalLoggingLevel(kHLTLogError);
	}
	
	AliHLTConfiguration src(
		"source",
		"FilePublisher",
		"",
		"-datatype SOMEDATA HLT -dataspec 0x0 -datafile corruptorInputTestFile.dat"
	);
	AliHLTConfiguration prc(
		"processor",
		"CorruptorComponent",
		"source",
		"-datatype SOMEDATA HLT -range 128 128 -errorcount 1 1 -replace 0x78563412/32 -remove 16 16 -range max max -insert removed"
	);
	AliHLTConfiguration snk(
		"sink",
		"FileWriter",
		"processor",
		"-datafile corruptorOutputTestFile_remove.dat"
	);
	
	sys.BuildTaskList("sink");
	sys.Run(1); // Run for 1 event.
}

/**
 * Checks if the output data is as expected for the remove command.
 */
bool CheckRemoveErrorsOutput(bool debug = false)
{
	char* buffer = NULL;
	size_t size = 0;
	if (! ReadFile("corruptorOutputTestFile_remove_0x00000000_0x00_HLT:SOMEDATA.dat", buffer, size, debug))
	{
		if (! debug)
		{
			cerr << "ERROR: Remove errors chain test did not generate correct"
				<< " output file. Run with 'debug = true' for more details."
				<< endl;
		}
		return false;
	}
	if (size != 256)
	{
		cerr << "ERROR: The output file is not the expected 256 bytes in size"
			<< " when testing the -remove option." << endl;
		delete [] buffer;
		return false;
	}
	// Check that the correct bits were removed and set.
	for (UInt_t j = 0; j < 128; ++j)
	{
		if (buffer[j] != 0)
		{
			cerr << "ERROR: When testing the -remove option we found bits set in byte"
				<< j << " but the byte should be zero."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	if (AliHLTUInt8_t(buffer[128]) != 0x56 || AliHLTUInt8_t(buffer[129]) != 0x78)
	{
		cerr << "ERROR: When testing the -remove option we found a value of "
			<< Form("0x%2.2X", int(AliHLTUInt8_t(buffer[128]))) << " for byte 128 and "
			<< Form("0x%2.2X", int(AliHLTUInt8_t(buffer[129]))) << " for byte 129,"
			<< " but we expected values of 0x56 and 0x78 respectively."
			<< endl;
		delete [] buffer;
		return false;
	}
	for (UInt_t j = 130; j < 254; ++j)
	{
		if (buffer[j] != 0)
		{
			cerr << "ERROR: When testing the -remove option we found bits set in byte"
				<< j << " but the byte should be zero."
				<< endl;
			delete [] buffer;
			return false;
		}
	}
	if (AliHLTUInt8_t(buffer[254]) != 0x12 || AliHLTUInt8_t(buffer[255]) != 0x34)
	{
		cerr << "ERROR: When testing the -remove option we found a value of "
			<< Form("0x%2.2X", int(AliHLTUInt8_t(buffer[254]))) << " for byte 128 and "
			<< Form("0x%2.2X", int(AliHLTUInt8_t(buffer[255]))) << " for byte 129,"
			<< " but we expected values of 0x12 and 0x34 respectively."
			<< endl;
		delete [] buffer;
		return false;
	}
	delete [] buffer;
	return true;
}

/**
 * Routine to check the relative address option for "-range".
 */
void RunChainToCheckRelativeAddress(bool debug = false)
{
	AliHLTSystem sys;
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	if (debug)
	{
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	else
	{
		sys.SetGlobalLoggingLevel(kHLTLogError);
	}
	
	AliHLTConfiguration src(
		"source",
		"FilePublisher",
		"",
		"-datatype SOMEDATA HLT -dataspec 0x0 -datafile corruptorInputTestFile.dat"
	);
	AliHLTConfiguration prc(
		"processor",
		"CorruptorComponent",
		"source",
		"-datatype SOMEDATA HLT -errorcount 1 1 -range max-4 max-4 -replace 0xFF/8 -range min+0:4 min+0:4 -replace 0xBA/8"
	);
	AliHLTConfiguration snk(
		"sink",
		"FileWriter",
		"processor",
		"-datafile corruptorOutputTestFile_reladdr.dat"
	);
	
	sys.BuildTaskList("sink");
	sys.Run(1); // Run for 1 event.
}

/**
 * Checks if the output data is as expected for the remove command.
 */
bool CheckRelativeAddressOutput(bool debug = false)
{
	char* buffer = NULL;
	size_t size = 0;
	if (! ReadFile("corruptorOutputTestFile_reladdr_0x00000000_0x00_HLT:SOMEDATA.dat", buffer, size, debug))
	{
		if (! debug)
		{
			cerr << "ERROR: Relative address chain test did not generate correct"
				<< " output file. Run with 'debug = true' for more details."
				<< endl;
		}
		return false;
	}
	if (size != 256)
	{
		cerr << "ERROR: The output file is not the expected 256 bytes in size"
			<< " when testing the relative address option." << endl;
		delete [] buffer;
		return false;
	}
	// Check that the correct bits were set.
	char refBuf[256];
	memset(refBuf, 0x0, sizeof(refBuf));
	refBuf[0] = 0xA0;
	refBuf[1] = 0x0B;
	refBuf[252] = 0xFF;
	for (UInt_t j = 0; j < 256; ++j)
	{
		if (buffer[j] != refBuf[j])
		{
			cerr << "ERROR: When testing the relative address option we find byte "
				<< j << " has a value of "
				<< Form("0x%2.2X", int(AliHLTUInt8_t(buffer[j])))
				<< ", but we expected a value of "
				<< Form("0x%2.2X", int(AliHLTUInt8_t(refBuf[j])))
				<< "." << endl;
			delete [] buffer;
			return false;
		}
	}
	delete [] buffer;
	return true;
}

/**
 * This is the top level testing method which calls individual tests.
 * \returns true if all tests succeeded and false otherwise.
 */
bool testAliHLTCorruptorComponent(bool debug = false)
{
	if (debug)
	{
		AliLog::SetGlobalLogLevel(AliLog::kMaxType);
	}
	else
	{
		// Done here to prevent output from AliHLTSystem.
		AliLog::SetGlobalLogLevel(AliLog::kError);
	}
	
	if (gClassTable->GetID("AliHLTCorruptorComponent") < 0)
	{
		gSystem->Load("libAliHLTUtil.so");
	}
	
	GenerateInputData(debug);
	RunChainToCheckFiltering(debug);
	if (! CheckFilteringOutput(debug)) return false;
	RunChainToCheckSingleFlips(debug);
	if (! CheckSingleFlipOutput(debug)) return false;
	RunChainToCheckBurstErrors(debug);
	if (! CheckBurstErrorOutput(debug)) return false;
	RunChainToCheckReplaceErrors(debug);
	if (! CheckReplaceErrorsOutput(debug)) return false;
	RunChainToCheckInsertErrors(debug);
	if (! CheckInsertErrorsOutput(debug)) return false;
	RunChainToCheckRemoveErrors(debug);
	if (! CheckRemoveErrorsOutput(debug)) return false;
	RunChainToCheckRelativeAddress(debug);
	if (! CheckRelativeAddressOutput(debug)) return false;
	
	// Cleanup all temporary files generated.
	gSystem->Exec("rm -f corruptorInputTestFile.dat corruptorOutputTestFile*.dat");
	return true;
}

#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	bool resultOk = testAliHLTCorruptorComponent();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__

