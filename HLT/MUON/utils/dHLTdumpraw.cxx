/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author:                                                                *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/**
 * @file   dHLTdumpraw.cxx
 * @author Artur Szostak <artursz@iafrica.com>,
 *         Seforo Mohlalisi <seforomohlalisi@yahoo.co.uk>
 * @date   
 * @brief  Command line utility to dump dHLT's internal raw data blocks.
 */

// We define NDEBUG for the AliHLTMUONDataBlockReader.h header file since this
// program by definition handles corrupt data. So we do not need the assertions
// in the AliHLTMUONDataBlockReader class to be checked.
#define NDEBUG
#include "AliHLTMUONDataBlockReader.h"
#undef NDEBUG
#include "AliHLTMUONUtils.h"

/*TODO: fix this. Need a platform independant way of checking the endian encoding.
 * This does not want to work on Apple Mac OS compiler: i686-apple-darw
#include <endian.h>
#ifndef LITTLE_ENDIAN
#error Handling of internal data for non little endian machines not yet implemented.
#endif
*/

#include <cstdlib>
#include <cassert>
#include <new>
#include <fstream>

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::showbase;
using std::noshowbase;
using std::hex;
using std::dec;

#include <iomanip>
using std::setw;
using std::left;
using std::internal;


#define CMDLINE_ERROR 1
#define PARSE_ERROR 2
#define SYSTEM_ERROR 3
#define FATAL_ERROR 4


void PrintRubbishData(AliHLTUInt32_t offset, const char* padByte, AliHLTUInt32_t padCount)
{
	if (padCount == 0) return;
	
	cerr << "ERROR: Found the following unexpected rubbish data at the"
		" end of the data block:" << endl;
	cerr << "Byte #\tValue\tCharacter" << endl;
	for (AliHLTUInt32_t i = 0; i < padCount; i++)
	{
		short value = short(padByte[i]) & 0xFF;
		char character = padByte[i];
		cerr << offset + i + 1 << "\t"
			<< noshowbase << hex << "0x" << value << dec << "\t"
			<< character << endl;
	}
}


template <typename FieldType>
int CheckHeaderField(
		FieldType& field, const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	const char* fieldptr = reinterpret_cast<const char*>(&field);
	const char* endptr = buffer + bufferSize;
	AliHLTUInt32_t bufferRemaining = endptr > fieldptr ? endptr - fieldptr : 0;
	
	if (bufferRemaining < sizeof(field))
	{
		cout << "..." << endl; // We may be half way through printing a line so end it.
		cerr << "ERROR: The data block is too short. The header is corrupt." << endl;
		if (continueParse)
		{
			AliHLTUInt32_t offset = fieldptr - buffer;
			PrintRubbishData(offset, fieldptr, bufferRemaining);
		}
		return PARSE_ERROR;
	}
	return EXIT_SUCCESS;
}


template <typename FieldType>
int CheckField(
		FieldType& field, const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	const char* fieldptr = reinterpret_cast<const char*>(&field);
	const char* endptr = buffer + bufferSize;
	AliHLTUInt32_t bufferRemaining = endptr > fieldptr ? endptr - fieldptr : 0;
	
	if (bufferRemaining < sizeof(field))
	{
		cout << "..." << endl; // We may be half way through printing a line so end it.
		cerr << "ERROR: The data block is too short. The data is corrupt." << endl;
		if (continueParse)
		{
			AliHLTUInt32_t offset = fieldptr - buffer;
			PrintRubbishData(offset, fieldptr, bufferRemaining);
		}
		return PARSE_ERROR;
	}
	return EXIT_SUCCESS;
}


template <typename BlockType>
int CheckCommonHeader(
		BlockType& block, const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	int result = EXIT_SUCCESS;

	// Check the fRecordWidth field in the common header.
	if (block.CommonBlockHeader().fRecordWidth !=
		sizeof(typename BlockType::ElementType))
	{
		cerr << "ERROR: The record width found in the header is incorrect."
			" Found a record width of "
			<< block.CommonBlockHeader().fRecordWidth << " bytes, but expected"
			" a value of " << sizeof(typename BlockType::ElementType)
			<< " bytes." << endl;
		result = PARSE_ERROR;
		if (not continueParse) return result;
	}
	
	if (not block.BufferSizeOk())
	{
		cerr << "ERROR: The size of the file is incorrect. It is "
			<< bufferSize << " bytes big, but according"
			" to the data block header it should be " << block.BytesUsed()
			<< " bytes." << endl;
		result = PARSE_ERROR;
		if (not continueParse) return result;
	}
	
	return result;
}


template <typename BlockType>
AliHLTUInt32_t CalculateNEntries(BlockType& block, unsigned long bufferSize)
{
	// Calculate how many entries we can display. If the buffer size is correct
	// we just use the number of entries the block specifies. Otherwise we need
	// to calculate it from the buffer size.
	AliHLTUInt32_t nentries;
	if (block.BytesUsed() == bufferSize)
	{
		nentries = block.Nentries();
	}
	else
	{
		AliHLTInt32_t dataSize = bufferSize
			- sizeof(typename BlockType::HeaderType);
		nentries = dataSize / sizeof(typename BlockType::ElementType);
		if (dataSize % sizeof(typename BlockType::ElementType) > 0)
			nentries++;
	}
	return nentries;
}


int DumpRecHitStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONRecHitStruct* hit,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(hit->fX, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << hit->fX << setw(0);

	result = CheckField(hit->fY, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << hit->fY << setw(0);

	result = CheckField(hit->fZ, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << hit->fZ << setw(0) << endl;

	return result;
}


int DumpRecHitsBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	int result = EXIT_SUCCESS;
	AliHLTMUONRecHitsBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	cout << " X (cm)     | Y (cm)     | Z (cm)" << endl;
	cout << "---------------------------------------" << endl;
	const AliHLTMUONRecHitStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		int subResult = DumpRecHitStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpTriggerRecordStruct(
               const char* buffer, unsigned long bufferSize, 
               const AliHLTMUONTriggerRecordStruct* record,
               bool continueParse
       )
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(record->fId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Trigger Record ID: " << record->fId <<endl;
	
	result = CheckField(record->fFlags, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Flags: " << showbase << hex << record->fFlags << dec;
		
	// Print the individual trigger bits.
	AliHLTMUONParticleSign sign;
	bool hitset[4];
	AliHLTMUONUtils::UnpackTriggerRecordFlags(record->fFlags, sign, hitset);
	cout << " [Sign: " << sign << ", Hits set on chambers: ";
	bool first = true;
	for (AliHLTUInt32_t i = 0; i < 4; i++)
	{
		if (hitset[i])
		{
			cout << (first ? "" : ", ") << i+11;
			first = false;
		}
	}
	cout << (first ? "none]" : "]") << endl;

	result = CheckField(record->fPx, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Momentum: (px = " << record->fPx << ", ";

	result = CheckField(record->fPy, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "py = " << record->fPy << ", ";

	result = CheckField(record->fPz, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "pz = " << record->fPz << ") GeV/c"<<endl;
	
	cout << "Track hits:" << endl;
	cout << "Chamber | X (cm)     | Y (cm)     | Z (cm)" << endl;
	cout << "------------------------------------------------" << endl;
	const AliHLTMUONRecHitStruct* hit = &record->fHit[0];
	for(AliHLTUInt32_t ch = 0; ch < 4; ch++)
	{
	        cout << setw(10) << left << ch + 11 << setw(0);
		result = DumpRecHitStruct(buffer, bufferSize, hit++, continueParse);
		if (result != EXIT_SUCCESS) return result;
	}

	return result;

}


int DumpTriggerRecordsBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	AliHLTMUONTriggerRecordsBlockReader block(buffer, bufferSize);
	
	int result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	const AliHLTMUONTriggerRecordStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		cout << "============================== Trigger Record number " << i+1
			<< " of " << nentries << " ==============================" << endl;
		int subResult = DumpTriggerRecordStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return EXIT_SUCCESS;
}


int DumpTrigRecInfoStruct(const char* buffer, unsigned long bufferSize, 
                          const AliHLTMUONTrigRecInfoStruct* debuginfo, 
                          bool continueParse
                          )
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(debuginfo->fTrigRecId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(22) << left << debuginfo->fTrigRecId << setw(0);

	result = CheckField(debuginfo->fDetElemId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(20) << left << debuginfo->fDetElemId << setw(0);
        
	result = CheckField(debuginfo->fZmiddle, buffer, bufferSize, continueParse);
	if(result != EXIT_SUCCESS) return result;
	cout << setw(30) << left << debuginfo->fZmiddle << setw(0);

	result = CheckField(debuginfo->fBl, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout <<debuginfo->fBl << setw(0) << endl;

	return result;
}


int DumpTrigRecsDebugBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	AliHLTMUONTrigRecsDebugBlockReader block(buffer, bufferSize);
	
	int result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	cout << "Trigger Record ID  | Detector ID  | Momentum X Component (Gev/c) | Integrated Magnetic Field (T.m)" << endl;
	cout << "--------------------------------------------------------------------------------------------------" << endl;
	const AliHLTMUONTrigRecInfoStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		int subResult = DumpTrigRecInfoStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return EXIT_SUCCESS;
}


int DumpTriggerChannelStruct(const char* buffer, unsigned long bufferSize, 
                      const AliHLTMUONTriggerChannelStruct* triggerchannel, 
                      bool continueParse
                      )
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(triggerchannel->fTrigRecId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(25) << left << triggerchannel->fTrigRecId << setw(0);

	result = CheckField(triggerchannel->fChamber, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << triggerchannel->fChamber << setw(0);

	result = CheckField(triggerchannel->fSignal, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(10) << left << triggerchannel->fSignal << setw(0);

	result = CheckField(triggerchannel->fRawDataWord, buffer, bufferSize, continueParse);
	if(result != EXIT_SUCCESS) return result;
	cout << showbase << hex << triggerchannel->fRawDataWord << dec << setw(0) << endl;
	return result;
}


int DumpTriggerChannelsBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
        int result = EXIT_SUCCESS;
	AliHLTMUONTriggerChannelsBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	cout << " Trigger Record ID   | Chamber    | Signal   | Raw Data Word " << endl;
	cout << "--------------------------------------------------------------" << endl;
	const AliHLTMUONTriggerChannelStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		int subResult = DumpTriggerChannelStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}

	return result;
}


int DumpClusterStruct(
		      const char* buffer, unsigned long bufferSize,
                      const AliHLTMUONClusterStruct* cluster,
                      bool continueParse
       )
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(cluster->fId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "cluster->fId: " << cluster->fId << "\t";

	result = CheckField(cluster->fDetElemId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "cluster->fDetElemId: " << cluster->fDetElemId << "\t";

	result = CheckField(cluster->fNchannels, buffer, bufferSize, continueParse);
	if(result != EXIT_SUCCESS) return result;
	cout << "cluster->fNchannels: " << cluster->fNchannels <<endl;

	cout << " Corresponding Hit: "<< endl;
	cout << " X (cm)     | Y (cm)     | Z (cm)" << endl;
	cout << "---------------------------------------" << endl;
	const AliHLTMUONRecHitStruct * hit = & cluster->fHit;
	result = DumpRecHitStruct(buffer, bufferSize, hit, continueParse);

	return result;
}


int DumpClustersBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
        int result = EXIT_SUCCESS;
	AliHLTMUONClustersBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	const AliHLTMUONClusterStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		cout << " ===================================================== Cluster Number "
			<< i+1 << "==================================================" << endl; 
		int subResult = DumpClusterStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}	

	return result;
}


int DumpChannelStruct(
               const char* buffer, unsigned long bufferSize,
               const AliHLTMUONChannelStruct* channel,
               bool continueParse                 
       )
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(channel->fClusterId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(16) << left << channel->fClusterId << setw(0);

	result = CheckField(channel->fManu, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(16) << left << channel->fManu << setw(0);

	result = CheckField(channel->fChannelAddress, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(19) << left << channel->fChannelAddress << setw(0);

	result = CheckField(channel->fSignal, buffer, bufferSize, continueParse);
	if(result != EXIT_SUCCESS) return result;
	cout << setw(16) << left << channel->fSignal << setw(0);

	result = CheckField(channel->fRawDataWord, buffer, bufferSize, continueParse);
	if(result != EXIT_SUCCESS) return result;
	cout << showbase << hex << channel->fRawDataWord << dec << setw(0) <<endl;

	return result;
}


int DumpChannelsBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
        int result = EXIT_SUCCESS;
	AliHLTMUONChannelsBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	cout << "Cluster Id  | Manu Address  | Channel Address  | Signal Value  | Raw Data Word " <<endl;
	cout << "-------------------------------------------------------------------------------" <<endl;
	const AliHLTMUONChannelStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{ 
		int subResult = DumpChannelStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	return EXIT_SUCCESS;
}


int DumpMansoTrackStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONMansoTrackStruct* track,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(track->fId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Track ID: " << track->fId << "\t";

	result = CheckField(track->fTrigRec, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Trigger Record ID: " << track->fTrigRec << endl;
	
	result = CheckField(track->fFlags, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Flags: " << showbase << hex << track->fFlags << dec;
	
	// Print the individual trigger bits.
	AliHLTMUONParticleSign sign;
	bool hitset[4];
	AliHLTMUONUtils::UnpackMansoTrackFlags(track->fFlags, sign, hitset);
	cout << " [Sign: " << sign << ", Hits set on chambers: ";
	bool first = true;
	for (AliHLTUInt32_t i = 0; i < 4; i++)
	{
		if (hitset[i])
		{
			cout << (first ? "" : ", ") << i+7;
			first = false;
		}
	}
	cout << (first ? "none]" : "]") << endl;

	result = CheckField(track->fPx, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Momentum: (px = " << track->fPx << ", ";

	result = CheckField(track->fPy, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "py = " << track->fPy << ", ";

	result = CheckField(track->fPz, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "pz = " << track->fPz << ") GeV/c\t";

	result = CheckField(track->fChi2, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Chi squared fit: " << track->fChi2 << endl;
	
	cout << "Track hits:" << endl;
	cout << "Chamber | X (cm)     | Y (cm)     | Z (cm)" << endl;
	cout << "------------------------------------------------" << endl;
	const AliHLTMUONRecHitStruct* hit = &track->fHit[0];
	for(AliHLTUInt32_t ch = 0; ch < 4; ch++)
	{
	        cout << setw(10) << left << ch + 7 << setw(0);
		result = DumpRecHitStruct(buffer, bufferSize, hit++, continueParse);
		if (result != EXIT_SUCCESS) return result;
	}

	return result;
}


int DumpMansoTracksBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	int result = EXIT_SUCCESS;
	AliHLTMUONMansoTracksBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	const AliHLTMUONMansoTrackStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		cout << "============================== Manso track number " << i+1
			<< " of " << nentries << " ==============================" << endl;
		int subResult = DumpMansoTrackStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpMansoRoIStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONMansoRoIStruct* roi,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(roi->fX, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << roi->fX << setw(0);

	result = CheckField(roi->fY, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << roi->fY << setw(0);

	result = CheckField(roi->fZ, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << roi->fZ << setw(0);

	result = CheckField(roi->fRadius, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << roi->fRadius << setw(0) << endl;

	return result;
}


int DumpMansoCandidateStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONMansoCandidateStruct* candidate,
		bool continueParse
	)
{
	int result = DumpMansoTrackStruct(buffer, bufferSize, &candidate->fTrack, continueParse);
	if (result != EXIT_SUCCESS) return result;
	
	cout << "Regions of interest:" << endl;
	cout << "Chamber | X (cm)     | Y (cm)     | Z (cm)     | Radius (cm)" << endl;
	cout << "-------------------------------------------------------------" << endl;
	const AliHLTMUONMansoRoIStruct* roi = &candidate->fRoI[0];
	for(AliHLTUInt32_t ch = 0; ch < 4; ch++)
	{
		cout << setw(10) << ch + 7;
		result = DumpMansoRoIStruct(buffer, bufferSize, roi++, continueParse);
		if (result != EXIT_SUCCESS) return result;
	}
	return result;
}


int DumpMansoCandidatesBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	int result = EXIT_SUCCESS;
	AliHLTMUONMansoCandidatesBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	const AliHLTMUONMansoCandidateStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		cout << "=========================== Manso track candidate number " << i+1
			<< " of " << nentries << " ===========================" << endl;
		int subResult = DumpMansoCandidateStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpSinglesDecisionBlockHeader(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONSinglesDecisionBlockStruct* header,
		bool continueParse
	)
{
	// Step through the header fields trying to print them.
	// At each step check if we have not overflowed the buffer, if we have
	// not then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckHeaderField(header->fNlowPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << " Number of low pt triggers: " << header->fNlowPt << endl;
	
	result = CheckHeaderField(header->fNhighPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Number of high pt triggers: " << header->fNhighPt << endl;

	return result;
}


int DumpTrackDecisionStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONTrackDecisionStruct* decision,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(decision->fTrackId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << decision->fTrackId << setw(0);
	
	result = CheckField(decision->fTriggerBits, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(12) << left << showbase << hex << decision->fTriggerBits
		<< setw(0) << dec;
		
	// Print the individual trigger bits.
	bool highPt, lowPt;
	AliHLTMUONUtils::UnpackTrackDecisionBits(decision->fTriggerBits, highPt, lowPt);
	cout << setw(7) << left << (highPt ? "yes" : "no");
	cout << setw(8) << left << (lowPt ? "yes" : "no");
	cout << setw(0) << endl;

	return result;
}


int DumpSinglesDecisionBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	int result = EXIT_SUCCESS;
	AliHLTMUONSinglesDecisionBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	// Dump the rest of the block header.
	const AliHLTMUONSinglesDecisionBlockStruct* header = &block.BlockHeader();
	int subResult = DumpSinglesDecisionBlockHeader(buffer, bufferSize, header, continueParse);
	if (subResult != EXIT_SUCCESS) return subResult;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	cout << "           |        Trigger Bits" << endl;
	cout << "Track ID   | Raw         HighPt LowPt" << endl;
	cout << "--------------------------------------" << endl;
	const AliHLTMUONTrackDecisionStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		subResult = DumpTrackDecisionStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpPairsDecisionBlockHeader(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONPairsDecisionBlockStruct* header,
		bool continueParse
	)
{
	// Step through the header fields trying to print them.
	// At each step check if we have not overflowed the buffer, if we have
	// not then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckHeaderField(header->fNunlikeAnyPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "      Number of unlike all pt triggers: " << header->fNunlikeAnyPt << endl;
	
	result = CheckHeaderField(header->fNunlikeLowPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "      Number of unlike low pt triggers: " << header->fNunlikeLowPt << endl;
	
	result = CheckHeaderField(header->fNunlikeHighPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "     Number of unlike high pt triggers: " << header->fNunlikeHighPt << endl;
	
	result = CheckHeaderField(header->fNlikeAnyPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "        Number of like any pt triggers: " << header->fNlikeAnyPt << endl;
	
	result = CheckHeaderField(header->fNlikeLowPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "        Number of like low pt triggers: " << header->fNlikeLowPt << endl;
	
	result = CheckHeaderField(header->fNlikeHighPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "       Number of like high pt triggers: " << header->fNlikeHighPt << endl;
	
	result = CheckHeaderField(header->fNmassAny, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << " Number of all invariant mass triggers: " << header->fNmassAny << endl;
	
	result = CheckHeaderField(header->fNmassLow, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << " Number of low invariant mass triggers: " << header->fNmassLow << endl;
	
	result = CheckHeaderField(header->fNmassHigh, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Number of high invariant mass triggers: " << header->fNmassHigh << endl;
	
	return result;
}


int DumpPairDecisionStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONPairDecisionStruct* decision,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(decision->fTrackAId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << decision->fTrackAId << setw(0);
	
	result = CheckField(decision->fTrackBId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << decision->fTrackBId << setw(0);
	
	result = CheckField(decision->fTriggerBits, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(12) << left << showbase << hex << decision->fTriggerBits
		<< setw(0) << dec;
		
	// Print the individual trigger bits.
	bool highMass, lowMass, unlike;
	AliHLTUInt8_t highPtCount, lowPtCount;
	AliHLTMUONUtils::UnpackPairDecisionBits(
			decision->fTriggerBits,
			highMass, lowMass, unlike, highPtCount, lowPtCount
		);
	cout << setw(7) << left << (highMass ? "yes" : "no");
	cout << setw(7) << left << (lowMass ? "yes" : "no");
	cout << setw(7) << left << (unlike ? "yes" : "no");
	cout << setw(6) << left << AliHLTUInt16_t(highPtCount);
	cout << setw(8) << left << AliHLTUInt16_t(lowPtCount);
	cout << setw(0);
	
	result = CheckField(decision->fInvMass, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << decision->fInvMass << endl;
	
	return EXIT_SUCCESS;
}


int DumpPairsDecisionBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	int result = EXIT_SUCCESS;
	AliHLTMUONPairsDecisionBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	// Dump the rest of the block header.
	const AliHLTMUONPairsDecisionBlockStruct* header = &block.BlockHeader();
	int subResult = DumpPairsDecisionBlockHeader(buffer, bufferSize, header, continueParse);
	if (subResult != EXIT_SUCCESS) return subResult;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	cout << "           |            |                Trigger Bits                  |" << endl;
	cout << "Track A ID | Track B ID | Raw         HiMass LoMass Unlike HiPt# LoPt# | Invariant mass" << endl;
	cout << "----------------------------------------------------------------------------------------" << endl;
	const AliHLTMUONPairDecisionStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		subResult = DumpPairDecisionStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpCommonHeader(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONDataBlockHeader* header, bool continueParse
	)
{
	// Step through the header fields trying to print them.
	// At each step check if we have not overflowed the buffer, if we have
	// not then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckHeaderField(header->fType, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	AliHLTMUONDataBlockType type = AliHLTMUONDataBlockType(header->fType);
	cout << "       Block type: " << type << endl;
	
	result = CheckHeaderField(header->fRecordWidth, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "     Record width: " << header->fRecordWidth << endl;
	
	result = CheckHeaderField(header->fNrecords, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Number of entries: " << header->fNrecords << endl;
	
	return result;
}


int ParseBuffer(
		const char* buffer, unsigned long bufferSize,
		bool continueParse, AliHLTMUONDataBlockType type
	)
{
	assert( buffer != NULL );
	int result = EXIT_SUCCESS;
	
	if (bufferSize < sizeof(AliHLTMUONDataBlockHeader))
	{
		cerr << "ERROR: The size of the file is too small to contain a"
			" valid data block." << endl;
		result = PARSE_ERROR;
		if (not continueParse) return result;
	}
	const AliHLTMUONDataBlockHeader* header =
		reinterpret_cast<const AliHLTMUONDataBlockHeader*>(buffer);

	int subResult = DumpCommonHeader(buffer, bufferSize, header, continueParse);
	if (subResult != EXIT_SUCCESS) return subResult;

	
	// Check if the block type in the header corresponds to the type given
	// by the '-type' command line parameter. If they do not then print an
	// error or big fat warning message and force interpretation of the data
	// block with the type given by '-type'.
	AliHLTMUONDataBlockType headerType = AliHLTMUONDataBlockType(header->fType);
	
	if (type == kUnknownDataBlock)
	{
		// -type not used in the command line so just use what is given
		// by the data block header.
		type = headerType;
	}
	else if (type != headerType)
	{
		cerr << "WARNING: The data block header indicates a type"
			" different from what was specified on the command line."
			" The data could be corrupt."
			<< endl;
		cerr << "WARNING: The type value in the file is "
			<< showbase << hex << header->fType
			<< " (" << headerType << "), but on the command line it is "
			<< showbase << hex << int(type) << dec
			<< " (" << type << ")."
			<< endl;
		cerr << "WARNING: Will force the interpretation of the data block"
			" with a type of " << type << "." << endl;
	}
	
	// Now we know what type the data block is supposed to be so we can
	// dump it to screen with the appropriate dump routine.
	switch (type)
	{
	case kTriggerRecordsDataBlock:
		subResult = DumpTriggerRecordsBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kTrigRecsDebugDataBlock:
		subResult = DumpTrigRecsDebugBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kTriggerChannelsDataBlock:
		subResult = DumpTriggerChannelsBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kRecHitsDataBlock:
		subResult = DumpRecHitsBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kClustersDataBlock:
		subResult = DumpClustersBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kChannelsDataBlock:
		return DumpChannelsBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kMansoTracksDataBlock:
		subResult = DumpMansoTracksBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kMansoCandidatesDataBlock:
		subResult = DumpMansoCandidatesBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kSinglesDecisionDataBlock:
		subResult = DumpSinglesDecisionBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kPairsDecisionDataBlock:
		return DumpPairsDecisionBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	default :
		cout << "ERROR: Unknown data block type. Found a type number of "
			<< showbase << hex << int(type) << dec
			<< " (" << int(type) << ")." << endl;
		result = PARSE_ERROR;
	}
	
	return result;
}


/**
 * The caller is responsible for freeing memory allocated for buffer with a call
 * to delete [] buffer.
 */
int ReadFile(const char* filename, char*& buffer, unsigned long& bufferSize)
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
	file.seekg(0, ios::end);
	if (not file)
	{
		cerr << "ERROR: Could not seek in the file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	bufferSize = file.tellg();
	if (not file)
	{
		cerr << "ERROR: Could not get file size for the file: " <<
			filename << endl;
		return SYSTEM_ERROR;
	}
	file.seekg(0, ios::beg);
	if (not file)
	{
		cerr << "ERROR: Could not seek in the file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	
	// Allocate the memory for the file.
	try
	{
		buffer = new char[bufferSize];
	}
	catch (const std::bad_alloc&)
	{
		cerr << "ERROR: Out of memory. Tried to allocate " << bufferSize
			<< " bytes." << endl;
		return SYSTEM_ERROR;
	}
	
	file.read(buffer, bufferSize);
	if (not file)
	{
		delete [] buffer;
		buffer = NULL;
		bufferSize = 0;
		cerr << "ERROR: Could not read from file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	
	file.close();
	if (not file)
	{
		delete [] buffer;
		buffer = NULL;
		bufferSize = 0;
		cerr << "ERROR: Could not close the file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	
	return EXIT_SUCCESS;
}

/**
 * Prints the command line usage of this program to standard error.
 */
void PrintUsage()
{
	cerr << "Usage: dHLTdumpraw [-help|-h] [-continue] [-type <typename>] <filename>" << endl;
	cerr << "Where <filename> is the name of a file containing a raw data block." << endl;
	cerr << "Options:" << endl;
	cerr << " -help | -h" << endl;
	cerr << "       Displays this message." << endl;
	cerr << " -continue" << endl;
	cerr << "       If specified, the program will try to continue parsing the data block" << endl;
	cerr << "       as much as possible rather than stopping at the first error." << endl;
	cerr << " -type <typename>" << endl;
	cerr << "       Forces the contents of the file to be interpreted as a specific" << endl;
	cerr << "       type of data block. Where <typename> can be one of:" << endl;
	cerr << "         trigrecs - trigger records data." << endl;
	cerr << "         trigrecsdebug - debugging information about trigger records." << endl;
	cerr << "         trigchannels - channel debugging in." << endl;
	cerr << "         rechits - reconstructed hits data." << endl;
	cerr << "         channels - channel debugging information from hit reconstruction." << endl;
	cerr << "         clusters - cluster debugging information from hit reconstruction." << endl;
	cerr << "         mansotracks - partial tracks from Manso algorithm." << endl;
	cerr << "         mansocandidates - track candidates considered in the Manso algorithm." << endl;
	cerr << "         singlesdecision - trigger decisions for single tracks." << endl;
	cerr << "         pairsdecision - trigger decisions for track pairs." << endl;
}

/**
 * Parse the string passed as the type of the block and return the corresponding
 * AliHLTMUONDataBlockType value.
 */
AliHLTMUONDataBlockType ParseCommandLineType(const char* type)
{
	if (strcmp(type, "trigrecs") == 0)
	{
		return kTriggerRecordsDataBlock;
	}
	else if (strcmp(type, "trigrecsdebug") == 0)
	{
		return kTrigRecsDebugDataBlock;
	}
	else if (strcmp(type, "trigchannels") == 0)
	{
		return kTriggerChannelsDataBlock;
	}
	else if (strcmp(type, "rechits") == 0)
	{      
		return kRecHitsDataBlock;
	}
	else if (strcmp(type,"channels") == 0)
	{
		return kChannelsDataBlock;
	}
	else if (strcmp(type,"clusters") == 0)
	{
		return kClustersDataBlock;
	}
	else if (strcmp(type, "mansotracks") == 0)
	{
		return kMansoTracksDataBlock;
	}
	else if (strcmp(type, "mansocandidates") == 0)
	{
		return kMansoCandidatesDataBlock;
	}
	else if (strcmp(type, "singlesdecision") == 0)
	{
		return kSinglesDecisionDataBlock;
	}
	else if (strcmp(type, "pairsdecision") == 0)
	{
		return kPairsDecisionDataBlock;
	}
	
	cerr << "ERROR: Invalid type name '" << type << "' specified for argument -type."
		<< endl << endl;
	PrintUsage();
	return kUnknownDataBlock;
}

/**
 * Parses the command line.
 * @param argc  Number of arguments as given in main().
 * @param argv  Array of arguments as given in main().
 * @param filename  Receives the pointer to the file name string.
 * @param type      Receives the type of the data block expected, i.e. the
 *                  value of the -type flag.
 * @return  A status flag suitable for returning from main(), containing either
 *          EXIT_SUCCESS or CMDLINE_ERROR.
 */
int ParseCommandLine(
		int argc, const char** argv,
		const char*& filename, bool& continueParse,
		AliHLTMUONDataBlockType& type
	)
{
	filename = NULL;
	continueParse = false;

	// Parse the command line.
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "-h") == 0)
		{
			PrintUsage();
			return EXIT_SUCCESS;
		}
		else if (strcmp(argv[i], "-continue") == 0)
		{
			continueParse = true;
		}
		else if (strcmp(argv[i], "-type") == 0)
		{
			// Now we need to parse the typename in the command line.
			type = ParseCommandLineType(argv[++i]);
			if (type == kUnknownDataBlock) return CMDLINE_ERROR;
		}
		else
		{
			if (filename != NULL)
			{
				cerr << "ERROR: Only one file can be specified, but got '"
					<< argv[i] << "', with '" << filename
					<< "' specified earlier." << endl << endl;
				PrintUsage();
				return CMDLINE_ERROR;
			}
			else
				filename = argv[i];
		}
	}
	
	// Now check that we have the filename and all the flags we need.
	if (filename == NULL)
	{
		cerr << "ERROR: Missing a file name. You must specify a file to process."
			<< endl << endl;
		PrintUsage();
		return CMDLINE_ERROR;
	}
	
	return EXIT_SUCCESS;
}


int main(int argc, const char** argv)
{
	const char* filename = NULL;
	bool continueParse = false;
	int returnCode = EXIT_SUCCESS;
	AliHLTMUONDataBlockType type = kUnknownDataBlock;
	char* buffer = NULL;

	try
	{
		returnCode = ParseCommandLine(argc, argv, filename, continueParse, type);

		if (returnCode == EXIT_SUCCESS and filename != NULL)
		{
			unsigned long bufferSize = 0;
			returnCode = ReadFile(filename, buffer, bufferSize);
			if (returnCode == EXIT_SUCCESS)
				returnCode = ParseBuffer(buffer, bufferSize, continueParse, type);
			if (buffer != NULL) delete [] buffer;
		}
		
	}
	catch (...)
	{
		cerr << "FATAL ERROR: An unknown exception occurred!" << endl << endl;
		returnCode = FATAL_ERROR;
		if (buffer != NULL) delete [] buffer;
	}

	return returnCode;
}
