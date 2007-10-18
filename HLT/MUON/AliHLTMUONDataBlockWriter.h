#ifndef ALIHLTMUONDATABLOCKWRITER_H
#define ALIHLTMUONDATABLOCKWRITER_H
/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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
 * @file   AliHLTMUONDataBlockWriter.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Definition of a writer class for internal dimuon HLT raw data blocks.
 */

#include "AliHLTMUONDataTypes.h"
#include <cassert>

#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include "AliHLTMUONTrigRecsDebugBlockStruct.h"
#include "AliHLTMUONTriggerChannelsBlockStruct.h"
#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONClustersBlockStruct.h"
#include "AliHLTMUONChannelsBlockStruct.h"
#include "AliHLTMUONMansoTracksBlockStruct.h"
#include "AliHLTMUONMansoCandidatesBlockStruct.h"
#include "AliHLTMUONSinglesDecisionBlockStruct.h"
#include "AliHLTMUONPairsDecisionBlockStruct.h"

/**
 * A light weight class for writing an internal dimuon HLT data block.
 * Suppose we are given a pointer 'buffer' to some empty memory buffer where we
 * can store our new data block and the size of the data block is given by the
 * variable 'size'. The data block is of type 'block_type', the data block entries
 * are of type 'entries_type' and the data block type code is 'type_code'.
 * The data block can be written in the following way:
 *
 *   void* buffer = somebuffer;
 *   AliHLTUInt32_t size = somebuffer_size;
 *   
 *   // Initialise the data block writer.
 *   AliHLTMUONDataBlockWriter<block_type, entries_type, type_code>
 *   block(buffer, size);
 *   
 *   // Initialise the block header
 *   if (not block.InitCommonHeader())
 *   {
 *      // handle error and exit...
 *   }
 *   
 *   // Tell the writer how many entries we are going to use.
 *   if (not block.SetNumberOfEntries(somevalue))
 *   {
 *      // handle error and exit...
 *   }
 *
 *   // Add all the entries to the data block.
 *   for (AliHLTUInt32_t i = 0; i < block.Nentries(); i++)
 *   {
 *      entries_type& entry = block[i];
 *      // fill the new entry...
 *      entry.somefield = somevalue;
 *   }
 *
 * The slightly slower but safer method is to do the following:
 *
 *   AliHLTMUONDataBlockWriter<block_type, entries_type, type_code>
 *   block(buffer, size);
 *   if (not block.InitCommonHeader())
 *   {
 *      // handle error and exit...
 *   }
 *   
 *   // For each new entry add it to the data block.
 *   while (HaveMoreEntries())
 *   {
 *      entries_type* entry = block.AddEntry();
 *      if (entry == NULL)
 *      {
 *         // handle buffer overflow and exit...
 *      }
 *      // fill the new entry...
 *      entry->somefield = somevalue;
 *   }
 */
template <
	class DataBlockType,
	class DataElementType,
	AliHLTMUONDataBlockType blockTypeCode
>
class AliHLTMUONDataBlockWriter
{
public:
	typedef DataBlockType HeaderType;
	typedef DataElementType ElementType;

	/**
	 * Constructor that sets the internal pointer to the start of the buffer
	 * space to write to and the total size of the buffer in bytes.
	 * @param buffer  The pointer to the first byte of the memory buffer.
	 * @param size    The total size of the buffer in bytes.
	 */
	AliHLTMUONDataBlockWriter(void* buffer, AliHLTUInt32_t size) :
		fSize(size),
		fMaxArraySize(size > sizeof(DataBlockType) ? size - sizeof(DataBlockType) : 0),
		fBlock(reinterpret_cast<DataBlockType*>(buffer)),
		fData(reinterpret_cast<DataElementType*>(
		       reinterpret_cast<DataBlockType*>(buffer) + 1
		      ))
	{
		assert( buffer != NULL );
	}
	
	/**
	 * Copy constructor that performs a shallow copy.
	 * Since this class does not take direct ownership of the buffer, never
	 * allocates or deallocates memory, this can be allowed.
	 */
	AliHLTMUONDataBlockWriter(const AliHLTMUONDataBlockWriter& writer)
	{
		fSize = writer.fSize;
		fMaxArraySize = writer.fMaxArraySize;
		fBlock = writer.fBlock;
		fData = writer.fData;
	}
	
	/**
	 * Assignment operator performs a shallow copy.
	 * This is OK because this class does not take direct ownership of the
	 * output memory buffer.
	 */
	AliHLTMUONDataBlockWriter& operator = (const AliHLTMUONDataBlockWriter& writer)
	{
		fSize = writer.fSize;
		fMaxArraySize = writer.fMaxArraySize;
		fBlock = writer.fBlock;
		fData = writer.fData;
		return *this;
	}

	/**
	 * Initialises the common data block header by setting the type and record
	 * width fields. If the buffer size was to small to create the header then
	 * this method returns false, otherwise true on success.
	 */
	bool InitCommonHeader() const
	{
		// The block size must be at least sizeof(DataBlockType) bytes.
		if (fSize < sizeof(DataBlockType)) return false;

		// Now fill the fields in the header.
		fBlock->fHeader.fType = blockTypeCode;
		fBlock->fHeader.fRecordWidth = sizeof(DataElementType);
		fBlock->fHeader.fNrecords = 0;
		return true;
	}
	
	/**
	 * Returns the common data block header.
	 */
	const AliHLTMUONDataBlockHeader& CommonBlockHeader() const
	{
		return fBlock->fHeader;
	}
	
	/**
	 * Returns the whole data block header.
	 */
	DataBlockType& BlockHeader()
	{
		return *fBlock;
	}
	
	const DataBlockType& BlockHeader() const
	{
		return *fBlock;
	}
	
	/**
	 * Returns a pointer to the next location where a data entry can be
	 * written and increments the number of entries.
	 * If the buffer is already full then NULL is returned and the number of
	 * entries is not changed.
	 */
	DataElementType* AddEntry() const
	{
		if ( (Nentries() + 1) * sizeof(DataElementType) > fMaxArraySize )
			return NULL;
		DataElementType* newentry = &fData[fBlock->fHeader.fNrecords];
		fBlock->fHeader.fNrecords++;
		return newentry;
	}
	
	/**
	 * Sets the number of entries that will be filled into the buffer.
	 * If the number of entries is to many to fit into the buffer then this
	 * method returns false, otherwise true.
	 */
	bool SetNumberOfEntries(AliHLTUInt32_t n) const
	{
		if (n * sizeof(DataElementType) > fMaxArraySize) return false;
		fBlock->fHeader.fNrecords = n;
		return true;
	}

	/**
	 * Returns the total number of entries already added to the data block.
	 */
	AliHLTUInt32_t Nentries() const { return fBlock->fHeader.fNrecords; }

	/**
	 * Returns a pointer to the i'th entry.
	 * If the index 'i' is out of bounds then NULL is returned.
	 * This is a safe access method because it does bounds checking but is
	 * a little slower than the array operators.
	 * @param i  The index number of the entry to be returned.
	 * @return  A pointer to the entry or NULL.
	 */
	DataElementType* Entry(AliHLTUInt32_t i)
	{
		return (i < Nentries()) ? &fData[i] : NULL;
	}
	
	const DataElementType* Entry(AliHLTUInt32_t i) const
	{
		return (i < Nentries()) ? &fData[i] : NULL;
	}

	/**
	 * Array operator for accessing the data entries directly.
	 * The index variable 'i' is not checked (except in debug compilations)
	 * so one should make sure they are within the valid range.
	 */
	DataElementType& operator [] (AliHLTUInt32_t i)
	{
		assert( i < Nentries() );
		return fData[i];
	}
	
	const DataElementType& operator [] (AliHLTUInt32_t i) const
	{
		assert( i < Nentries() );
		return fData[i];
	}

	/**
	 * Returns a pointer to the array of elements in the data block.
	 * Care must be taken not to read beyond the array limits given by
	 * Nentries().
	 */
	DataElementType* GetArray() { return fData; }
	const DataElementType* GetArray() const { return fData; }
	
	/**
	 * Calculates the number of bytes used for the data block in the buffer.
	 * This value will only make sense if a call to InitCommonHeader() was
	 * made and it returned true.
	 */
	AliHLTUInt32_t BytesUsed() const
	{
		assert( sizeof(DataElementType) == fBlock->fHeader.fRecordWidth);
		return sizeof(DataBlockType) + Nentries() * sizeof(DataElementType);
	}
	
	/**
	 * Calculates the maximum number of entries that will fit into the
	 * memory buffer.
	 */
	AliHLTUInt32_t MaxNumberOfEntries() const
	{
		return fMaxArraySize / sizeof(DataElementType);
	}
	
	AliHLTUInt32_t BufferSize() { return fSize; }
	
private:

	AliHLTUInt32_t fSize;   // Size of the buffer in bytes.
	AliHLTUInt32_t fMaxArraySize; // Maximum size of the fData array in bytes.
	DataBlockType* fBlock; // Pointer to the start of the data block.
	DataElementType* fData; // Pointer to the start of the data array.
};


// We now define the writer classes for the various data block types from the
// template class AliHLTMUONDataBlockWriter.

typedef AliHLTMUONDataBlockWriter<
		AliHLTMUONTriggerRecordsBlockStruct,
		AliHLTMUONTriggerRecordStruct,
		kTriggerRecordsDataBlock
	> AliHLTMUONTriggerRecordsBlockWriter;

typedef AliHLTMUONDataBlockWriter<
		AliHLTMUONTrigRecsDebugBlockStruct,
		AliHLTMUONTrigRecInfoStruct,
		kTrigRecsDebugDataBlock
	> AliHLTMUONTrigRecsDebugBlockWriter;

typedef AliHLTMUONDataBlockWriter<
		AliHLTMUONTriggerChannelsBlockStruct,
		AliHLTMUONTriggerChannelStruct,
		kTriggerChannelsDataBlock
	> AliHLTMUONTriggerChannelsBlockWriter;

typedef AliHLTMUONDataBlockWriter<
		AliHLTMUONRecHitsBlockStruct,
		AliHLTMUONRecHitStruct,
		kRecHitsDataBlock
	> AliHLTMUONRecHitsBlockWriter;

typedef AliHLTMUONDataBlockWriter<
		AliHLTMUONClustersBlockStruct,
		AliHLTMUONClusterStruct,
		kClustersDataBlock
	> AliHLTMUONClustersBlockWriter;

typedef AliHLTMUONDataBlockWriter<
		AliHLTMUONChannelsBlockStruct,
		AliHLTMUONChannelStruct,
		kChannelsDataBlock
	> AliHLTMUONChannelsBlockWriter;

typedef AliHLTMUONDataBlockWriter<
		AliHLTMUONMansoTracksBlockStruct,
		AliHLTMUONMansoTrackStruct,
		kMansoTracksDataBlock
	> AliHLTMUONMansoTracksBlockWriter;

typedef AliHLTMUONDataBlockWriter<
		AliHLTMUONMansoCandidatesBlockStruct,
		AliHLTMUONMansoCandidateStruct,
		kMansoCandidatesDataBlock
	> AliHLTMUONMansoCandidatesBlockWriter;
	
typedef AliHLTMUONDataBlockWriter<
		AliHLTMUONSinglesDecisionBlockStruct,
		AliHLTMUONTrackDecisionStruct,
		kSinglesDecisionDataBlock
	> AliHLTMUONSinglesDecisionBlockWriter;
	
typedef AliHLTMUONDataBlockWriter<
		AliHLTMUONPairsDecisionBlockStruct,
		AliHLTMUONPairDecisionStruct,
		kPairsDecisionDataBlock
	> AliHLTMUONPairsDecisionBlockWriter;

#endif // ALIHLTMUONDATABLOCKWRITER_H
