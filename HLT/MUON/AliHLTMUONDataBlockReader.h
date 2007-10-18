#ifndef ALIHLTMUONDATABLOCKREADER_H
#define ALIHLTMUONDATABLOCKREADER_H
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
 * @file   AliHLTMUONDataBlockReader.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Definition of a reader class for internal dimuon HLT raw data blocks.
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
 * A light weight class for reading the contents of an internal dimuon HLT
 * data block.
 * Suppose we are given a pointer 'buffer' to the buffer where a data block is
 * stored in memory and the size of the data block is given by the variable 'size'.
 * The data block is of type 'block_type' and the data block entries are of type
 * 'entries_type'. The data block can be accessed in the following way:
 *
 *   void* buffer = somebuffer;
 *   AliHLTUInt32_t size = somebuffer_size;
 *   
 *   // Initialise the data block reader.
 *   AliHLTMUONDataBlockReader<block_type, entries_type> block(buffer, size);
 *   
 *   // Check that the buffer has the expected size.
 *   if (not block.BufferSizeOk())
 *   {
 *      // handle error...
 *   }
 *   
 *   // Find the number of entries in the data block.
 *   AliHLTUInt32_t nentries = block.Nentries();
 *   
 *   // Loop over all entries in the data block.
 *   for (AliHLTUInt32_t i = 0; i < nentries; i++)
 *   {
 *      const entries_type& entry = block[i];
 *      // Do something with the entry...
 *   }
 */
template <class DataBlockType, class DataElementType>
class AliHLTMUONDataBlockReader
{
public:
	typedef DataBlockType HeaderType;
	typedef DataElementType ElementType;

	/**
	 * Constructor that sets the internal pointer to the start of the data
	 * block and the total size of the block in bytes.
	 * @param buffer  The pointer to the first byte of the block in memory.
	 * @param size    The total size of the data block in bytes.
	 */
	AliHLTMUONDataBlockReader(const void* buffer, AliHLTUInt32_t size) :
		fSize(size),
		fBlock(reinterpret_cast<const DataBlockType*>(buffer)),
		fData(reinterpret_cast<const DataElementType*>(
		       reinterpret_cast<const DataBlockType*>(buffer) + 1
		      ))
	{
		assert( buffer != NULL );
	}
	
	/**
	 * Copy constructor that performs a shallow copy.
	 * Since this class does not take direct ownership of the buffer, never
	 * allocates or deallocates memory, this can be allowed.
	 */
	AliHLTMUONDataBlockReader(const AliHLTMUONDataBlockReader& reader)
	{
		fSize = reader.fSize;
		fBlock = reader.fBlock;
		fData = reader.fData;
	}
	
	/**
	 * Assignment operator performs a shallow copy.
	 * This is OK because this class does not take direct ownership of the
	 * output memory buffer.
	 */
	AliHLTMUONDataBlockReader& operator = (const AliHLTMUONDataBlockReader& reader)
	{
		fSize = reader.fSize;
		fBlock = reader.fBlock;
		fData = reader.fData;
		return *this;
	}

	/**
	 * Checks that the size of the buffer storing the data block is correct.
	 * Basic sanity checks are performed such as seeing if the data block
	 * size corresponds to the number of reconstructed hits stored and that
	 * the size of the buffer is at least sizeof(DataBlockType) bytes big.
	 */
	bool BufferSizeOk() const
	{
		// The block size must be at least sizeof(DataBlockType) bytes.
		// Do not try read the header otherwise, because we could get a
		// seg fault.
		if (fSize < sizeof(DataBlockType)) return false;

		// Now check if the size of the data block corresponds to the
		// number of entries it claims to contain.
		AliHLTUInt32_t arraysize = fSize - sizeof(DataBlockType);
		return arraysize == Nentries() * sizeof(DataElementType);
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
	const DataBlockType& BlockHeader() const
	{
		return *fBlock;
	}

	/**
	 * Returns the total number of entries in the data block.
	 */
	AliHLTUInt32_t Nentries() const { return fBlock->fHeader.fNrecords; }

	/**
	 * Returns a pointer to the i'th entry.
	 * If the index 'i' is out of bounds then NULL is returned.
	 * This is a safe access method because it does bounds checking but is
	 * a little slower than the array operator.
	 * @param i  The index number of the entry to be returned.
	 * @return  A pointer to the entry or NULL.
	 */
	const DataElementType* Entry(AliHLTUInt32_t i) const
	{
		return (i < Nentries()) ? &fData[i] : NULL;
	}

	/**
	 * Array operator for accessing the data entries directly.
	 * The index variable 'i' is not checked (except in debug compilations)
	 * so one should make sure they are within the valid range.
	 */
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
	const DataElementType* GetArray() const { return fData; }

	/**
	 * Calculates the number of bytes used for the data block in the buffer.
	 * This value should be the same as what is returned by BufferSize()
	 * unless too much buffer space was allocated.
	 */
	AliHLTUInt32_t BytesUsed() const
	{
		assert( sizeof(DataElementType)	== fBlock->fHeader.fRecordWidth);
		return sizeof(DataBlockType) + Nentries() * sizeof(DataElementType);
	}

	AliHLTUInt32_t BufferSize() { return fSize; }
	
private:

	AliHLTUInt32_t fSize;   // Size of the data block in bytes.
	const DataBlockType* fBlock; // Pointer to the data block buffer.
	const DataElementType* fData; // Pointer to the data array.
};


// We now define the reader classes for the various data block types from the
// template class AliHLTMUONDataBlockReader.

typedef AliHLTMUONDataBlockReader<
		AliHLTMUONTriggerRecordsBlockStruct,
		AliHLTMUONTriggerRecordStruct
	> AliHLTMUONTriggerRecordsBlockReader;

typedef AliHLTMUONDataBlockReader<
		AliHLTMUONTrigRecsDebugBlockStruct,
		AliHLTMUONTrigRecInfoStruct
	> AliHLTMUONTrigRecsDebugBlockReader;

typedef AliHLTMUONDataBlockReader<
		AliHLTMUONTriggerChannelsBlockStruct,
		AliHLTMUONTriggerChannelStruct
	> AliHLTMUONTriggerChannelsBlockReader;

typedef AliHLTMUONDataBlockReader<
		AliHLTMUONRecHitsBlockStruct,
		AliHLTMUONRecHitStruct
	> AliHLTMUONRecHitsBlockReader;

typedef AliHLTMUONDataBlockReader<
		AliHLTMUONClustersBlockStruct,
		AliHLTMUONClusterStruct
	> AliHLTMUONClustersBlockReader;

typedef AliHLTMUONDataBlockReader<
		AliHLTMUONChannelsBlockStruct,
		AliHLTMUONChannelStruct
	> AliHLTMUONChannelsBlockReader;

typedef AliHLTMUONDataBlockReader<
		AliHLTMUONMansoTracksBlockStruct,
		AliHLTMUONMansoTrackStruct
	> AliHLTMUONMansoTracksBlockReader;

typedef AliHLTMUONDataBlockReader<
		AliHLTMUONMansoCandidatesBlockStruct,
		AliHLTMUONMansoCandidateStruct
	> AliHLTMUONMansoCandidatesBlockReader;

typedef AliHLTMUONDataBlockReader<
		AliHLTMUONSinglesDecisionBlockStruct,
		AliHLTMUONTrackDecisionStruct
	> AliHLTMUONSinglesDecisionBlockReader;

typedef AliHLTMUONDataBlockReader<
		AliHLTMUONPairsDecisionBlockStruct,
		AliHLTMUONPairDecisionStruct
	> AliHLTMUONPairsDecisionBlockReader;

#endif // ALIHLTMUONDATABLOCKREADER_H
