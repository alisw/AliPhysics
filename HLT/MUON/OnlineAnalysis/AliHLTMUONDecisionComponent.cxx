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

///
///  @file   AliHLTMUONDecisionComponent.cxx
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   30 April 2008
///  @brief  Implementation of the decision component for dimuon HLT triggering.
///
// class documentation is in the header file.

#include "AliHLTMUONDecisionComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONCalculations.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONDataBlockWriter.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TObjString.h"
#include "TString.h"
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cmath>
#include <new>


// Helper type for memory allocation.
typedef const AliHLTMUONMansoTrackStruct* AliHLTMUONMansoTrackStructP;


ClassImp(AliHLTMUONDecisionComponent);


AliHLTMUONDecisionComponent::AliHLTMUONDecisionComponent() :
	AliHLTProcessor(),
	fMaxTracks(1),
	fTrackCount(0),
	fTracks(new AliHLTMUONMansoTrackStructP[fMaxTracks]),
	fLowPtCut(1.),  // 1 GeV/c cut
	fHighPtCut(2.),  // 2 GeV/c cut
	fLowMassCut(2.5),  // 2.7 GeV/c^2 cut
	fHighMassCut(7.),  // 8 GeV/c^2 cut
	fWarnForUnexpecedBlock(false)
{
	///
	/// Default constructor.
	///
}


AliHLTMUONDecisionComponent::~AliHLTMUONDecisionComponent()
{
	///
	/// Default destructor deletes the fTracks array.
	///
	
	assert(fTracks != NULL);
	delete [] fTracks;
}


const char* AliHLTMUONDecisionComponent::GetComponentID()
{
	///
	/// Inherited from AliHLTComponent. Returns the component ID.
	///
	
	return AliHLTMUONConstants::DecisionComponentId();
}


void AliHLTMUONDecisionComponent::GetInputDataTypes(
		vector<AliHLTComponentDataType>& list
	)
{
	///
	/// Inherited from AliHLTProcessor. Returns the list of expected input data types.
	///
	
	assert( list.empty() );
	list.push_back( AliHLTMUONConstants::MansoTracksBlockDataType() );
}


AliHLTComponentDataType AliHLTMUONDecisionComponent::GetOutputDataType()
{
	/// Inherited from AliHLTComponent. Returns kAliHLTMultipleDataType
	/// refer to GetOutputDataTypes for all returned data types.
	
	return kAliHLTMultipleDataType;
}


int AliHLTMUONDecisionComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
	/// Inherited from AliHLTComponent. Returns the output data types.
	
	assert( list.empty() );
	list.push_back( AliHLTMUONConstants::SinglesDecisionBlockDataType() );
	list.push_back( AliHLTMUONConstants::PairsDecisionBlockDataType() );
	return 1;
}


void AliHLTMUONDecisionComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	///
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	///
	
	constBase = sizeof(AliHLTMUONSinglesDecisionBlockStruct);
	constBase += sizeof(AliHLTMUONPairsDecisionBlockStruct);
	inputMultiplier = 1;
}


AliHLTComponent* AliHLTMUONDecisionComponent::Spawn()
{
	///
	/// Inherited from AliHLTComponent. Creates a new object instance.
	///
	
	return new AliHLTMUONDecisionComponent;
}


int AliHLTMUONDecisionComponent::DoInit(int argc, const char** argv)
{
	///
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	///
	
	HLTInfo("Initialising dHLT trigger decision component.");
	
	bool lowPtCutSet = false;
	bool highPtCutSet = false;
	bool lowMassCutSet = false;
	bool highMassCutSet = false;
	fWarnForUnexpecedBlock = false;
	
	for (int i = 0; i < argc; i++)
	{
		if (strcmp( argv[i], "-lowptcut" ) == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The value for the low pT cut was not specified.");
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			double num = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a floating point number.", argv[i+1]);
				return -EINVAL;
			}
			fLowPtCut = (AliHLTFloat32_t)num;
			lowPtCutSet = true;
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-highptcut" ) == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The value for the high pT cut was not specified.");
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			double num = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a floating point number.", argv[i+1]);
				return -EINVAL;
			}
			fHighPtCut = (AliHLTFloat32_t)num;
			highPtCutSet = true;
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-lowmasscut" ) == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The value for the low invariant mass cut was not specified.");
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			double num = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a floating point number.", argv[i+1]);
				return -EINVAL;
			}
			fLowMassCut = (AliHLTFloat32_t)num;
			lowMassCutSet = true;
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-highmasscut" ) == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The value for the high invariant mass cut was not specified.");
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			double num = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a floating point number.", argv[i+1]);
				return -EINVAL;
			}
			fHighMassCut = (AliHLTFloat32_t)num;
			highMassCutSet = true;
			
			i++;
			continue;
		}
		
		if (strcmp(argv[i], "-warn_on_unexpected_block") == 0)
		{
			fWarnForUnexpecedBlock = true;
			continue;
		}

		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}
	
	// Read cut parameters from CDB if they were not specified on the command line.
	if (not lowPtCutSet or not highPtCutSet or not lowMassCutSet or not highMassCutSet)
	{
		int result = ReadConfigFromCDB(
				NULL,
				not lowPtCutSet, not highPtCutSet,
				not lowMassCutSet, not highMassCutSet
			);
		if (result != 0) return result;
	}
	
	HLTDebug("Using the following cut parameters:");
	HLTDebug("              Low pT cut = %f GeV/c", fLowPtCut);
	HLTDebug("             High pT cut = %f GeV/c", fHighPtCut);
	HLTDebug("  Low invariant mass cut = %f GeV/c^2", fLowMassCut);
	HLTDebug(" High invariant mass cut = %f GeV/c^2", fHighMassCut);
	
	return 0;
}


int AliHLTMUONDecisionComponent::DoDeinit()
{
	///
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	///
	
	HLTInfo("Deinitialising dHLT trigger decision component.");
	return 0;
}


int AliHLTMUONDecisionComponent::Reconfigure(const char* cdbEntry, const char* componentId)
{
	/// Inherited from AliHLTComponent. Reconfigures the component from CDB.
	
	if (strcmp(componentId, GetComponentID()) == 0)
	{
		HLTInfo("Reading new entries for cut parameters from CDB.");
		int result = ReadConfigFromCDB(cdbEntry);
		HLTDebug("Using the following new cut parameters:");
		HLTDebug("              Low pT cut = %f GeV/c", fLowPtCut);
		HLTDebug("             High pT cut = %f GeV/c", fHighPtCut);
		HLTDebug("  Low invariant mass cut = %f GeV/c^2", fLowMassCut);
		HLTDebug(" High invariant mass cut = %f GeV/c^2", fHighMassCut);
		return result;
	}
	else
		return 0;
}


int AliHLTMUONDecisionComponent::DoEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks,
		AliHLTComponentTriggerData& /*trigData*/,
		AliHLTUInt8_t* outputPtr,
		AliHLTUInt32_t& size,
		std::vector<AliHLTComponentBlockData>& outputBlocks
	)
{
	///
	/// Inherited from AliHLTProcessor. Processes the new event data.
	///
	
	AliHLTUInt32_t specification = 0;  // Contains the output data block spec bits.
	
	// Loop over all input blocks in the event with track data and add pointers
	// to the tracks into the tracks array. These will be used later by the
	// trigger algorithm to get to the individual tracks.
	fTrackCount = 0; // reset number of tracks in array.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
		if (blocks[n].fDataType == AliHLTMUONConstants::MansoTracksBlockDataType())
		{
			// Build up the specification which indicates what DDLs
			// contributed to the output data.
			specification |= blocks[n].fSpecification;
			
			AliHLTMUONMansoTracksBlockReader inblock(blocks[n].fPtr, blocks[n].fSize);
			if (not inblock.BufferSizeOk())
			{
				size_t headerSize = sizeof(AliHLTMUONMansoTracksBlockReader::HeaderType);
				if (blocks[n].fSize < headerSize)
				{
					HLTError("Received a manso tracks data block with a size of %d bytes,"
						" which is smaller than the minimum valid header size of %d bytes."
						" The block must be corrupt.",
						blocks[n].fSize, headerSize
					);
					continue;
				}
				
				size_t expectedWidth = sizeof(AliHLTMUONMansoTracksBlockReader::ElementType);
				if (inblock.CommonBlockHeader().fRecordWidth != expectedWidth)
				{
					HLTError("Received a manso tracks data block with a record"
						" width of %d bytes, but the expected value is %d bytes."
						" The block might be corrupt.",
						inblock.CommonBlockHeader().fRecordWidth, expectedWidth
					);
					continue;
				}
				
				HLTError("Received a manso tracks data block with a size of %d bytes,"
					" but the block header claims the block should be %d bytes."
					" The block might be corrupt.",
					blocks[n].fSize, inblock.BytesUsed()
				);
				continue;
			}
			
			for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
			{
				int result = AddTrack(&inblock[i]);
				if (result != 0)
				{
					size = 0; // Important to tell framework that nothing was generated.
					return result;
				}
			}
		}
		else if (blocks[n].fDataType != AliHLTMUONConstants::MansoTracksBlockDataType())
		{
			// Log a message indicating that we got a data block that we
			// do not know how to handle.
			char id[kAliHLTComponentDataTypefIDsize+1];
			for (int i = 0; i < kAliHLTComponentDataTypefIDsize; i++)
				id[i] = blocks[n].fDataType.fID[i];
			id[kAliHLTComponentDataTypefIDsize] = '\0';
			char origin[kAliHLTComponentDataTypefOriginSize+1];
			for (int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++)
				origin[i] = blocks[n].fDataType.fOrigin[i];
			origin[kAliHLTComponentDataTypefOriginSize] = '\0';
			
			if (fWarnForUnexpecedBlock)
				HLTWarning("Received a data block of a type we cannot handle: %s origin: %s",
					static_cast<char*>(id), static_cast<char*>(origin)
				);
			else
				HLTDebug("Received a data block of a type we cannot handle: %s origin: %s",
					static_cast<char*>(id), static_cast<char*>(origin)
				);
		}
	}
	
	// Now we can create our two new output data blocks for the single tracks
	// and track pairs.
	AliHLTMUONSinglesDecisionBlockWriter singlesBlock(outputPtr, size);
	
	if (not singlesBlock.InitCommonHeader())
	{
		Logging(kHLTLogError,
			"AliHLTMUONDecisionComponent::DoEvent",
			"Buffer overflow",
			"The buffer is only %d bytes in size. We need a minimum of"
			" %d bytes for the singles output data block.",
			size, sizeof(AliHLTMUONSinglesDecisionBlockWriter::HeaderType)
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -ENOBUFS;
	}
	
	if (not singlesBlock.SetNumberOfEntries(fTrackCount))
	{
		AliHLTUInt32_t bytesneeded = sizeof(AliHLTMUONSinglesDecisionBlockWriter::HeaderType)
			+ fTrackCount * sizeof(AliHLTMUONSinglesDecisionBlockWriter::ElementType);
		HLTError("The buffer is only %d bytes in size. We need a minimum of"
			" %d bytes for the singles output data block.",
			size, bytesneeded
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -ENOBUFS;
	}
	
	AliHLTMUONPairsDecisionBlockWriter pairsBlock(
			outputPtr + singlesBlock.BytesUsed(),
			size - singlesBlock.BytesUsed()
		);
	
	if (not pairsBlock.InitCommonHeader())
	{
		Logging(kHLTLogError,
			"AliHLTMUONDecisionComponent::DoEvent",
			"Buffer overflow",
			"The buffer is only %d bytes in size. We need a minimum of"
			" %d bytes for the pairs output data block.",
			size,
			sizeof(AliHLTMUONPairsDecisionBlockWriter::HeaderType) + singlesBlock.BytesUsed()
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -ENOBUFS;
	}
	
	AliHLTUInt32_t numOfPairs = fTrackCount * (fTrackCount-1) / 2;
	if (not pairsBlock.SetNumberOfEntries(numOfPairs))
	{
		AliHLTUInt32_t bytesneeded = sizeof(AliHLTMUONPairsDecisionBlockWriter::HeaderType)
			+ fTrackCount * sizeof(AliHLTMUONPairsDecisionBlockWriter::ElementType)
			+ singlesBlock.BytesUsed();
		HLTError("The buffer is only %d bytes in size. We need a minimum of"
			" %d bytes for the pairs output data block.",
			size, bytesneeded
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -ENOBUFS;
	}
	
	ApplyTriggerAlgorithm(
			singlesBlock.BlockHeader(),
			singlesBlock.GetArray(),
			pairsBlock.BlockHeader(),
			pairsBlock.GetArray()
		);
	
	AliHLTComponentBlockData sbd;
	FillBlockData(sbd);
	sbd.fPtr = outputPtr;
	sbd.fOffset = 0;
	sbd.fSize = singlesBlock.BytesUsed();
	sbd.fDataType = AliHLTMUONConstants::SinglesDecisionBlockDataType();
	sbd.fSpecification = specification;
	outputBlocks.push_back(sbd);
	size = singlesBlock.BytesUsed();
	
	AliHLTComponentBlockData pbd;
	FillBlockData(pbd);
	pbd.fPtr = outputPtr;
	pbd.fOffset = singlesBlock.BytesUsed();
	pbd.fSize = pairsBlock.BytesUsed();
	pbd.fDataType = AliHLTMUONConstants::PairsDecisionBlockDataType();
	pbd.fSpecification = specification;
	outputBlocks.push_back(pbd);
	size += pairsBlock.BytesUsed();
	
	return 0;
}


int AliHLTMUONDecisionComponent::ReadConfigFromCDB(
		const char* path,
		bool setLowPtCut, bool setHighPtCut,
		bool setLowMassCut, bool setHighMassCut
	)
{
	/// Reads the cut parameters from the CDB.
	
	assert(AliCDBManager::Instance() != NULL);
	
	const char* pathToEntry = AliHLTMUONConstants::DecisionComponentCDBPath();
	if (path != NULL)
		pathToEntry = path;
	
	AliCDBEntry* entry = AliCDBManager::Instance()->Get(pathToEntry);
	if (entry == NULL)
	{
		HLTError("Could not get the CDB entry for \"%s\".", pathToEntry);
		return -EIO;
	}
	
	TObject* obj = entry->GetObject();
	if (obj == NULL)
	{
		HLTError("Configuration object for \"%s\" is missing.", pathToEntry);
		return -ENOENT;
	}
	
	if (obj->IsA() != TMap::Class())
	{
		HLTError("Wrong type for configuration object in \"%s\". Found a %s but we need a TMap.",
			pathToEntry, obj->ClassName()
		);
		return -EPROTO;
	}
	TMap* map = dynamic_cast<TMap*>(obj);
	
	if (setLowPtCut)
	{
		TPair* pair = static_cast<TPair*>(map->FindObject("lowptcut"));
		if (pair == NULL)
		{
			HLTError("Configuration object for \"%s\" does not contain the low pT cut value.",
				pathToEntry
			);
			return -ENOENT;
		}
		TObject* valueObj = pair->Value();
		if (valueObj->IsA() != TObjString::Class())
		{
			HLTError("The low pT cut parameter found in configuration object \"%s\""
				" is not a TObjString. Found an object of type %s instead.",
				pathToEntry, valueObj->ClassName()
			);
			return -EPROTO;
		}
		TString value = dynamic_cast<TObjString*>(valueObj)->GetString();
		if (not value.IsFloat())
		{
			HLTError("The low pT cut parameter found in configuration object \"%s\""
				"is not a floating point string; found \"%s\".",
				pathToEntry, value.Data()
			);
			return -EPROTO;
		}
		fLowPtCut = (AliHLTFloat32_t) value.Atof();
	}
	
	if (setHighPtCut)
	{
		TPair* pair = static_cast<TPair*>(map->FindObject("highptcut"));
		if (pair == NULL)
		{
			HLTError("Configuration object for \"%s\" does not contain the high pT cut value.",
				pathToEntry
			);
			return -ENOENT;
		}
		TObject* valueObj = pair->Value();
		if (valueObj->IsA() != TObjString::Class())
		{
			HLTError("The high pT cut parameter found in configuration object \"%s\""
				" is not a TObjString. Found an object of type %s instead.",
				pathToEntry, valueObj->ClassName()
			);
			return -EPROTO;
		}
		TString value = dynamic_cast<TObjString*>(valueObj)->GetString();
		if (not value.IsFloat())
		{
			HLTError("The high pT cut parameter found in configuration object \"%s\""
				"is not a floating point string; found \"%s\".",
				pathToEntry, value.Data()
			);
			return -EPROTO;
		}
		fHighPtCut = (AliHLTFloat32_t) value.Atof();
	}
	
	if (setLowMassCut)
	{
		TPair* pair = static_cast<TPair*>(map->FindObject("lowmasscut"));
		if (pair == NULL)
		{
			HLTError("Configuration object for \"%s\" does not contain the low invariant mass cut value.",
				pathToEntry
			);
			return -ENOENT;
		}
		TObject* valueObj = pair->Value();
		if (valueObj->IsA() != TObjString::Class())
		{
			HLTError("The low invariant mass cut parameter found in configuration object \"%s\""
				" is not a TObjString. Found an object of type %s instead.",
				pathToEntry, valueObj->ClassName()
			);
			return -EPROTO;
		}
		TString value = dynamic_cast<TObjString*>(valueObj)->GetString();
		if (not value.IsFloat())
		{
			HLTError("The low invariant mass cut parameter found in configuration object \"%s\""
				"is not a floating point string; found \"%s\".",
				pathToEntry, value.Data()
			);
			return -EPROTO;
		}
		fLowMassCut = (AliHLTFloat32_t) value.Atof();
	}
	
	if (setHighMassCut)
	{
		TPair* pair = static_cast<TPair*>(map->FindObject("highmasscut"));
		if (pair == NULL)
		{
			HLTError("Configuration object for \"%s\" does not contain the high invariant mass cut value.",
				pathToEntry
			);
			return -ENOENT;
		}
		TObject* valueObj = pair->Value();
		if (valueObj->IsA() != TObjString::Class())
		{
			HLTError("The high invariant mass cut parameter found in configuration object \"%s\""
				" is not a TObjString. Found an object of type %s instead.",
				pathToEntry, valueObj->ClassName()
			);
			return -EPROTO;
		}
		TString value = dynamic_cast<TObjString*>(valueObj)->GetString();
		if (not value.IsFloat())
		{
			HLTError("The high invariant mass cut parameter found in configuration object \"%s\""
				"is not a floating point string; found \"%s\".",
				pathToEntry, value.Data()
			);
			return -EPROTO;
		}
		fHighMassCut = (AliHLTFloat32_t) value.Atof();
	}
	
	return 0;
}


int AliHLTMUONDecisionComponent::AddTrack(const AliHLTMUONMansoTrackStruct* track)
{
	/// Adds a track to the internal track list for future reference in
	/// ApplyTriggerAlgorithm when we actually apply the trigger algorithm.

	assert(fTrackCount <= fMaxTracks);
	assert(fTracks != NULL);
	
	if (fTrackCount == fMaxTracks)
	{
		// Buffer full so we need to resize it.
		const AliHLTMUONMansoTrackStruct** tmp = NULL;
		try
		{
			tmp = new AliHLTMUONMansoTrackStructP[fMaxTracks+1];
		}
		catch (const std::bad_alloc&)
		{
			HLTError("Could not allocate more memory for the track array.");
			return -ENOMEM;
		}
		
		// Copy over the exisiting data and then delete the old array.
		memcpy(tmp, fTracks, sizeof(AliHLTMUONMansoTrackStructP)*fTrackCount);
		delete [] fTracks;
		fTracks = tmp;
		fMaxTracks = fMaxTracks+1;
	}
	
	fTracks[fTrackCount] = track;
	fTrackCount++;
	return 0;
}


void AliHLTMUONDecisionComponent::ApplyTriggerAlgorithm(
		AliHLTMUONSinglesDecisionBlockStruct& singlesHeader,
		AliHLTMUONTrackDecisionStruct* singlesDecision,
		AliHLTMUONPairsDecisionBlockStruct& pairsHeader,
		AliHLTMUONPairDecisionStruct* pairsDecision
	)
{
	/// This method applies the dHLT trigger decision algorithm to all the
	/// tracks found in the input data.

	// Zero the trigger counters for single tracks.
	singlesHeader.fNlowPt = 0;
	singlesHeader.fNhighPt = 0;

	// Zero the trigger counters for pairs.
	pairsHeader.fNunlikeAnyPt = 0;
	pairsHeader.fNunlikeLowPt = 0;
	pairsHeader.fNunlikeHighPt = 0;
	pairsHeader.fNlikeAnyPt = 0;
	pairsHeader.fNlikeLowPt = 0;
	pairsHeader.fNlikeHighPt = 0;
	pairsHeader.fNmassAny = 0;
	pairsHeader.fNmassLow = 0;
	pairsHeader.fNmassHigh = 0;
	
	// For the single tracks we check if a track has pT larger than either
	// the low or high pT cut. If it does then we increment the appropriate
	// counters in the header.
	for (AliHLTUInt32_t n = 0; n < fTrackCount; n++)
	{
		const AliHLTMUONMansoTrackStruct* track = fTracks[n];
		AliHLTMUONTrackDecisionStruct& decision = singlesDecision[n];
		
		bool passedHighPtCut = false;
		bool passedLowPtCut = false;
		
		AliHLTFloat32_t pt = sqrt(track->fPx * track->fPx + track->fPy * track->fPy);
		
		if (pt > fHighPtCut)
		{
			passedHighPtCut = true;
			singlesHeader.fNlowPt++;
		}
		if (pt > fLowPtCut)
		{
			passedLowPtCut = true;
			singlesHeader.fNhighPt++;
		}
		
		decision.fTrackId = track->fId;
		decision.fTriggerBits = AliHLTMUONUtils::PackTrackDecisionBits(
				passedHighPtCut, passedLowPtCut
			);
		decision.fPt = pt;
	}
	
	// Now we generate all the possible pairs of tracks and fill in the
	// trigger information. This will consist of calculating the invariant
	// mass for the pair, checking if it passes the low or high mass cut
	// and incrementing the appropriate statistics.
	AliHLTUInt32_t currentPair = 0;
	for (AliHLTUInt32_t i = 0; i < fTrackCount; i++)
	for (AliHLTUInt32_t j = i+1; j < fTrackCount; j++)
	{
		const AliHLTMUONMansoTrackStruct* tracki = fTracks[i];
		const AliHLTMUONMansoTrackStruct* trackj = fTracks[j];
		const AliHLTMUONTrackDecisionStruct& trackidecision = singlesDecision[i];
		const AliHLTMUONTrackDecisionStruct& trackjdecision = singlesDecision[j];
		AliHLTMUONPairDecisionStruct& decision = pairsDecision[currentPair];
		
		AliHLTFloat32_t muMass = 0.1056583568; // muon mass in GeV/c^2
		
		AliHLTFloat32_t mass = AliHLTMUONCalculations::ComputeMass(
				muMass, tracki->fPx, tracki->fPy, tracki->fPz,
				muMass, trackj->fPx, trackj->fPy, trackj->fPz
			);
		
		AliHLTMUONParticleSign signi, signj;
		bool hitset[4];
		AliHLTMUONUtils::UnpackMansoTrackFlags(tracki->fFlags, signi, hitset);
		AliHLTMUONUtils::UnpackMansoTrackFlags(trackj->fFlags, signj, hitset);
		
		AliHLTUInt8_t highPtCount = 0;
		if (trackidecision.fPt > fHighPtCut) highPtCount++;
		if (trackjdecision.fPt > fHighPtCut) highPtCount++;
		AliHLTUInt8_t lowPtCount = 0;
		if (trackidecision.fPt > fLowPtCut) lowPtCount++;
		if (trackjdecision.fPt > fLowPtCut) lowPtCount++;
		
		bool unlikeSign = (signi == kSignMinus and signj == kSignPlus) or
		                  (signi == kSignPlus  and signj == kSignMinus);
		
		bool passedHighMassCut = false;
		bool passedLowMassCut = false;
		if (unlikeSign)
		{
			pairsHeader.fNunlikeAnyPt++;
			if (lowPtCount == 2) pairsHeader.fNunlikeLowPt++;
			if (highPtCount == 2) pairsHeader.fNunlikeHighPt++;
			
			if (mass > fHighMassCut)
			{
				passedHighMassCut = true;
				if (highPtCount == 2) pairsHeader.fNmassHigh++;
			}
			if (mass > fLowMassCut)
			{
				passedLowMassCut = true;
				pairsHeader.fNmassAny++;
				if (lowPtCount == 2) pairsHeader.fNmassLow++;
			}
		}
		else
		{
			pairsHeader.fNlikeAnyPt++;
			if (lowPtCount == 2) pairsHeader.fNlikeLowPt++;
			if (highPtCount == 2) pairsHeader.fNlikeHighPt++;
		}
		
		decision.fTrackAId = tracki->fId;
		decision.fTrackBId = trackj->fId;
		decision.fTriggerBits = AliHLTMUONUtils::PackPairDecisionBits(
				passedHighMassCut, passedLowMassCut, unlikeSign,
				highPtCount, lowPtCount
			);
		decision.fInvMass = mass;
		
		currentPair++;
	}
	
	assert( currentPair == fTrackCount * (fTrackCount-1) / 2 );
}

