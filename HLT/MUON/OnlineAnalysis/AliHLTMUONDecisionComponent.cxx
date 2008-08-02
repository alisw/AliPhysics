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
	AliHLTMUONProcessor(),
	fMaxTracks(1),
	fTrackCount(0),
	fTracks(new AliHLTMUONMansoTrackStructP[fMaxTracks]),
	fLowPtCut(1.),  // 1 GeV/c cut
	fHighPtCut(2.),  // 2 GeV/c cut
	fLowMassCut(2.5),  // 2.7 GeV/c^2 cut
	fHighMassCut(7.),  // 8 GeV/c^2 cut
	fWarnForUnexpecedBlock(false),
	fLowPtCutSet(false),
	fHighPtCutSet(false),
	fLowMassCutSet(false),
	fHighMassCutSet(false),
	fFillSinglesDetail(false),
	fFillPairsDetail(false)
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


void AliHLTMUONDecisionComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
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
	inputMultiplier = 100;
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
	
	// Inherit the parents functionality.
	int result = AliHLTMUONProcessor::DoInit(argc, argv);
	if (result != 0) return result;

	fWarnForUnexpecedBlock = false;
	fLowPtCutSet = false;
	fHighPtCutSet = false;
	fLowMassCutSet = false;
	fHighMassCutSet = false;
	fFillSinglesDetail = true;
	fFillPairsDetail = true;
	
	for (int i = 0; i < argc; i++)
	{
		if (ArgumentAlreadyHandled(i, argv[i])) continue;

		if (strcmp( argv[i], "-lowptcut" ) == 0)
		{
			if (fLowPtCutSet)
			{
				HLTWarning("Low pT cut parameter was already specified."
					" Will replace previous value given by -lowptcut."
				);
			}
			
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
			fLowPtCutSet = true;
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-highptcut" ) == 0)
		{
			if (fHighPtCutSet)
			{
				HLTWarning("High pT cut parameter was already specified."
					" Will replace previous value given by -highptcut."
				);
			}
			
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
			fHighPtCutSet = true;
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-lowmasscut" ) == 0)
		{
			if (fLowMassCutSet)
			{
				HLTWarning("Low invariant mass cut parameter was already specified."
					" Will replace previous value given by -lowmasscut."
				);
			}
			
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
			fLowMassCutSet = true;
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-highmasscut" ) == 0)
		{
			if (fHighMassCutSet)
			{
				HLTWarning("High invariant mass cut parameter was already specified."
					" Will replace previous value given by -highmasscut."
				);
			}
			
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
			fHighMassCutSet = true;
			
			i++;
			continue;
		}
		
		if (strcmp(argv[i], "-warn_on_unexpected_block") == 0)
		{
			fWarnForUnexpecedBlock = true;
			continue;
		}
		
		if (strcmp( argv[i], "-no_singles_detail" ) == 0)
		{
			fFillSinglesDetail = false;
			continue;
		}
		
		if (strcmp( argv[i], "-no_pairs_detail" ) == 0)
		{
			fFillPairsDetail = false;
			continue;
		}

		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}
	
	if (not DelaySetup())
	{
		// Read cut parameters from CDB if they were not specified on the command line.
		if (not fLowPtCutSet or not fHighPtCutSet or not fLowMassCutSet or not fHighMassCutSet)
		{
			HLTInfo("Loading cut parameters from CDB.");
			int result = ReadConfigFromCDB(
					not fLowPtCutSet, not fHighPtCutSet,
					not fLowMassCutSet, not fHighMassCutSet
				);
			if (result != 0) return result;
		}
		else
		{
			// Print the debug messages here since ReadConfigFromCDB does not get called,
			// in-which the debug messages would have been printed.
			HLTDebug("Using the following cut parameters:");
			HLTDebug("              Low pT cut = %f GeV/c", fLowPtCut);
			HLTDebug("             High pT cut = %f GeV/c", fHighPtCut);
			HLTDebug("  Low invariant mass cut = %f GeV/c^2", fLowMassCut);
			HLTDebug(" High invariant mass cut = %f GeV/c^2", fHighMassCut);
		}
	}
	
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
	
	/// Inherited from AliHLTComponent. This method will reload CDB configuration
	/// entries for this component from the CDB.
	/// \param cdbEntry If this is NULL or equals "HLT/ConfigMUON/DecisionComponent"
	///     then new configuration parameters are loaded, otherwise nothing is done.
	/// \param componentId  The name of the component in the current chain.
	
	bool givenConfigPath = strcmp(cdbEntry, AliHLTMUONConstants::DecisionComponentCDBPath()) == 0;
	
	if (cdbEntry == NULL or givenConfigPath)
	{
		HLTInfo("Reading new configuration entries from CDB for component '%s'.", componentId);
		int result = ReadConfigFromCDB();
		if (result != 0) return result;
	}
	
	return 0;
}


int AliHLTMUONDecisionComponent::ReadPreprocessorValues(const char* modules)
{
	/// Inherited from AliHLTComponent. 
	/// Updates the configuration of this component if HLT or ALL has been
	/// specified in the 'modules' list.

	TString mods = modules;
	if (mods.Contains("ALL"))
	{
		return Reconfigure(NULL, GetComponentID());
	}
	if (mods.Contains("HLT"))
	{
		return Reconfigure(AliHLTMUONConstants::DecisionComponentCDBPath(), GetComponentID());
	}
	return 0;
}


int AliHLTMUONDecisionComponent::DoEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr,
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks
	)
{
	///
	/// Inherited from AliHLTProcessor. Processes the new event data.
	///
	
	// Initialise the cut parameters from CDB if we were requested to
	// initialise only when the first event was received.
	if (DelaySetup())
	{
		// Load the cut paramters from CDB if they have not been given
		// on the command line.
		if (not fLowPtCutSet or not fHighPtCutSet or not fLowMassCutSet or not fHighMassCutSet)
		{
			HLTInfo("Loading cut parameters from CDB.");
			int result = ReadConfigFromCDB(
					not fLowPtCutSet, not fHighPtCutSet,
					not fLowMassCutSet, not fHighMassCutSet
				);
			if (result != 0) return result;
		}
		
		DoneDelayedSetup();
	}
	
	AliHLTUInt32_t specification = 0;  // Contains the output data block spec bits.
	
	// Loop over all input blocks in the event with track data and add pointers
	// to the tracks into the tracks array. These will be used later by the
	// trigger algorithm to get to the individual tracks.
	fTrackCount = 0; // reset number of tracks in array.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
		HLTDebug("Handling block: %u, with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
			n, DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fPtr, blocks[n].fSize
		);
		
		if (blocks[n].fDataType == AliHLTMUONConstants::MansoTracksBlockDataType())
		{
			// Build up the specification which indicates what DDLs
			// contributed to the output data.
			specification |= blocks[n].fSpecification;
			
			AliHLTMUONMansoTracksBlockReader inblock(blocks[n].fPtr, blocks[n].fSize);
			if (not BlockStructureOk(inblock))
			{
				if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
				continue;
			}
			
			for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
			{
				int result = AddTrack(&inblock[i]);
				if (result != 0)
				{
					if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
					size = 0; // Important to tell framework that nothing was generated.
					return result;
				}
			}
		}
		else
		{
			// Log a message indicating that we got a data block that we
			// do not know how to handle.
			if (fWarnForUnexpecedBlock)
				HLTWarning("Received a data block of a type we cannot handle: '%s', spec: 0x%X",
					DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fSpecification
				);
			else
				HLTDebug("Received a data block of a type we cannot handle: '%s', spec: 0x%X",
					DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fSpecification
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
		if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
		size = 0; // Important to tell framework that nothing was generated.
		return -ENOBUFS;
	}

	AliHLTUInt32_t numOfTracks = fTrackCount;
	if (not fFillSinglesDetail) numOfTracks = 0;
	if (not singlesBlock.SetNumberOfEntries(numOfTracks))
	{
		AliHLTUInt32_t bytesneeded = sizeof(AliHLTMUONSinglesDecisionBlockWriter::HeaderType)
			+ numOfTracks * sizeof(AliHLTMUONSinglesDecisionBlockWriter::ElementType);
		HLTError("The buffer is only %d bytes in size. We need a minimum of"
			" %d bytes for the singles output data block.",
			size, bytesneeded
		);
		if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
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
		if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
		size = 0; // Important to tell framework that nothing was generated.
		return -ENOBUFS;
	}

	AliHLTUInt32_t numOfPairs = fTrackCount * (fTrackCount-1) / 2;
	if (not fFillPairsDetail) numOfPairs = 0;
	if (not pairsBlock.SetNumberOfEntries(numOfPairs))
	{
		AliHLTUInt32_t bytesneeded = sizeof(AliHLTMUONPairsDecisionBlockWriter::HeaderType)
			+ numOfPairs * sizeof(AliHLTMUONPairsDecisionBlockWriter::ElementType)
			+ singlesBlock.BytesUsed();
		HLTError("The buffer is only %d bytes in size. We need a minimum of"
			" %d bytes for the pairs output data block.",
			size, bytesneeded
		);
		if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
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
		bool setLowPtCut, bool setHighPtCut,
		bool setLowMassCut, bool setHighMassCut
	)
{
	/// Reads the cut parameters from the CDB.
	/// \param setLowPtCut  Indicates if the low pT cut should be set (default true).
	/// \param setHighPtCut  Indicates if the high pT cut should be set (default true).
	/// \param setLowMassCut  Indicates if the low invariant mass cut should be set (default true).
	/// \param setHighMassCut  Indicates if the high invariant mass cut should be set (default true).
	/// \return 0 is returned on success and a non-zero value to indicate failure.
	
	assert(AliCDBManager::Instance() != NULL);
	
	const char* pathToEntry = AliHLTMUONConstants::DecisionComponentCDBPath();
	
	TMap* map = NULL;
	int result = FetchTMapFromCDB(pathToEntry, map);
	if (result != 0) return result;
	
	if (setLowPtCut)
	{
		Double_t value = 0;
		result = GetFloatFromTMap(map, "lowptcut", value, pathToEntry, "low pT cut");
		if (result != 0) return result;
		fLowPtCut = (AliHLTFloat32_t) value;
	}
	
	if (setHighPtCut)
	{
		Double_t value = 0;
		result = GetFloatFromTMap(map, "highptcut", value, pathToEntry, "high pT cut");
		if (result != 0) return result;
		fHighPtCut = (AliHLTFloat32_t) value;
	}
	
	if (setLowMassCut)
	{
		Double_t value = 0;
		result = GetFloatFromTMap(map, "lowmasscut", value, pathToEntry, "low invariant mass cut");
		if (result != 0) return result;
		fLowMassCut = (AliHLTFloat32_t) value;
	}
	
	if (setHighMassCut)
	{
		Double_t value = 0;
		result = GetFloatFromTMap(map, "highmasscut", value, pathToEntry, "high invariant mass cut");
		if (result != 0) return result;
		fHighMassCut = (AliHLTFloat32_t) value;
	}
	
	HLTDebug("Using the following cut parameters:");
	HLTDebug("              Low pT cut = %f GeV/c", fLowPtCut);
	HLTDebug("             High pT cut = %f GeV/c", fHighPtCut);
	HLTDebug("  Low invariant mass cut = %f GeV/c^2", fLowMassCut);
	HLTDebug(" High invariant mass cut = %f GeV/c^2", fHighMassCut);
	
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


int AliHLTMUONDecisionComponent::ApplyTriggerAlgorithm(
		AliHLTMUONSinglesDecisionBlockStruct& singlesHeader,
		AliHLTMUONTrackDecisionStruct* singlesDecision,
		AliHLTMUONPairsDecisionBlockStruct& pairsHeader,
		AliHLTMUONPairDecisionStruct* pairsDecision
	)
{
	/// This method applies the dHLT trigger decision algorithm to all the
	/// tracks found in the input data.
	/// @return zero on success and -ENOMEM if out of memory.

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

	// Allocate a temporary memory buffer for the calculated pT values if
	// we are not storing them to shared memory as part of the new block.
	AliHLTFloat32_t* ptValues = NULL;
	if (not fFillSinglesDetail)
	{
		try
		{
			ptValues = new AliHLTFloat32_t[fTrackCount];
		}
		catch(const std::bad_alloc&)
		{
			HLTError("Could not allocate memory buffer for pT values.");
			return -ENOMEM;
		}
	}
	
	// For the single tracks we check if a track has pT larger than either
	// the low or high pT cut. If it does then we increment the appropriate
	// counters in the header.
	for (AliHLTUInt32_t n = 0; n < fTrackCount; n++)
	{
		const AliHLTMUONMansoTrackStruct* track = fTracks[n];
		
		bool passedHighPtCut = false;
		bool passedLowPtCut = false;
		
		AliHLTFloat32_t pt = sqrt(track->fPx * track->fPx + track->fPy * track->fPy);
		
		if (pt > fHighPtCut)
		{
			passedHighPtCut = true;
			singlesHeader.fNhighPt++;
		}
		if (pt > fLowPtCut)
		{
			passedLowPtCut = true;
			singlesHeader.fNlowPt++;
		}
		
		if (fFillSinglesDetail)
		{
			AliHLTMUONTrackDecisionStruct& decision = singlesDecision[n];
			decision.fTrackId = track->fId;
			decision.fTriggerBits = AliHLTMUONUtils::PackTrackDecisionBits(
					passedHighPtCut, passedLowPtCut
				);
			decision.fPt = pt;
		}
		else
		{
			ptValues[n] = pt;
		}
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
		AliHLTUInt8_t lowPtCount = 0;
		if (fFillSinglesDetail)
		{
			const AliHLTMUONTrackDecisionStruct& trackidecision = singlesDecision[i];
			const AliHLTMUONTrackDecisionStruct& trackjdecision = singlesDecision[j];
			if (trackidecision.fPt > fHighPtCut) highPtCount++;
			if (trackjdecision.fPt > fHighPtCut) highPtCount++;
			if (trackidecision.fPt > fLowPtCut) lowPtCount++;
			if (trackjdecision.fPt > fLowPtCut) lowPtCount++;
		}
		else
		{
			if (ptValues[i] > fHighPtCut) highPtCount++;
			if (ptValues[j] > fHighPtCut) highPtCount++;
			if (ptValues[i] > fLowPtCut) lowPtCount++;
			if (ptValues[j] > fLowPtCut) lowPtCount++;
		}
		
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
		
		if (fFillPairsDetail)
		{
			AliHLTMUONPairDecisionStruct& decision = pairsDecision[currentPair];

			decision.fTrackAId = tracki->fId;
			decision.fTrackBId = trackj->fId;
			decision.fTriggerBits = AliHLTMUONUtils::PackPairDecisionBits(
					passedHighMassCut, passedLowMassCut, unlikeSign,
					highPtCount, lowPtCount
				);
			decision.fInvMass = mass;
		
			currentPair++;
		}
	}
	
	assert( fFillPairsDetail == false or (fFillPairsDetail == true and currentPair == fTrackCount * (fTrackCount-1) / 2) );

	if (not fFillSinglesDetail)
	{
		assert(ptValues != NULL);
		delete [] ptValues;
	}

	return 0;
}

