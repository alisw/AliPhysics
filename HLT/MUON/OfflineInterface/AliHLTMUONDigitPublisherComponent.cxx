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

/* $Id: AliHLTMUONDigitPublisherComponent.cxx 26179 2008-05-29 22:27:27Z aszostak $ */

///
/// @file   AliHLTMUONDigitPublisherComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 May 2008
/// @brief  Implementation of the dHLT digit publisher component.
///
/// This component is used to publish simulated or reconstructed digits from
/// the digits trees as DDL raw data. The data is converted into DDL format
/// on the fly.
///

#include "AliHLTMUONDigitPublisherComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTLogging.h"
#include "AliHLTSystem.h"
#include "AliHLTDefinitions.h"
#include "AliRawDataHeader.h"
#include "AliMUONTrackerDDLDecoderEventHandler.h"
#include "AliMUONConstants.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONDataInterface.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMpExMap.h"
#include "AliMpCDB.h"
#include "AliMpDDL.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpTriggerCrate.h"
#include "AliMpConstants.h"
#include "AliBitPacking.h"
#include "AliMUONBlockHeader.h"
#include "AliMUONBusStruct.h"
#include "AliMUONConstants.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONVDigit.h"
#include "AliMUONDspHeader.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONRegHeader.h"
#include "AliRunLoader.h"
#include "AliCentralTrigger.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <cassert>


ClassImp(AliHLTMUONDigitPublisherComponent)


AliHLTMUONDigitPublisherComponent::AliHLTMUONDigitPublisherComponent() :
	AliHLTOfflineDataSource(),
	fDDL(-1),
	fCurrentEventIndex(0),
	fMakeScalars(false),
	fMCDataInterface(NULL),
	fDataInterface(NULL),
	fChamberExclusionList(0),
	fDetElemExclusionList(0)
{
	/// Default constructor.
}


AliHLTMUONDigitPublisherComponent::~AliHLTMUONDigitPublisherComponent()
{
	/// Default destructor.
	
	if (fMCDataInterface != NULL) delete fMCDataInterface;
	if (fDataInterface != NULL) delete fDataInterface;
}

const char* AliHLTMUONDigitPublisherComponent::GetComponentID()
{
	/// Inherited from AliHLTComponent. Returns the component ID.
	
	return AliHLTMUONConstants::DigitPublisherId();
}


AliHLTComponentDataType AliHLTMUONDigitPublisherComponent::GetOutputDataType()
{
	/// Inherited from AliHLTComponent. Returns the raw DDL data type.
	
	return AliHLTMUONConstants::DDLRawDataType();
}


void AliHLTMUONDigitPublisherComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	/// Inherited from AliHLTComponent.
	/// Returns an estimate of the expected output data size.
	
	// estimated as max number of channels * raw data word size + max headers size.
	constBase = sizeof(AliRawDataHeader) + 65536*sizeof(UInt_t)
		+ sizeof(AliMUONBlockHeaderStruct)*2 + sizeof(AliMUONDSPHeaderStruct)*10
		+ sizeof(AliMUONBusPatchHeaderStruct) * 50;
	inputMultiplier = 0;
}


AliHLTComponent* AliHLTMUONDigitPublisherComponent::Spawn()
{
	/// Inherited from AliHLTComponent. Creates a new object instance.
	
	return new AliHLTMUONDigitPublisherComponent;
}


int AliHLTMUONDigitPublisherComponent::ParseChamberString(const char* str)
{
	/// Parses a string with the following format:
	///   <number>|<number>-<number>[,<number>|<number>-<number>,...]
	/// For example: 1  1,2,3  1-2   1,2-4,5  etc...
	/// Chamber numbers must be in the range [1..10] for tracking chambers.
	/// All valid tracking chamber numbers will added to fChamberExclusionList.
	/// @param str  The string to parse.
	/// @return  Zero on success and EINVAL if there is a parse error.
	
	char* end = const_cast<char*>(str);
	long lastChamber = -1;
	do
	{
		// Parse the next number.
		char* current = end;
		long chamber = strtol(current, &end, 0);
		
		// Check for parse errors of the number.
		if (current == end)
		{
			HLTError("Expected a number in the range [1..%d] but got '%s'.",
				AliMUONConstants::NTrackingCh(), current
			);
			return -EINVAL;
		}
		if (chamber < 1 or AliMUONConstants::NTrackingCh() < chamber)
		{
			HLTError("Received the chamber number %d, which is outside the valid range of [1..%d].",
				chamber, AliMUONConstants::NTrackingCh()
			);
			return -EINVAL;
		}
		
		// Skip any whitespace after the number
		while (*end != '\0' and (*end == ' ' or *end == '\t' or *end == '\r' or *end == '\n')) end++;
		
		// Check if we are dealing with a list or range, or if we are at
		// the end of the string.
		if (*end == '-')
		{
			lastChamber = chamber;
			end++;
			continue;
		}
		else if (*end == ',')
		{
			assert( 1 <= chamber and chamber <= AliMUONConstants::NTrackingCh() );
			Int_t size = fChamberExclusionList.GetSize();
			fChamberExclusionList.Set(size+1);
			fChamberExclusionList[size] = chamber-1;
			end++;
		}
		else if (*end == '\0')
		{
			assert( 1 <= chamber and chamber <= AliMUONConstants::NTrackingCh() );
			Int_t size = fChamberExclusionList.GetSize();
			fChamberExclusionList.Set(size+1);
			fChamberExclusionList[size] = chamber-1;
		}
		else
		{
			HLTError("Could not understand parameter list '%s'. Expected '-', ','"
				  " or end of line, but received '%c' at character %d.",
				str, *end, (int)(end - str) +1
			);
			return -EINVAL;
		}
		
		// Set the range of chambers to publish for.
		if (lastChamber > 0)
		{
			Int_t min, max;
			if (lastChamber < chamber)
			{
				min = lastChamber;
				max = chamber;
			}
			else
			{
				min = chamber;
				max = lastChamber;
			}
			assert( min >= 1 );
			assert( max <= AliMUONConstants::NTrackingCh() );
			for (Int_t i = min; i <= max; i++)
			{
				Int_t size = fChamberExclusionList.GetSize();
				fChamberExclusionList.Set(size+1);
				fChamberExclusionList[size] = i-1;
			}
		}
		lastChamber = -1;
	}
	while (*end != '\0');
	return 0;
}


int AliHLTMUONDigitPublisherComponent::ParseDetElemString(const char* str)
{
	/// Parses a string with the following format:
	///   <number>|<number>-<number>[,<number>|<number>-<number>,...]
	/// For example: 100  100,201,208  100-104   105,202-204,503  etc...
	/// Detector element numbers must be in the range [100..1099] for tracking stations.
	/// All valid detector element numbers will added to fDetElemExclusionList.
	/// @param str  The string to parse.
	/// @return  Zero on success and EINVAL if there is a parse error.
	
	char* end = const_cast<char*>(str);
	long lastDetElem = -1;
	do
	{
		// Parse the next number.
		char* current = end;
		long detElem = strtol(current, &end, 0);
		
		// Check for parse errors of the number.
		if (current == end)
		{
			HLTError("Expected a number in the range [100..1099] but got '%s'.",
				current
			);
			return -EINVAL;
		}
		if (detElem < 100 or 1099 < detElem)
		{
			HLTError("Received the detector element ID number of %d,"
				" which is outside the valid range of [100..1099].",
				detElem
			);
			return -EINVAL;
		}
		
		// Skip any whitespace after the number
		while (*end != '\0' and (*end == ' ' or *end == '\t' or *end == '\r' or *end == '\n')) end++;
		
		// Check if we are dealing with a list or range, or if we are at
		// the end of the string.
		if (*end == '-')
		{
			lastDetElem = detElem;
			end++;
			continue;
		}
		else if (*end == ',')
		{
			assert( 100 <= detElem and detElem <= 1099 );
			Int_t size = fDetElemExclusionList.GetSize();
			fDetElemExclusionList.Set(size+1);
			fDetElemExclusionList[size] = detElem-1;
			end++;
		}
		else if (*end == '\0')
		{
			assert( 100 <= detElem and detElem <= 1099 );
			Int_t size = fDetElemExclusionList.GetSize();
			fDetElemExclusionList.Set(size+1);
			fDetElemExclusionList[size] = detElem-1;
		}
		else
		{
			HLTError("Could not understand parameter list '%s'. Expected '-', ','"
				  " or end of line, but received '%c' at character %d.",
				str, *end, (int)(end - str) +1
			);
			return -EINVAL;
		}
		
		// Set the range of detector elements to publish for.
		if (lastDetElem > 0)
		{
			Int_t min, max;
			if (lastDetElem < detElem)
			{
				min = lastDetElem;
				max = detElem;
			}
			else
			{
				min = detElem;
				max = lastDetElem;
			}
			assert( min >= 100 );
			assert( max <= 1099 );
			for (Int_t i = min; i <= max; i++)
			{
				Int_t size = fDetElemExclusionList.GetSize();
				fDetElemExclusionList.Set(size+1);
				fDetElemExclusionList[size] = i-1;
			}
		}
		lastDetElem = -1;
	}
	while (*end != '\0');
	return 0;
}


int AliHLTMUONDigitPublisherComponent::DoInit(int argc, const char** argv)
{
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	
	HLTInfo("Initialising dHLT digit publisher component.");

	if (fMCDataInterface != NULL)
	{
		delete fMCDataInterface;
		fMCDataInterface = NULL;
	}
	if (fDataInterface != NULL)
	{
		delete fDataInterface;
		fDataInterface = NULL;
	}
	
	// Initialise with default values.
	fDDL = -1;
	fCurrentEventIndex = 0;
	fMakeScalars = false;
	fChamberExclusionList.Set(0);
	fDetElemExclusionList.Set(0);
	bool simdata = false;
	bool recdata = false;
	bool firstEventSet = false;
	bool eventNumLitSet = false;

	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-makescalars") == 0)
		{
			fMakeScalars = true;
			continue;
		}
		if (strcmp(argv[i], "-simdata") == 0)
		{
			simdata = true;
			continue;
		}
		if (strcmp(argv[i], "-recdata") == 0)
		{
			recdata = true;
			continue;
		}
		if (strcmp(argv[i], "-ddl") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("DDL number not specified. It must be in the range [1..22]" );
				return -EINVAL;
			}
		
			char* cpErr = NULL;
			unsigned long num = strtoul(argv[i+1], &cpErr, 0);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a DDL number.", argv[i+1]);
				return -EINVAL;
			}
			if (num < 1 or 22 < num)
			{
				HLTError("The DDL number must be in the range [1..22].");
				return -EINVAL;
			}
			fDDL = num - 1; // Convert to DDL number in the range 0..21
			
			i++;
			continue;
		}
		if (strcmp(argv[i], "-ddlid") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("DDL equipment ID number not specified."
					" It must be in the range [2560..2579] or [2816..2817]."
				);
				return -EINVAL;
			}
		
			char* cpErr = NULL;
			unsigned long num = strtoul(argv[i+1], &cpErr, 0);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a DDL equipment ID number.", argv[i+1]);
				return -EINVAL;
			}
			fDDL = AliHLTMUONUtils::EquipIdToDDLNumber(num); // Convert to DDL number in the range 0..21
			if (fDDL < 0 or 21 < fDDL)
			{
				HLTError("The DDL equipment ID number must be in the range"
					" [2560..2579] or [2816..2817]."
				);
				return -EINVAL;
			}
			
			i++;
			continue;
		}
		if (strcmp(argv[i], "-firstevent") == 0)
		{
			if (eventNumLitSet)
			{
				HLTWarning("The -firstevent flag is overridden by a"
					" previous use of -event_number_literal."
				);
			}
			if (++i >= argc)
			{
				HLTError("Expected a positive number after -firstevent.");
				return -EINVAL;
			}
			char* end = NULL;
			long num = strtol(argv[i], &end, 0);
			if ((end != NULL and *end != '\0') or num < 0) // Check if the conversion is OK.
			{
				HLTError("Expected a positive number after -firstevent"
					" but got: %s", argv[i]
				);
				return -EINVAL;
			}
			fCurrentEventIndex = Int_t(num);
			firstEventSet = true;
			continue;
		}
		if (strcmp(argv[i], "-event_number_literal") == 0)
		{
			if (firstEventSet)
			{
				HLTWarning("The -event_number_literal option will"
					" override -firstevent."
				);
			}
			fCurrentEventIndex = -1;
			eventNumLitSet = true;
			continue;
		}
		if (strcmp(argv[i], "-exclude_chamber") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("Expected a chamber number, a range eg. '1-10', or a list eg."
					" '1,2,3' after '-exclude_chamber'."
				);
				return -EINVAL;
			}
			
			int result = ParseChamberString(argv[i+1]);
			if (result != 0) return result;
			i++;
			continue;
		}
		if (strcmp(argv[i], "-exclude_detelem") == 0)
		{
			if (argc <= i+1)
			{
				HLTError("Expected a detector element ID number, a range eg. '100-108',"
					" or a list eg. '100,102,301' after '-exclude_detelem'."
				);
				return -EINVAL;
			}
			
			int result = ParseDetElemString(argv[i+1]);
			if (result != 0) return result;
			i++;
			continue;
		}
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}
	
	if (fDDL == -1)
	{
		HLTError("DDL number must be set with the -ddl option, but it was not.");
		return -EINVAL;
	}
	
	// Must load the mapping data if it is not already loaded.
	if (AliMpDDLStore::Instance(false) == NULL)
	{
		AliMpCDB::LoadDDLStore();
		if (AliMpDDLStore::Instance(false) == NULL)
		{
			HLTError("Could not load the DDL mapping store from CDB.");
			return -EFAULT;
		}
	}
	
	// Now we can initialise the data interface objects and loaders.
	if (simdata)
	{
		HLTDebug("Loading simulated digits with AliMUONMCDataInterface.");

		try
		{
			fMCDataInterface = new AliMUONMCDataInterface("galice.root");
		}
		catch (const std::bad_alloc&)
		{
			HLTError("Not enough memory to allocate AliMUONMCDataInterface.");
			return -ENOMEM;
		}
	}
	else if (recdata)
	{
		HLTDebug("Loading reconstructed digits with AliMUONDataInterface.");
		
		try
		{
			fDataInterface = new AliMUONDataInterface("galice.root");
		}
		catch (const std::bad_alloc&)
		{
			HLTError("Not enough memory to allocate AliMUONDataInterface.");
			return -ENOMEM;
		}
	}
	
	// Check that the fCurrentEventIndex number falls within the correct range.
	UInt_t maxevent = 0;
	if (fMCDataInterface != NULL)
		maxevent = UInt_t(fMCDataInterface->NumberOfEvents());
	else if (fDataInterface != NULL)
		maxevent = UInt_t(fDataInterface->NumberOfEvents());
	if (fCurrentEventIndex != -1 and UInt_t(fCurrentEventIndex) >= maxevent and maxevent != 0)
	{
		fCurrentEventIndex = 0;
		HLTWarning("The selected first event number (%d) was larger than"
			" the available number of events (%d). Resetting the event"
			" counter to zero.", fCurrentEventIndex, maxevent
		);
	}

	return 0;
}


int AliHLTMUONDigitPublisherComponent::DoDeinit()
{
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	
	HLTInfo("Deinitialising dHLT digit publisher component.");
	
	if (fMCDataInterface != NULL)
	{
		delete fMCDataInterface;
		fMCDataInterface = NULL;
	}
	if (fDataInterface != NULL)
	{
		delete fDataInterface;
		fDataInterface = NULL;
	}
	return 0;
}


int AliHLTMUONDigitPublisherComponent::GetEvent(
		const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& /*trigData*/,
		AliHLTUInt8_t* outputPtr,
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks
	)
{
	/// Inherited from AliHLTOfflineDataSource.
	
	assert( fMCDataInterface != NULL or fDataInterface != NULL );
	
	if (not IsDataEvent()) return 0;  // ignore non data events.

	// Check the size of the event descriptor structure.
	if (evtData.fStructSize < sizeof(AliHLTComponentEventData))
	{
		HLTError(kHLTLogError,
			"The event descriptor (AliHLTComponentEventData) size is"
			  " smaller than expected. It claims to be %d bytes, but"
			  " we expect it to be %d bytes.",
			evtData.fStructSize,
			sizeof(AliHLTComponentEventData)
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -EINVAL;
	}
	
	// Use the fEventID as the event number to load if fCurrentEventIndex == -1,
	// check it and load that event with the runloader.
	// If fCurrentEventIndex is a positive number then use it instead and
	// increment it.
	UInt_t eventnumber = UInt_t(evtData.fEventID);
	UInt_t maxevent = 0;
	if (fMCDataInterface != NULL)
		maxevent = UInt_t(fMCDataInterface->NumberOfEvents());
	else if (fDataInterface != NULL)
		maxevent = UInt_t(fDataInterface->NumberOfEvents());
	if (fCurrentEventIndex != -1)
	{
		eventnumber = UInt_t(fCurrentEventIndex);
		fCurrentEventIndex++;
		if (UInt_t(fCurrentEventIndex) >= maxevent)
			fCurrentEventIndex = 0;
	}
	if ( eventnumber >= maxevent )
	{
		HLTError("The event number (%d) is larger than the available number"
			  " of events on file (%d).",
			eventnumber, maxevent
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -EINVAL;
	}
	
	const AliMUONVDigitStore* digitStore = NULL;
	const AliMUONVTriggerStore* triggerStore = NULL;
	
	if (fMCDataInterface != NULL)
	{
		HLTDebug("Filling data block with simulated digits for event %d.", eventnumber);
		
		if (fDDL < 20)
		{
			digitStore = fMCDataInterface->DigitStore(eventnumber);
		}
		else
		{
			triggerStore = fMCDataInterface->TriggerStore(eventnumber);
		}
	}
	else if (fDataInterface != NULL)
	{
		HLTDebug("Filling data block with reconstructed digits for event %d.", eventnumber);
		
		if (fDDL < 20)
		{
			digitStore = fDataInterface->DigitStore(eventnumber);
		}
		else
		{
			triggerStore = fDataInterface->TriggerStore(eventnumber);
		}
	}
	else
	{
		HLTError("Neither AliMUONDataInterface nor AliMUONMCDataInterface were created.");
		size = 0; // Important to tell framework that nothing was generated.
		return -EFAULT;
	}
	
	// Make sure we have the correct CTP trigger loaded.
#ifndef HAVE_NOT_ALIRUNLOADER30859
	AliRunLoader* runloader = AliRunLoader::Instance();
#else
	// the old way before rev 30859
	AliRunLoader *runloader = AliRunLoader::GetRunLoader();
#endif
	if (runloader != NULL)
	{
		if (runloader->GetTrigger() == NULL)
			runloader->LoadTrigger();
		runloader->GetEvent(eventnumber);
	}
	
	if (fDDL < 20 and digitStore != NULL)
	{
		int result = WriteTrackerDDL(digitStore, fDDL, outputPtr, size);
		if (result != 0)
		{
			size = 0; // Important to tell framework that nothing was generated.
			return result;
		}
	}
	else if (triggerStore != NULL)
	{
		int result = WriteTriggerDDL(triggerStore, fDDL-20, outputPtr, size, fMakeScalars);
		if (result != 0)
		{
			size = 0; // Important to tell framework that nothing was generated.
			return result;
		}
	}
	else
	{
		size = 0; // Important to tell framework that nothing was generated.
		return 0;
	}
	
	AliHLTComponentBlockData bd;
	FillBlockData(bd);
	bd.fPtr = outputPtr;
	bd.fOffset = 0;
	bd.fSize = size;
	bd.fDataType = AliHLTMUONConstants::DDLRawDataType();
	bd.fSpecification = AliHLTMUONUtils::DDLNumberToSpec(fDDL);
	outputBlocks.push_back(bd);
	
	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////
// Methods copied from AliMUONRawWriter.
//TODO: This is not ideal. We should have AliMUONRawWriter re-factored so that
// we can have raw data generated into a memory resident buffer, rather than
// always written to a file on disk, as it is now. But this will take some time
// since people need to be convinced of this fact.

//____________________________________________________________________
void  AliHLTMUONDigitPublisherComponent::LocalWordPacking(UInt_t& word, UInt_t locId, UInt_t locDec, 
					 UInt_t trigY, UInt_t posY, UInt_t posX, 
					 UInt_t sdevX, UInt_t devX)
{
/// pack local trigger word

    AliBitPacking::PackWord(locId,word,19,22); //card id number in crate
    AliBitPacking::PackWord(locDec,word,15,18);
    AliBitPacking::PackWord(trigY,word,14,14);
    AliBitPacking::PackWord(posY,word,10,13);
    AliBitPacking::PackWord(sdevX,word,9,9);
    AliBitPacking::PackWord(devX,word,5,8);
    AliBitPacking::PackWord(posX,word,0,4);
}

//______________________________________________________________________________
void 
AliHLTMUONDigitPublisherComponent::Digits2BusPatchMap(
		const AliMUONVDigitStore& digitStore,
		AliMpExMap& busPatchMap, Int_t iDDL
	)
{
  /// Create bus patch structures corresponding to digits in the store
  
  AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
  assert(ddlStore != NULL);
  
  AliMpDDL* ddl = ddlStore->GetDDL(iDDL);
  busPatchMap.SetSize(ddl->GetNofBusPatches());
  
  if (ddl->GetNofDEs() <= 0) return;
  Int_t minDetElem = ddl->GetDEId(0);
  Int_t maxDetElem = ddl->GetDEId(0);
  for (Int_t i = 1; i < ddl->GetNofDEs(); i++)
  {
    if (ddl->GetDEId(i) < minDetElem) minDetElem = ddl->GetDEId(i);
    if (ddl->GetDEId(i) > maxDetElem) maxDetElem = ddl->GetDEId(i);
  }
  
  static const Int_t kMAXADC = (1<<12)-1; // We code the charge on a 12 bits ADC.
    
  // DDL event one per half chamber
  
  // raw data
  Char_t parity = 0x4;
  UShort_t manuId = 0;
  UChar_t channelId = 0;
  UShort_t charge = 0;
  Int_t busPatchId = 0;
  Int_t currentBusPatchId = -1;
  UInt_t word;
  
  AliMUONBusStruct* busStruct(0x0);
  
  TIter next(digitStore.CreateIterator(minDetElem, maxDetElem));
  AliMUONVDigit* digit;
  
  while ( ( digit = static_cast<AliMUONVDigit*>(next()) ) )
  {
    // Check if we should exclude digits from a particular chamber or detector element.
    bool excludeDigit = false;
    for (Int_t i = 0; i < fDetElemExclusionList.GetSize(); i++)
    {
      if (digit->DetElemId() == fDetElemExclusionList[i])
      {
        excludeDigit = true;
        break;
      }
    }
    for (Int_t i = 0; i < fChamberExclusionList.GetSize(); i++)
    {
      if (AliMpDEManager::GetChamberId(digit->DetElemId()) == fChamberExclusionList[i])
      {
        excludeDigit = true;
        break;
      }
    }
    if (excludeDigit) continue;
  
    charge = digit->ADC();
    if ( charge > kMAXADC )
    {
      // This is most probably an error in the digitizer (which should insure
      // the adc is below kMAXADC), so make it a (non-fatal) error indeed.
      HLTError("ADC value %d above 0x%x for DE %d . Setting to 0x%x. Digit is:",
                    charge,kMAXADC,digit->DetElemId(),kMAXADC);
      charge = kMAXADC;
    }
    
    // inverse mapping
    busPatchId = ddlStore->GetBusPatchId(digit->DetElemId(), digit->ManuId());

    if (busPatchId<0) continue;
    
    if ( digit->ManuId() > 0x7FF ||
         digit->ManuChannel() > 0x3F )
    {
      HLTFatal("<%s>: ID %12u DE %4d Cath %d (Ix,Iy)=(%3d,%3d) (Manu,Channel)=(%4d,%2d)"
               ", Charge=%7.2f\nManuId,ManuChannel are invalid for this digit.",
               digit->ClassName(), digit->GetUniqueID(),
               digit->DetElemId(), digit->Cathode(), digit->PadX(), digit->PadY(),
               digit->ManuId(), digit->ManuChannel(), digit->Charge()
      );
    }
    
    manuId = ( digit->ManuId() & 0x7FF ); // 11 bits
    channelId = ( digit->ManuChannel() & 0x3F ); // 6 bits
    
    //packing word
    word = 0;
    AliBitPacking::PackWord((UInt_t)manuId,word,18,28);
    AliBitPacking::PackWord((UInt_t)channelId,word,12,17);
    AliBitPacking::PackWord((UInt_t)charge,word,0,11);
    
    // parity word
    parity = word & 0x1;
    for (Int_t i = 1; i <= 30; ++i) 
    {
      parity ^=  ((word >> i) & 0x1);
    }
    AliBitPacking::PackWord((UInt_t)parity,word,31,31);

    if ( currentBusPatchId != busPatchId ) 
    {
      busStruct = 
        static_cast<AliMUONBusStruct*>(busPatchMap.GetValue(busPatchId));
      currentBusPatchId = busPatchId;
    }
    
    if (!busStruct)
    {
      busStruct = new AliMUONBusStruct;
      busStruct->SetDataKey(busStruct->GetDefaultDataKey());
      busStruct->SetBusPatchId(busPatchId);
      busStruct->SetLength(0);
      busPatchMap.Add(busPatchId,busStruct);
    }
    
    // set sub Event
    busStruct->AddData(word);
  }
}

//______________________________________________________________________________
int AliHLTMUONDigitPublisherComponent::WriteTrackerDDL(
		const AliMUONVDigitStore* digitStore, Int_t iDDL,
		AliHLTUInt8_t* outBuffer, AliHLTUInt32_t& outBufferSize
	)
{
  /// Write DDL file for one tracker DDL
  
  assert(0 <= iDDL and iDDL <= 19);
  
  AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
  assert(ddlStore != NULL);
  
  if (ddlStore->GetDDL(iDDL) == NULL)
  {
  	HLTError("Could not find DDL mapping for DDL %d.", iDDL+1);
  	return -EFAULT;
  }
  
  AliMpExMap busPatchMap;
  Digits2BusPatchMap(*digitStore,busPatchMap,iDDL);
  
  AliMUONBlockHeader blockHeader;
  AliMUONDspHeader dspHeader;
  blockHeader.SetDataKey(blockHeader.GetDefaultDataKey());
  dspHeader.SetDataKey(dspHeader.GetDefaultDataKey());
  
  if (outBufferSize < sizeof(AliRawDataHeader))
  {
  	HLTError("The output buffer size is too small to write output."
  		" It is only %d bytes, but we need at least %d bytes.",
  		outBufferSize, sizeof(AliRawDataHeader)
  	);
  	return -ENOBUFS;
  }
  AliRawDataHeader* header = reinterpret_cast<AliRawDataHeader*>(outBuffer);
  // Fill header with default values.
  *header = AliRawDataHeader();
#ifndef HAVE_NOT_ALIRUNLOADER30859
  AliRunLoader* runloader = AliRunLoader::Instance();
#else
  // the old way before rev 30859
  AliRunLoader *runloader = AliRunLoader::GetRunLoader();
#endif
  if (runloader != NULL)
  {
    if (runloader->GetTrigger() != NULL)
    {
      AliCentralTrigger *aCTP = runloader->GetTrigger();
      ULong64_t mask = aCTP->GetClassMask();
      header->SetTriggerClass(mask);
    }
  }
  
  Int_t* buffer = reinterpret_cast<Int_t*>(header+1);
  Int_t endOfBuffer = (outBufferSize - sizeof(AliRawDataHeader)) / sizeof(Int_t);
  
  // buffer size (max'ed out)
  // (((43 manus max per bus patch *64 channels + 4 bus patch words) * 5 bus patch 
  //   + 10 dsp words)*5 dsps + 8 block words)*2 blocks 
  
  AliMpDDL* ddl = ddlStore->GetDDL(iDDL);
  Int_t iDspMax = ddl->GetMaxDsp();
  Int_t iBusPerDSP[5]; //number of bus patches per DSP
  ddl->GetBusPerDsp(iBusPerDSP);
  Int_t busIter = 0;
  
  Int_t totalDDLLength = 0;
  
  Int_t index = 0;
  
  // two blocks A and B per DDL
  for (Int_t iBlock = 0; iBlock < 2; ++iBlock) 
  {
    Int_t length = blockHeader.GetHeaderLength();
    if (index + length >= endOfBuffer)
    {
      HLTError("The output buffer size is too small to write output."
               " It is only %d bytes, but we need at least %d bytes.",
               outBufferSize,
               sizeof(AliRawDataHeader) + (index+length)*sizeof(UInt_t)
      );
      return -ENOBUFS;
    }
  
    // block header
    memcpy(&buffer[index],blockHeader.GetHeader(),length*4);
    Int_t indexBlk = index;
    index += length; 
    
    // 5 DSP's max per block
    for (Int_t iDsp = 0; iDsp < iDspMax; ++iDsp) 
    {
      Int_t dspHeaderLength = dspHeader.GetHeaderLength();
      if (index + dspHeaderLength >= endOfBuffer)
      {
        HLTError("The output buffer size is too small to write output."
                 " It is only %d bytes, but we need at least %d bytes.",
                 outBufferSize,
                 sizeof(AliRawDataHeader) + (index+dspHeaderLength)*sizeof(UInt_t)
        );
        return -ENOBUFS;
      }
      
      // DSP header
      memcpy(&buffer[index],dspHeader.GetHeader(),dspHeaderLength*4);
      Int_t indexDsp = index;
      index += dspHeaderLength; 
      
      // 5 buspatches max per DSP
      for (Int_t i = 0; i < iBusPerDSP[iDsp]; ++i) 
      {
        Int_t iBusPatch = ddl->GetBusPatchId(busIter++);
        
        // iteration over bus patch in DDL
        if (iBusPatch == -1) 
        {
          AliWarning(Form("Error in bus itr in DDL %d\n", iDDL));
          continue;
        }
        
        AliMUONBusStruct* busStructPtr = static_cast<AliMUONBusStruct*>(busPatchMap.GetValue(iBusPatch));
        
        Int_t busHeaderLength = busStructPtr->GetHeaderLength();
        if (index + busHeaderLength >= endOfBuffer)
        {
          HLTError("The output buffer size is too small to write output."
                   " It is only %d bytes, but we need at least %d bytes.",
                   outBufferSize,
                   sizeof(AliRawDataHeader) + (index+busHeaderLength)*sizeof(UInt_t)
          );
          return -ENOBUFS;
        }
      
        // check if buspatchid has digit
        if (busStructPtr) 
        {
          // add bus patch structure header
          memcpy(&buffer[index],busStructPtr->GetHeader(),busHeaderLength*4);
          index += busHeaderLength;
          
          Int_t busLength = busStructPtr->GetLength();
          if (index + busLength >= endOfBuffer)
          {
            HLTError("The output buffer size is too small to write output."
                     " It is only %d bytes, but we need at least %d bytes.",
                     outBufferSize,
                     sizeof(AliRawDataHeader) + (index+busLength)*sizeof(UInt_t)
            );
            return -ENOBUFS;
          }
          
          // add bus patch data
          memcpy(&buffer[index],busStructPtr->GetData(),busLength*4);
          index += busLength;
        } 
        else 
        {
          // writting anyhow buspatch structure (empty ones)
          buffer[index++] = busStructPtr->GetDefaultDataKey(); // fill it also for empty data size
          buffer[index++] = busStructPtr->GetHeaderLength(); // header length
          buffer[index++] = 0; // raw data length
          buffer[index++] = iBusPatch; // bus patch
        }
      } // bus patch
      
      if (index + 1 >= endOfBuffer)
      {
        HLTError("The output buffer size is too small to write output."
                 " It is only %d bytes, but we need at least %d bytes.",
                 outBufferSize,
                 sizeof(AliRawDataHeader) + (index+1)*sizeof(UInt_t)
        );
        return -ENOBUFS;
      }
      
      // check if totalLength even
      // set padding word in case
      // Add one word 0xBEEFFACE at the end of DSP structure
      Int_t totalDspLength  = index - indexDsp;
      if ((totalDspLength % 2) == 1) 
      { 
        buffer[indexDsp + dspHeader.GetHeaderLength() - 2] = 1;
        buffer[index++] = dspHeader.GetDefaultPaddingWord();
        totalDspLength++;
      }
      
      Int_t dspLength     = totalDspLength - dspHeader.GetHeaderLength();
      
      buffer[indexDsp+1] = totalDspLength; // dsp total length
      buffer[indexDsp+2] = dspLength; // data length  
      
    } // dsp
    
    Int_t totalBlkLength  = index - indexBlk;
    Int_t blkLength       = totalBlkLength - blockHeader.GetHeaderLength();
    totalDDLLength       += totalBlkLength;
    
    buffer[indexBlk+1] = totalBlkLength; // total block length
    buffer[indexBlk+2] = blkLength;
        
  } // block
  
  if (index + 2 >= endOfBuffer)
  {
    HLTError("The output buffer size is too small to write output."
             " It is only %d bytes, but we need at least %d bytes.",
             outBufferSize,
             sizeof(AliRawDataHeader) + (index+2)*sizeof(UInt_t)
    );
    return -ENOBUFS;
  }
  
  // add twice the end of CRT structure data key
  // hope it's good placed (ChF)
  buffer[index++] = blockHeader.GetDdlDataKey();
  buffer[index++] = blockHeader.GetDdlDataKey();
  totalDDLLength  += 2;
  
  header->fSize = (totalDDLLength) * sizeof(Int_t) + sizeof(AliRawDataHeader);
  outBufferSize = header->fSize;
  
  return 0;
}

//______________________________________________________________________________
int AliHLTMUONDigitPublisherComponent::WriteTriggerDDL(
		const AliMUONVTriggerStore* triggerStore, Int_t iDDL,
		AliHLTUInt8_t* outBuffer, AliHLTUInt32_t& outBufferSize,
		bool scalarEvent
	)
{
  /// Write trigger DDL
  
  assert(0 <= iDDL and iDDL <= 19);
  
  AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
  assert(ddlStore != NULL);
  
  if (outBufferSize < sizeof(AliRawDataHeader))
  {
  	HLTError("The output buffer size is too small to write output."
  		" It is only %d bytes, but we need at least %d bytes.",
  		outBufferSize, sizeof(AliRawDataHeader)
  	);
  	return -ENOBUFS;
  }
  AliRawDataHeader* header = reinterpret_cast<AliRawDataHeader*>(outBuffer);
  // Fill header with default values.
  *header = AliRawDataHeader();
#ifndef HAVE_NOT_ALIRUNLOADER30859
  AliRunLoader* runloader = AliRunLoader::Instance();
#else
  // the old way before rev 30859
  AliRunLoader *runloader = AliRunLoader::GetRunLoader();
#endif
  if (runloader != NULL)
  {
    if (runloader->GetTrigger() != NULL)
    {
      AliCentralTrigger *aCTP = runloader->GetTrigger();
      ULong64_t mask = aCTP->GetClassMask();
      header->SetTriggerClass(mask);
    }
  }

  // global trigger for trigger pattern
  AliMUONGlobalTrigger* gloTrg = triggerStore->Global();
  if (!gloTrg) 
  {
    return 0;
  }
  
  Int_t gloTrigResp = gloTrg->GetGlobalResponse();

  UInt_t word;
  Int_t* buffer = reinterpret_cast<Int_t*>(header+1);
  Int_t index;
  Int_t locCard;
  UChar_t locDec, trigY, posY, posX, regOut;
  UInt_t regInpLpt;
  UInt_t regInpHpt;

  UInt_t devX;
  UChar_t sdevX;
  UInt_t version = 1; // software version
  UInt_t eventPhys = 1; // trigger type: 1 for physics, 0 for software
  UInt_t serialNb = 0xF; // serial nb of card: all bits on for the moment
  Int_t globalFlag = 0; // set to 1 if global info present in DDL else set to 0

  AliMUONDarcHeader  darcHeader;
  AliMUONRegHeader   regHeader;
  AliMUONLocalStruct localStruct;

  // size of headers
  static const Int_t kDarcHeaderLength   = darcHeader.GetDarcHeaderLength();
  static const Int_t kGlobalHeaderLength = darcHeader.GetGlobalHeaderLength();
  static const Int_t kDarcScalerLength   = darcHeader.GetDarcScalerLength();
  static const Int_t kGlobalScalerLength = darcHeader.GetGlobalScalerLength();
  static const Int_t kRegHeaderLength    = regHeader.GetHeaderLength();
  static const Int_t kRegScalerLength    = regHeader.GetScalerLength();
  static const Int_t kLocHeaderLength    = localStruct.GetLength();
  static const Int_t kLocScalerLength    = localStruct.GetScalerLength();

  // [16(local)*6 words + 6 words]*8(reg) + 8 words = 824 
  static const Int_t kBufferSize = (16 * (kLocHeaderLength+1) +  (kRegHeaderLength+1))* 8 
      +  kDarcHeaderLength + kGlobalHeaderLength + 2;

  // [16(local)*51 words + 16 words]*8(reg) + 8 + 10 + 8 words scaler event 6682 words
  static const Int_t kScalerBufferSize = (16 * (kLocHeaderLength +  kLocScalerLength +1) +  
					 (kRegHeaderLength + kRegScalerLength +1))* 8 +
                                         (kDarcHeaderLength + kDarcScalerLength + 
					  kGlobalHeaderLength + kGlobalScalerLength + 2);
  if(scalarEvent) {
    eventPhys = 0; //set to generate scaler events
    header->fWord2 |= (0x1 << 14); // set L1SwC bit on
  }
  if(scalarEvent)
  {
    if (outBufferSize < sizeof(AliRawDataHeader) + kScalerBufferSize)
    {
      HLTError("The output buffer size is too small to write output."
               " It is only %d bytes, but we need at least %d bytes.",
               outBufferSize, sizeof(AliRawDataHeader) + kScalerBufferSize
      );
      return -ENOBUFS;
    }
  }
  else
  {
    if (outBufferSize < sizeof(AliRawDataHeader) + kBufferSize)
    {
      HLTError("The output buffer size is too small to write output."
               " It is only %d bytes, but we need at least %d bytes.",
               outBufferSize, sizeof(AliRawDataHeader) + kBufferSize
      );
      return -ENOBUFS;
    }
  }

    index = 0; 

    if (iDDL == 0) // suppose global info in DDL one
      globalFlag = 1;
    else 
      globalFlag = 0;

    word = 0;
    // set darc status word
    // see AliMUONDarcHeader.h for details
    AliBitPacking::PackWord((UInt_t)eventPhys,word,30,30);
    AliBitPacking::PackWord((UInt_t)serialNb,word,20,23);
    AliBitPacking::PackWord((UInt_t)globalFlag,word,10,10);
    AliBitPacking::PackWord((UInt_t)version,word,12,19);
    darcHeader.SetWord(word);

    memcpy(&buffer[index], darcHeader.GetHeader(), (kDarcHeaderLength)*4); 
    index += kDarcHeaderLength;

    // no global input for the moment....
    if (iDDL == 0)
     darcHeader.SetGlobalOutput(gloTrigResp);
    else 
     darcHeader.SetGlobalOutput(0);

    if (scalarEvent) {
      // 6 DARC scaler words
      memcpy(&buffer[index], darcHeader.GetDarcScalers(),kDarcScalerLength*4);
      index += kDarcScalerLength;
    }
    // end of darc word
    buffer[index++] = darcHeader.GetEndOfDarc();

    // 4 words of global board input + Global board output
    memcpy(&buffer[index], darcHeader.GetGlobalInput(), (kGlobalHeaderLength)*4); 
    index += kGlobalHeaderLength; 

    if (scalarEvent) {
      // 10 Global scaler words
      memcpy(darcHeader.GetGlobalScalers(), &buffer[index], kGlobalScalerLength*4);
      index += kGlobalScalerLength;
    }

    // end of global word
    buffer[index++] = darcHeader.GetEndOfGlobal();
    const AliMpRegionalTrigger* reg = ddlStore->GetRegionalTrigger(); 

    Int_t nCrate = reg->GetNofTriggerCrates()/2;
    // 8 regional cards per DDL
    for (Int_t iReg = 0; iReg < nCrate; ++iReg) {

        // crate info
      AliMpTriggerCrate* crate = ddlStore->GetTriggerCrate(iDDL, iReg);

      if (!crate) 
	AliWarning(Form("Missing crate number %d in DDL %d\n", iReg, iDDL));

      // regional info tree, make sure that no reg card missing
      AliMUONRegionalTrigger* regTrg  = triggerStore->FindRegional(crate->GetId());
      if (!regTrg) 
        AliError(Form("Missing regional board %d in trigger Store\n", crate->GetId()));
    
      // Regional card header
      word = 0;

      // set darc status word
      regHeader.SetDarcWord(word);

      regOut    = regTrg->GetOutput();
      regInpHpt = regTrg->GetLocalOutput(0);
      regInpLpt = regTrg->GetLocalOutput(1);

      // fill darc word, not darc status for the moment (empty)
      //see  AliMUONRegHeader.h for details
      AliBitPacking::PackWord((UInt_t)eventPhys,word,31,31); 
      AliBitPacking::PackWord((UInt_t)serialNb,word,20,25); 
      AliBitPacking::PackWord((UInt_t)version,word,8,15);
      AliBitPacking::PackWord((UInt_t)crate->GetId(),word,16,19);
      AliBitPacking::PackWord((UInt_t)regOut,word,0,7);
      regHeader.SetWord(word);


      // fill header later, need local response
      Int_t indexReg = index;
      index += kRegHeaderLength;

      // 11 regional scaler word
      if (scalarEvent) {
	memcpy(&buffer[index], regHeader.GetScalers(), kRegScalerLength*4);
	index += kRegScalerLength;
      }

      // end of regional word
      buffer[index++] = regHeader.GetEndOfReg();
      
      // 16 local card per regional board
      //      UShort_t localMask = 0x0;
      
      Int_t nLocalBoard = AliMpConstants::LocalBoardNofChannels();

      for (Int_t iLoc = 0; iLoc < nLocalBoard; iLoc++) {
	  
	// slot zero for Regional card
	Int_t localBoardId = crate->GetLocalBoardId(iLoc);

	if (localBoardId) { // if not empty slot
	  AliMpLocalBoard* localBoard = ddlStore->GetLocalBoard(localBoardId);

	  if (localBoard->IsNotified()) {// if notified board 
	    AliMUONLocalTrigger* locTrg = triggerStore->FindLocal(localBoardId);

	    locCard = locTrg->LoCircuit();
	    locDec  = locTrg->GetLoDecision();
	    trigY   = locTrg->LoTrigY();
	    posY    = locTrg->LoStripY();
	    posX    = locTrg->LoStripX();
	    devX    = locTrg->LoDev();
	    sdevX   = locTrg->LoSdev();
		  
	    AliDebug(4,Form("loctrg %d, posX %d, posY %d, devX %d\n", 
			    locTrg->LoCircuit(),locTrg->LoStripX(),locTrg->LoStripY(),locTrg->LoDev()));  
	    //packing word
	    word = 0;
	    LocalWordPacking(word, (UInt_t)iLoc, (UInt_t)locDec, (UInt_t)trigY, (UInt_t)posY, 
			     (UInt_t)posX, (UInt_t)sdevX, (UInt_t)devX);

	    buffer[index++] = (locTrg->GetX1Pattern() | (locTrg->GetX2Pattern() << 16));
	    buffer[index++] = (locTrg->GetX3Pattern() | (locTrg->GetX4Pattern() << 16));
	    buffer[index++] = (locTrg->GetY1Pattern() | (locTrg->GetY2Pattern() << 16));
	    buffer[index++] = (locTrg->GetY3Pattern() | (locTrg->GetY4Pattern() << 16));
	    buffer[index++] = (Int_t)word; // data word
		      
		
	  }
	  // fill copy card X-Y inputs from the notified cards 
	  if (localBoard->GetInputXfrom() && localBoard->GetInputYfrom()) 
	  {
	    // not triggered
	    locDec = 0; trigY = 1; posY = 15; 	 
	    posX   = 0; devX  = 0; sdevX = 1;
	    LocalWordPacking(word, (UInt_t)iLoc, (UInt_t)locDec, (UInt_t)trigY, (UInt_t)posY, 
			     (UInt_t)posX, (UInt_t)sdevX, (UInt_t)devX);

	    Int_t localFromId = localBoard->GetInputXfrom();
	    AliMUONLocalTrigger* locTrgfrom  = triggerStore->FindLocal(localFromId);

	    buffer[index++] = 0; // copy only X3-4 & Y1-4
	    buffer[index++] = (locTrgfrom->GetX3Pattern() | (locTrgfrom->GetX4Pattern() << 16));
	    buffer[index++] = (locTrgfrom->GetY1Pattern() | (locTrgfrom->GetY2Pattern() << 16));
	    buffer[index++] = (locTrgfrom->GetY3Pattern() | (locTrgfrom->GetY4Pattern() << 16));
	    buffer[index++] = word;
	  }

	} else { 
	  // fill with 10CDEAD word for empty slots
	  for (Int_t i = 0; i < localStruct.GetLength(); i++)
	      buffer[index++] = localStruct.GetDisableWord(); 
	}// condition localBoard
	  
	// 45 regional scaler word
	if (scalarEvent) {
	  memcpy(&buffer[index], localStruct.GetScalers(), kLocScalerLength*4);
	  index += kLocScalerLength;
	}
	  
	// end of local structure words
	buffer[index++] = localStruct.GetEndOfLocal();
	  
      } // local card 
      // fill regional header with local output
      regHeader.SetInput(regInpHpt, 0);
      regHeader.SetInput(regInpHpt, 1);
      memcpy(&buffer[indexReg],regHeader.GetHeader(),kRegHeaderLength*4);
      
    } // Regional card

  header->fSize = index * sizeof(Int_t) + sizeof(AliRawDataHeader);
  outBufferSize = header->fSize;

  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////
