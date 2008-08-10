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
///  @file   AliHLTMUONRawDataHistoComponent.cxx
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   30 April 2008
///  @brief  Implementation of the raw data histogramming component for dHLT.
///
/// The class implements 

#include "AliHLTMUONRawDataHistoComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTDataTypes.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliRawDataHeader.h"
#include "TTimeStamp.h"
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cmath>
#include <new>


// Helper type for memory allocation.
typedef const AliHLTMUONMansoTrackStruct* AliHLTMUONMansoTrackStructP;


ClassImp(AliHLTMUONRawDataHistoComponent);


AliHLTMUONRawDataHistoComponent::AliHLTMUONRawDataHistoComponent() :
	AliHLTMUONProcessor(),
	fTrackerDecoder(),
	fTriggerDecoder(),
	fLastPublishTime(-1),
	fCurrentEventTime(-1),
	fPublishDelay(1),
	fSuppressEmptyHists(false),
	fProcessDataEventsOnly(false)
{
	/// Default constructor initialises all histogram object pointers to NULL.
	
	for (int i = 0; i < 22; i++)
	{
		fErrorHist[i] = NULL;
	}
	for (int i = 0; i < 20; i++)
	{
		fManuHist[i] = NULL;
		fSignalHist[i] = NULL;
	}
	
	fTrackerDecoder.ExitOnError(false);
	fTrackerDecoder.TryRecover(false);
	fTrackerDecoder.SendDataOnParityError(true);
	fTrackerDecoder.AutoDetectTrailer(true);
	fTrackerDecoder.CheckForTrailer(true);
	
	fTriggerDecoder.ExitOnError(false);
	fTriggerDecoder.TryRecover(false);
	fTriggerDecoder.AutoDetectScalars(false);
}


AliHLTMUONRawDataHistoComponent::~AliHLTMUONRawDataHistoComponent()
{
	/// Default destructor deletes all histogram objects if they are still allocated.
	
}


const char* AliHLTMUONRawDataHistoComponent::GetComponentID()
{
	/// Inherited from AliHLTComponent. Returns the component ID.
	
	return AliHLTMUONConstants::RawDataHistogrammerId();
}


void AliHLTMUONRawDataHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
	/// Inherited from AliHLTProcessor. Returns the list of expected input data types.
	
	assert( list.empty() );
	list.push_back( AliHLTMUONConstants::DDLRawDataType() );
}


AliHLTComponentDataType AliHLTMUONRawDataHistoComponent::GetOutputDataType()
{
	/// Inherited from AliHLTComponent. Returns kAliHLTHistogramDataTypeID.
	
	return AliHLTMUONConstants::HistogramDataType();
}


void AliHLTMUONRawDataHistoComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	
	constBase = sizeof(TH1D) * 1024*1024;
	inputMultiplier = 0;
}


AliHLTComponent* AliHLTMUONRawDataHistoComponent::Spawn()
{
	/// Inherited from AliHLTComponent. Creates a new object instance.
	
	return new AliHLTMUONRawDataHistoComponent;
}


bool AliHLTMUONRawDataHistoComponent::IgnoreArgument(const char* arg) const
{
	/// Return true if the argument is one of -cdbpath -run or -delaysetup
	/// to prevent the parent class from parsing these arguments in DoInit.
	
	if (strcmp(arg, "-cdbpath") == 0 or strcmp(arg, "-run") == 0 or
	    strcmp(arg, "-delaysetup") == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}


int AliHLTMUONRawDataHistoComponent::DoInit(int argc, const char** argv)
{
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	
	HLTInfo("Initialising dHLT raw data histogrammer component.");
	
	// Inherit the parents functionality.
	int result = AliHLTMUONProcessor::DoInit(argc, argv);
	if (result != 0) return result;

	fLastPublishTime = fCurrentEventTime = -1;
	fPublishDelay = 1;
	bool pubDelaySet = false;
	fSuppressEmptyHists = false;
	fProcessDataEventsOnly = false;
	fTrackerDecoder.TryRecover(false);
	fTriggerDecoder.TryRecover(false);
	
	for (int i = 0; i < argc; i++)
	{
		if (ArgumentAlreadyHandled(i, argv[i])) continue;

		if (strcmp(argv[i], "-pubdelay") == 0)
		{
			if (pubDelaySet)
			{
				HLTWarning("The publishing delay value was already specified."
					" Will replace previous value given by -pubdelay."
				);
			}
			
			if (argc <= i+1)
			{
				HLTError("The value for the publishing delay was not specified.");
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			double num = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0' or num < 0)
			{
				HLTError("Cannot convert '%s' to a positive floating point number.",
					argv[i+1]
				);
				return -EINVAL;
			}
			fPublishDelay = num;
			pubDelaySet = true;
			
			i++;
			continue;
		}
		
		if (strcmp(argv[i], "-noemptyhists") == 0)
		{
			fSuppressEmptyHists = true;
			continue;
		}
		
		if (strcmp(argv[i], "-onlydataevents") == 0)
		{
			fProcessDataEventsOnly = true;
			continue;
		}
		
		if (strcmp(argv[i], "-tryrecover") == 0)
		{
			fTrackerDecoder.TryRecover(true);
			fTriggerDecoder.TryRecover(true);
			continue;
		}

		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}
	
	try
	{
		char name[256];
		char title[1024];
		
		// Do not add to current directory to prevent memory leak warning.
		// We will not be leaking any memory if we dont add to the directory.
		TH1::AddDirectory(kFALSE);
		
		for (int i = 0; i < 22; i++)
		{
			AliHLTInt32_t equipId = AliHLTMUONUtils::DDLNumberToEquipId(i);
			sprintf(name, "rawDataErrors_%d", equipId);
			sprintf(title, "Distribution of errors found in raw data from DDL %d.", equipId);
			fErrorHist[i] = new TH1D(name, title, 40, 0.5, 40.5);
			fErrorHist[i]->SetXTitle("Error code");
			fErrorHist[i]->SetYTitle("Number of errors");
		}
		for (int i = 0; i < 20; i++)
		{
			AliHLTInt32_t equipId = AliHLTMUONUtils::DDLNumberToEquipId(i);
			sprintf(name, "manuDistrib_%d", equipId);
			sprintf(title, "Distribution of MANUs containing raw data in DDL %d.", equipId);
			fManuHist[i] = new TH1D(name, title, 2048, -0.5, 2047.5);
			fManuHist[i]->SetXTitle("MANU number (as seen in raw data)");
			fManuHist[i]->SetYTitle("Number of signals read.");
			sprintf(name, "signalDistrib_%d", equipId);
			sprintf(title, "Distribution of signals in raw data from DDL %d.", equipId);
			fSignalHist[i] = new TH1D(name, title, 4096, -0.5, 4095.5);
			fSignalHist[i]->SetXTitle("Channels");
			fSignalHist[i]->SetYTitle("dN/dChannel");
		}
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for histogram objects.");
		FreeObjects();
		return -ENOMEM;
	}
	
	return 0;
}


int AliHLTMUONRawDataHistoComponent::DoDeinit()
{
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	/// Will delete all histogram objects.
	
	HLTInfo("Deinitialising dHLT raw data histogrammer component.");
	fCurrentEventTime = -1;
	FreeObjects();
	return 0;
}


int AliHLTMUONRawDataHistoComponent::DoEvent(
		const AliHLTComponentEventData& /*evtData*/,
		AliHLTComponentTriggerData& /*trigData*/
	)
{
	/// Inherited from AliHLTProcessor.
	/// Processes the new event data and generates summary histograms.
	
	if (fProcessDataEventsOnly and not IsDataEvent()) return 0;  // Only process data events.
	
	fCurrentEventTime = TTimeStamp().AsDouble();
	
	const AliHLTComponentBlockData* block = GetFirstInputBlock(AliHLTMUONConstants::DDLRawDataType());
	for ( ; block != NULL; block = GetNextInputBlock())
	{
		HLTDebug("Handling block with fDataType = '%s', fPtr = %p,"
			" fSize = %u bytes and fSpecification = 0x%8.8X.",
			DataType2Text(block->fDataType).c_str(), block->fPtr,
			block->fSize, block->fSpecification
		);

		if (AliHLTMUONUtils::IsTrackerDDL(block->fSpecification))
		{
			ProcessTrackerDDL(block);
		}
		else if (AliHLTMUONUtils::IsTriggerDDL(block->fSpecification))
		{
			ProcessTriggerDDL(block);
		}
		else
		{
			HLTError("Received a raw data block with an invalid specification of"
				" 0x%8.8X. Expected raw data only from one DDL and not multiple"
				" DDLs as indicated by the specification.",
				block->fSpecification
			);
		}
	}
	
	// See if 'fPublishDelay' number of seconds has elapsed or this is the first event,
	// in that case publish the histograms. Do not publish histograms that are empty
	// if the fSuppressEmptyHists flag is set.
	if (fLastPublishTime == -1 or fCurrentEventTime - fLastPublishTime >= fPublishDelay)
	{
		for (int i = 0; i < 22; i++)
		{
			if (fSuppressEmptyHists and fErrorHist[i]->GetEntries() == 0) continue;
			PushBack(fErrorHist[i],
				AliHLTMUONConstants::HistogramDataType(),
				AliHLTMUONUtils::DDLNumberToSpec(i)
			);
		}
		for (int i = 0; i < 20; i++)
		{
			AliHLTUInt32_t spec = AliHLTMUONUtils::DDLNumberToSpec(i);
			if (not (fSuppressEmptyHists and fManuHist[i]->GetEntries() == 0))
			{
				PushBack(fManuHist[i], AliHLTMUONConstants::HistogramDataType(), spec);
			}
			if (not (fSuppressEmptyHists and fSignalHist[i]->GetEntries() == 0))
			{
				PushBack(fSignalHist[i], AliHLTMUONConstants::HistogramDataType(), spec);
			}
		}
		fLastPublishTime = fCurrentEventTime;
	}
	
	return 0;
}


void AliHLTMUONRawDataHistoComponent::ProcessTrackerDDL(const AliHLTComponentBlockData* block)
{
	/// Processes a raw data block from the tracker stations.
	
	AliHLTInt32_t ddl = AliHLTMUONUtils::SpecToDDLNumber(block->fSpecification);
	assert(0 <= ddl and ddl < 20);
	
	fTrackerDecoder.GetHandler().ErrorHist(fErrorHist[ddl]);
	fTrackerDecoder.GetHandler().ManuHist(fManuHist[ddl]);
	fTrackerDecoder.GetHandler().SignalHist(fSignalHist[ddl]);
	
	if (block->fSize >= sizeof(AliRawDataHeader))
	{
		AliHLTUInt8_t* payload = reinterpret_cast<AliHLTUInt8_t*>(block->fPtr)
			+ sizeof(AliRawDataHeader);
		UInt_t payloadSize = UInt_t(block->fSize) - sizeof(AliRawDataHeader);
		fTrackerDecoder.Decode(payload, payloadSize);
	}
	else
	{
		HLTError("Received a raw data block that is too short to be valid."
			" Its size is only %d bytes",
			block->fSize
		);
		fErrorHist[ddl]->Fill(40);
	}
}


void AliHLTMUONRawDataHistoComponent::ProcessTriggerDDL(const AliHLTComponentBlockData* block)
{
	/// Processes a raw data block from the trigger stations.
	
	AliHLTInt32_t ddl = AliHLTMUONUtils::SpecToDDLNumber(block->fSpecification);
	assert(21 <= ddl and ddl < 22);
	
	fTriggerDecoder.GetHandler().ErrorHist(fErrorHist[ddl]);
	
	if (block->fSize >= sizeof(AliRawDataHeader))
	{
		AliRawDataHeader* header = reinterpret_cast<AliRawDataHeader*>(block->fPtr);
		AliHLTUInt8_t* payload = reinterpret_cast<AliHLTUInt8_t*>(header+1);
		UInt_t payloadSize = UInt_t(block->fSize) - sizeof(AliRawDataHeader);
		bool scalarEvent = ((header->GetL1TriggerMessage() & 0x1) == 0x1);
		fTriggerDecoder.Decode(payload, payloadSize, scalarEvent);
	}
	else
	{
		HLTError("Received a raw data block that is too short to be valid."
			" Its size is only %d bytes",
			block->fSize
		);
		fErrorHist[ddl]->Fill(40);
	}
}


void AliHLTMUONRawDataHistoComponent::FreeObjects()
{
	/// Deletes all the histogram objects that were allocated.
	
	for (int i = 0; i < 22; i++)
	{
		if (fErrorHist[i] != NULL)
		{
			delete fErrorHist[i];
			fErrorHist[i] = NULL;
		}
	}
	for (int i = 0; i < 20; i++)
	{
		if (fManuHist[i] != NULL)
		{
			delete fManuHist[i];
			fManuHist[i] = NULL;
		}
		if (fSignalHist[i] != NULL)
		{
			delete fSignalHist[i];
			fSignalHist[i] = NULL;
		}
	}
}

