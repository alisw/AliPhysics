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
#include "AliHLTLogging.h"
#include "AliHLTSystem.h"
#include "AliHLTDefinitions.h"
#include "AliRawDataHeader.h"
#include "AliMUONConstants.h"
#include "AliMUONTrackerDDLDecoderEventHandler.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <cassert>


ClassImp(AliHLTMUONDigitPublisherComponent)


AliHLTMUONDigitPublisherComponent::AliHLTMUONDigitPublisherComponent() :
	AliHLTOfflineDataSource(),
	fDDL(-1),
	fCurrentEventIndex(0)
{
	/// Default constructor.
}


AliHLTMUONDigitPublisherComponent::~AliHLTMUONDigitPublisherComponent()
{
	/// Default destructor.
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


int AliHLTMUONDigitPublisherComponent::DoInit(int argc, const char** argv)
{
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	
	HLTInfo("Initialising dHLT digit publisher component.");

	// Initialise with default values.
	fDDL = -1;
	fCurrentEventIndex = 0;
	bool simdata = false;
	bool recdata = false;
	bool firstEventSet = false;
	bool eventNumLitSet = false;

	for (int i = 0; i < argc; i++)
	{
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
				HLTError("DDL number not specified. It must be in the range [21..22]" );
				return -EINVAL;
			}
		
			char* cpErr = NULL;
			unsigned long num = strtoul(argv[i+1], &cpErr, 0);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a DDL Number.", argv[i+1]);
				return -EINVAL;
			}
			if (num < 1 or 22 < num)
			{
				HLTError("The DDL ID number must be in the range [1..22].");
				return -EINVAL;
			}
			fDDL = num - 1; // Convert to DDL number in the range 0..21
			
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
			char* end = "";
			long num = strtol(argv[i], &end, 0);
			if (*end != '\0' or num < 0) // Check if the conversion is OK.
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
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}
	
	if (fDDL == -1)
	{
		HLTError("DDL number must be set with the -ddl option, but it was not.");
		return -EINVAL;
	}
	
	//TODO initialise data source

	return 0;
}


int AliHLTMUONDigitPublisherComponent::DoDeinit()
{
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	
	HLTInfo("Deinitialising dHLT digit publisher component.");
	return 0;
}


int AliHLTMUONDigitPublisherComponent::GetEvent(
		const AliHLTComponentEventData& /*evtData*/,
		AliHLTComponentTriggerData& /*trigData*/,
		AliHLTUInt8_t* /*outputPtr*/,
		AliHLTUInt32_t& size,
		vector<AliHLTComponentBlockData>& /*outputBlocks*/
	)
{
	/// Inherited from AliHLTOfflineDataSource.
	
	//TODO
	
	size = 0;
	return 0;
}

