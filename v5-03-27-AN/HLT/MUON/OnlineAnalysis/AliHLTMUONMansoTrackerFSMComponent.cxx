/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$

///
///  @file   AliHLTMUONMansoTrackerFSMComponent.cxx
///  @author Artur Szostak <artursz@iafrica.com>,
///          Indranil Das <indra.das@saha.ac.in>
///  @date   18 Sep 2007
///  @brief  Implementation of AliHLTMUONMansoTrackerFSMComponent class.
///

#include "AliHLTMUONMansoTrackerFSMComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONCalculations.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONMansoTrackerFSM.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONDataBlockWriter.h"
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <new>

ClassImp(AliHLTMUONMansoTrackerFSMComponent);


AliHLTMUONMansoTrackerFSMComponent::AliHLTMUONMansoTrackerFSMComponent() :
	AliHLTMUONProcessor(),
	AliHLTMUONMansoTrackerFSMCallback(),
	fTracker(NULL),
	fTrackCount(0),
	fBlock(NULL),
	fRecHitBlockArraySize(0),
	fWarnForUnexpecedBlock(false),
	fCanLoadZmiddle(true),
	fCanLoadBL(true)
{
	///
	/// Default constructor.
	///
	
	for (int i = 0; i < 4; i++)
	{
		fRecHitBlockCount[i] = 0;
		fRecHitBlock[i] = NULL;
	}
	
	ResetCanLoadFlags();
}


AliHLTMUONMansoTrackerFSMComponent::~AliHLTMUONMansoTrackerFSMComponent()
{
	///
	/// Default destructor.
	///
	
	// Should never have the following 2 pointers non-NULL since DoDeinit
	// should have been called before, but handle this case anyway.
	if (fTracker != NULL) delete fTracker;
	
	// Remember that only fRecHitBlock[0] stores the pointer to the allocated
	// memory. The other pointers are just reletive to this.
	if (fRecHitBlock[0] != NULL) delete [] fRecHitBlock[0];
}


const char* AliHLTMUONMansoTrackerFSMComponent::GetComponentID()
{
	///
	/// Inherited from AliHLTComponent. Returns the component ID.
	///
	
	return AliHLTMUONConstants::MansoTrackerFSMId();
}


void AliHLTMUONMansoTrackerFSMComponent::GetInputDataTypes(
		AliHLTComponentDataTypeList& list
	)
{
	///
	/// Inherited from AliHLTProcessor. Returns the list of expected input data types.
	///
	
	assert( list.empty() );
	list.push_back( AliHLTMUONConstants::TriggerRecordsBlockDataType() );
	list.push_back( AliHLTMUONConstants::RecHitsBlockDataType() );
}


AliHLTComponentDataType AliHLTMUONMansoTrackerFSMComponent::GetOutputDataType()
{
	///
	/// Inherited from AliHLTComponent. Returns the output data type.
	///
	
	return kAliHLTMultipleDataType;
}


int AliHLTMUONMansoTrackerFSMComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
	/// Inherited from AliHLTComponent. Returns the output data types.
	
	assert( list.empty() );
	list.push_back( AliHLTMUONConstants::MansoTracksBlockDataType() );
	list.push_back( AliHLTMUONConstants::MansoCandidatesBlockDataType() );
	return list.size();
}


void AliHLTMUONMansoTrackerFSMComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	///
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	///
	
	constBase = sizeof(AliHLTMUONMansoTracksBlockStruct) + 1024*1024;
	inputMultiplier = 1;
}


AliHLTComponent* AliHLTMUONMansoTrackerFSMComponent::Spawn()
{
	///
	/// Inherited from AliHLTComponent. Creates a new object instance.
	///
	
	return new AliHLTMUONMansoTrackerFSMComponent;
}


int AliHLTMUONMansoTrackerFSMComponent::DoInit(int argc, const char** argv)
{
	///
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	///
	
	HLTInfo("Initialising dHLT manso tracker FSM component.");
	
	// Inherit the parents functionality.
	int result = AliHLTMUONProcessor::DoInit(argc, argv);
	if (result != 0) return result;

	// Just in case for whatever reason we still have some of the internal
	// object allocated previously still hanging around delete them now.
	FreeMemory();
	
	fWarnForUnexpecedBlock = false;
	bool makeCandidates = false;
	ResetCanLoadFlags();
	double zmiddle = 0;
	double bfieldintegral = 0;
	double roiA[4] = {0, 0, 0, 0};
	double roiB[4] = {0, 0, 0, 0};
	double chamberZ[6] = {0, 0, 0, 0, 0, 0};
	
	for (int i = 0; i < argc; i++)
	{
		if (ArgumentAlreadyHandled(i, argv[i])) continue;

		if (strcmp( argv[i], "-zmiddle" ) == 0)
		{
			if (not fCanLoadZmiddle)
			{
				HLTWarning("The Z coordinate for the middle of the dipole was already specified."
					" Will replace previous value given by -zmiddle."
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("The Z coordinate for the middle of the dipole was not specified." );
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			zmiddle = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a valid floating point number.",
					argv[i+1]
				);
				return -EINVAL;
			}
			
			fCanLoadZmiddle = false;  // Prevent loading from CDB.
			i++;
			continue;
		}
	
		if (strcmp( argv[i], "-bfieldintegral" ) == 0)
		{
			if (not fCanLoadBL)
			{
				HLTWarning("The magnetic field integral was already specified."
					" Will replace previous value given by -bfieldintegral."
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("The magnetic field integral was not specified." );
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			bfieldintegral = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a valid floating point number.",
					argv[i+1]
				);
				return -EINVAL;
			}
			
			fCanLoadBL = false;  // Prevent loading from CDB.
			i++;
			continue;
		}
	
		if (strcmp(argv[i], "-a7") == 0 or strcmp(argv[i], "-a8") == 0 or
		    strcmp(argv[i], "-a9") == 0 or strcmp(argv[i], "-a10") == 0
		   )
		{
			int chamber = 7; int chamberIndex = 0;
			switch (argv[i][2])
			{
			case '7': chamber = 7; chamberIndex = 0; break;
			case '8': chamber = 8; chamberIndex = 1; break;
			case '9': chamber = 9; chamberIndex = 2; break;
			case '1': chamber = 10; chamberIndex = 3; break;
			}
			
			if (not fCanLoadA[chamberIndex])
			{
				HLTWarning("The region of interest parameter 'A' for chamber %d was"
					" already specified. Will replace previous value given by -a%d.",
					chamber, chamber
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("The region of interest parameter was not specified." );
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			roiA[chamberIndex] = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a valid floating point number.",
					argv[i+1]
				);
				return -EINVAL;
			}
			
			fCanLoadA[chamberIndex] = false;  // Prevent loading from CDB.
			i++;
			continue;
		}
	
		if (strcmp(argv[i], "-b7") == 0 or strcmp(argv[i], "-b8") == 0 or
		    strcmp(argv[i], "-b9") == 0 or strcmp(argv[i], "-b10") == 0
		   )
		{
			int chamber = 7; int chamberIndex = 0;
			switch (argv[i][2])
			{
			case '7': chamber = 7; chamberIndex = 0; break;
			case '8': chamber = 8; chamberIndex = 1; break;
			case '9': chamber = 9; chamberIndex = 2; break;
			case '1': chamber = 10; chamberIndex = 3; break;
			}
			
			if (not fCanLoadB[chamberIndex])
			{
				HLTWarning("The region of interest parameter 'B' for chamber %d was"
					" already specified. Will replace previous value given by -b%d.",
					chamber, chamber
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("The region of interest parameter was not specified." );
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			roiB[chamberIndex] = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a valid floating point number.",
					argv[i+1]
				);
				return -EINVAL;
			}
			
			fCanLoadB[chamberIndex] = false;  // Prevent loading from CDB.
			i++;
			continue;
		}
		
		if (strcmp(argv[i], "-z7") == 0 or strcmp(argv[i], "-z8") == 0 or
		    strcmp(argv[i], "-z9") == 0 or strcmp(argv[i], "-z10") == 0 or
		    strcmp(argv[i], "-z11") == 0 or strcmp(argv[i], "-z13") == 0
		   )
		{
			int chamber = 7; int chamberIndex = 0;
			switch (argv[i][2])
			{
			case '7': chamber = 7; chamberIndex = 0; break;
			case '8': chamber = 8; chamberIndex = 1; break;
			case '9': chamber = 9; chamberIndex = 2; break;
			case '1':
				switch (argv[i][3])
				{
				case '0': chamber = 10; chamberIndex = 3; break;
				case '1': chamber = 11; chamberIndex = 4; break;
				case '3': chamber = 13; chamberIndex = 5; break;
				}
				break;
			}
			
			if (not fCanLoadZ[chamberIndex])
			{
				HLTWarning("The nominal Z coordinate of chamber %d was already"
					" specified. Will replace previous value given by -z%d.",
					chamber, chamber
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("The region of interest parameter was not specified." );
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			chamberZ[chamberIndex] = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a valid floating point number.",
					argv[i+1]
				);
				return -EINVAL;
			}
			
			fCanLoadZ[chamberIndex] = false;  // Prevent loading from CDB.
			i++;
			continue;
		}
		
		if (strcmp(argv[i], "-warn_on_unexpected_block") == 0)
		{
			fWarnForUnexpecedBlock = true;
			continue;
		}
		
		if (strcmp(argv[i], "-makecandidates") == 0)
		{
			makeCandidates = true;
			continue;
		}
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}
	
	try
	{
		fTracker = new AliHLTMUONMansoTrackerFSM();
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the tracker component.");
		return -ENOMEM;
	}
	fTracker->SetCallback(this);
	fTracker->MakeCandidates(makeCandidates);
	
	// Set all the parameters that were found on the command line.
	if (not fCanLoadZmiddle) AliHLTMUONCalculations::Zf(zmiddle);
	if (not fCanLoadBL) AliHLTMUONCalculations::QBL(bfieldintegral);
	if (not fCanLoadA[0]) fTracker->SetA7(roiA[0]);
	if (not fCanLoadA[1]) fTracker->SetA8(roiA[1]);
	if (not fCanLoadA[2]) fTracker->SetA9(roiA[2]);
	if (not fCanLoadA[3]) fTracker->SetA10(roiA[3]);
	if (not fCanLoadB[0]) fTracker->SetB7(roiB[0]);
	if (not fCanLoadB[1]) fTracker->SetB8(roiB[1]);
	if (not fCanLoadB[2]) fTracker->SetB9(roiB[2]);
	if (not fCanLoadB[3]) fTracker->SetB10(roiB[3]);
	if (not fCanLoadZ[0]) fTracker->SetZ7(chamberZ[0]);
	if (not fCanLoadZ[1]) fTracker->SetZ8(chamberZ[1]);
	if (not fCanLoadZ[2]) fTracker->SetZ9(chamberZ[2]);
	if (not fCanLoadZ[3]) fTracker->SetZ10(chamberZ[3]);
	if (not fCanLoadZ[4]) fTracker->SetZ11(chamberZ[4]);
	if (not fCanLoadZ[5]) fTracker->SetZ13(chamberZ[5]);
	
	if (not DelaySetup())
	{
		if (AtLeastOneCanLoadFlagsIsSet())
		{
			HLTInfo("Loading configuration parameters from CDB.");
			
			result = ReadConfigFromCDB();
			if (result != 0)
			{
				// Error messages already generated in ReadConfigFromCDB.
				FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
				return result;
			}
		}
		else
		{
			// Print the debug messages here since ReadConfigFromCDB does not get called,
			// in-which the debug messages would have been printed.
			HLTDebug("Using the following configuration parameters:");
			HLTDebug("                    Middle of dipole Z coordinate = %f cm", AliHLTMUONCalculations::Zf());
			HLTDebug("                          Magnetic field integral = %f T.m", AliHLTMUONCalculations::QBL());
			HLTDebug("   Region of interest parameter 'A' for chamber 7 = %f",    fTracker->GetA7());
			HLTDebug("   Region of interest parameter 'B' for chamber 7 = %f cm", fTracker->GetB7());
			HLTDebug("   Region of interest parameter 'A' for chamber 8 = %f",    fTracker->GetA8());
			HLTDebug("   Region of interest parameter 'B' for chamber 8 = %f cm", fTracker->GetB8());
			HLTDebug("   Region of interest parameter 'A' for chamber 9 = %f",    fTracker->GetA9());
			HLTDebug("   Region of interest parameter 'B' for chamber 9 = %f cm", fTracker->GetB9());
			HLTDebug("  Region of interest parameter 'A' for chamber 10 = %f",    fTracker->GetA10());
			HLTDebug("  Region of interest parameter 'B' for chamber 10 = %f cm", fTracker->GetB10());
			HLTDebug("               Nominal Z coordinate for chamber 7 = %f cm", fTracker->GetZ7());
			HLTDebug("               Nominal Z coordinate for chamber 8 = %f cm", fTracker->GetZ8());
			HLTDebug("               Nominal Z coordinate for chamber 9 = %f cm", fTracker->GetZ9());
			HLTDebug("              Nominal Z coordinate for chamber 10 = %f cm", fTracker->GetZ10());
			HLTDebug("              Nominal Z coordinate for chamber 11 = %f cm", fTracker->GetZ11());
			HLTDebug("              Nominal Z coordinate for chamber 13 = %f cm", fTracker->GetZ13());
		}
		
		ResetCanLoadFlags();  // From this point read all parameters from CDB.
	}
	
	const int initArraySize = 10;
	// Allocate some initial memory for the reconstructed hit arrays.
	try
	{
		fRecHitBlock[0] = new AliRecHitBlockInfo[initArraySize*4];
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the reconstructed hit arrays.");
		FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
		return -ENOMEM;
	}
	// Only set the arrays' size once we have successfully allocated the memory for the arrays.
	fRecHitBlockArraySize = initArraySize;
	// Now we need to set the pointers fRecHitBlock[i] {i>0} relative to fRecHitBlock[0].
	for (Int_t i = 1; i < 4; i++)
	{
		fRecHitBlock[i] = fRecHitBlock[i-1] + fRecHitBlockArraySize;
	}
	// And reset the number of records actually stored in the arrays.
	for (Int_t i = 0; i < 4; i++)
	{
		fRecHitBlockCount[i] = 0;
	}
	
	return 0;
}


int AliHLTMUONMansoTrackerFSMComponent::Reconfigure(
		const char* cdbEntry, const char* componentId
	)
{
	/// Inherited from AliHLTComponent. This method will reload CDB configuration
	/// entries for this component from the CDB.
	/// \param cdbEntry If this is NULL or equals "HLT/ConfigMUON/MansoTrackerFSM"
	///     then new configuration parameters are loaded, otherwise nothing is done.
	/// \param componentId  The name of the component in the current chain.
	
	TString path = cdbEntry;
	bool givenConfigPath = (path == AliHLTMUONConstants::MansoTrackerFSMCDBPath());
	
	if (cdbEntry == NULL or givenConfigPath)
	{
		HLTInfo("Reading new configuration entries from CDB for component '%s'.", componentId);
		ResetCanLoadFlags();  // Make sure to allow reading all values from CDB.
		int result = ReadConfigFromCDB();
		if (result != 0) return result;
	}
	
	return 0;
}


int AliHLTMUONMansoTrackerFSMComponent::ReadPreprocessorValues(const char* modules)
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
		return Reconfigure(AliHLTMUONConstants::MansoTrackerFSMCDBPath(), GetComponentID());
	}
	return 0;
}


int AliHLTMUONMansoTrackerFSMComponent::ReadConfigFromCDB()
{
	/// Reads this component's configuration parameters from the CDB.
	/// These include the middle of the dipole Z coordinate (zmiddle), the
	/// integrated magnetic field of the dipole, Z coordinates of the chambers
	/// and the region of interest parameters used during the tracking.
	/// \param setZmiddle Indicates if the zmiddle parameter should be set
	///       (default true).
	/// \param setBL Indicates if the integrated magnetic field parameter should
	///       be set (default true).
	/// \return 0 if no errors occured and negative error code compatible with
	///       the HLT framework on errors.

	const char* pathToEntry = AliHLTMUONConstants::MansoTrackerFSMCDBPath();
	
	TMap* map = NULL;
	int result = FetchTMapFromCDB(pathToEntry, map);
	if (result != 0) return result;
	
	Double_t value = 0;
	if (fCanLoadZmiddle)
	{
		result = GetFloatFromTMap(map, "zmiddle", value, pathToEntry, "dipole middle Z coordinate");
		if (result != 0) return result;
		AliHLTMUONCalculations::Zf(value);
	}
	
	if (fCanLoadBL)
	{
		Double_t bfieldintegral;
		result = FetchFieldIntegral(bfieldintegral);
		if (result == 0)
		{
			AliHLTMUONCalculations::QBL(bfieldintegral);
		}
		else
		{
			HLTWarning("Failed to load the magnetic field integral from GRP information.");
			result = GetFloatFromTMap(map, "bfieldintegral", value, pathToEntry, "integrated magnetic field");
			if (result != 0) return result;
			HLTWarning(Form("Using deprecated magnetic field integral value of %f T.m.", value));
			AliHLTMUONCalculations::QBL(value);
		}
	}
	
	if (fCanLoadA[0])
	{
		result = GetFloatFromTMap(map, "roi_paramA_chamber7", value, pathToEntry, "chamber 7 region of interest 'A'");
		if (result != 0) return result;
		fTracker->SetA7(value);
	}
	if (fCanLoadA[1])
	{
		result = GetFloatFromTMap(map, "roi_paramA_chamber8", value, pathToEntry, "chamber 8 region of interest 'A'");
		if (result != 0) return result;
		fTracker->SetA8(value);
	}
	if (fCanLoadA[2])
	{
		result = GetFloatFromTMap(map, "roi_paramA_chamber9", value, pathToEntry, "chamber 9 region of interest 'A'");
		if (result != 0) return result;
		fTracker->SetA9(value);
	}
	if (fCanLoadA[3])
	{
		result = GetFloatFromTMap(map, "roi_paramA_chamber10", value, pathToEntry, "chamber 10 region of interest 'A'");
		if (result != 0) return result;
		fTracker->SetA10(value);
	}
	
	if (fCanLoadB[0])
	{
		result = GetFloatFromTMap(map, "roi_paramB_chamber7", value, pathToEntry, "chamber 7 region of interest 'B'");
		if (result != 0) return result;
		fTracker->SetB7(value);
	}
	if (fCanLoadB[1])
	{
		result = GetFloatFromTMap(map, "roi_paramB_chamber8", value, pathToEntry, "chamber 8 region of interest 'B'");
		if (result != 0) return result;
		fTracker->SetB8(value);
	}
	if (fCanLoadB[2])
	{
		result = GetFloatFromTMap(map, "roi_paramB_chamber9", value, pathToEntry, "chamber 9 region of interest 'B'");
		if (result != 0) return result;
		fTracker->SetB9(value);
	}
	if (fCanLoadB[3])
	{
		result = GetFloatFromTMap(map, "roi_paramB_chamber10", value, pathToEntry, "chamber 10 region of interest 'B'");
		if (result != 0) return result;
		fTracker->SetB10(value);
	}
	
	if (fCanLoadZ[0])
	{
		result = GetFloatFromTMap(map, "chamber7postion", value, pathToEntry, "nominal chamber 7 Z coordinate");
		if (result != 0) return result;
		fTracker->SetZ7(value);
	}
	if (fCanLoadZ[1])
	{
		result = GetFloatFromTMap(map, "chamber8postion", value, pathToEntry, "nominal chamber 8 Z coordinate");
		if (result != 0) return result;
		fTracker->SetZ8(value);
	}
	if (fCanLoadZ[2])
	{
		result = GetFloatFromTMap(map, "chamber9postion", value, pathToEntry, "nominal chamber 9 Z coordinate");
		if (result != 0) return result;
		fTracker->SetZ9(value);
	}
	if (fCanLoadZ[3])
	{
		result = GetFloatFromTMap(map, "chamber10postion", value, pathToEntry, "nominal chamber 10 Z coordinate");
		if (result != 0) return result;
		fTracker->SetZ10(value);
	}
	if (fCanLoadZ[4])
	{
		result = GetFloatFromTMap(map, "chamber11postion", value, pathToEntry, "nominal chamber 11 Z coordinate");
		if (result != 0) return result;
		fTracker->SetZ11(value);
	}
	if (fCanLoadZ[5])
	{
		result = GetFloatFromTMap(map, "chamber13postion", value, pathToEntry, "nominal chamber 13 Z coordinate");
		if (result != 0) return result;
		fTracker->SetZ13(value);
	}
	
	HLTDebug("Using the following configuration parameters:");
	HLTDebug("                    Middle of dipole Z coordinate = %f cm", AliHLTMUONCalculations::Zf());
	HLTDebug("                          Magnetic field integral = %f T.m", AliHLTMUONCalculations::QBL());
	HLTDebug("   Region of interest parameter 'A' for chamber 7 = %f",    fTracker->GetA7());
	HLTDebug("   Region of interest parameter 'B' for chamber 7 = %f cm", fTracker->GetB7());
	HLTDebug("   Region of interest parameter 'A' for chamber 8 = %f",    fTracker->GetA8());
	HLTDebug("   Region of interest parameter 'B' for chamber 8 = %f cm", fTracker->GetB8());
	HLTDebug("   Region of interest parameter 'A' for chamber 9 = %f",    fTracker->GetA9());
	HLTDebug("   Region of interest parameter 'B' for chamber 9 = %f cm", fTracker->GetB9());
	HLTDebug("  Region of interest parameter 'A' for chamber 10 = %f",    fTracker->GetA10());
	HLTDebug("  Region of interest parameter 'B' for chamber 10 = %f cm", fTracker->GetB10());
	HLTDebug("               Nominal Z coordinate for chamber 7 = %f cm", fTracker->GetZ7());
	HLTDebug("               Nominal Z coordinate for chamber 8 = %f cm", fTracker->GetZ8());
	HLTDebug("               Nominal Z coordinate for chamber 9 = %f cm", fTracker->GetZ9());
	HLTDebug("              Nominal Z coordinate for chamber 10 = %f cm", fTracker->GetZ10());
	HLTDebug("              Nominal Z coordinate for chamber 11 = %f cm", fTracker->GetZ11());
	HLTDebug("              Nominal Z coordinate for chamber 13 = %f cm", fTracker->GetZ13());
	
	return 0;
}


int AliHLTMUONMansoTrackerFSMComponent::DoDeinit()
{
	///
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	///
	
	HLTInfo("Deinitialising dHLT manso tracker FSM component.");
	FreeMemory();
	return 0;
}


int AliHLTMUONMansoTrackerFSMComponent::DoEvent(
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
	
	// Initialise the configuration parameters from CDB if we were
	// requested to initialise only when the first event was received.
	if (DelaySetup())
	{
		// Load the configuration paramters from CDB if they have not
		// been given on the command line.
		if (AtLeastOneCanLoadFlagsIsSet())
		{
			HLTInfo("Loading configuration parameters from CDB.");
			int result = ReadConfigFromCDB();
			if (result != 0) return result;
		}
		
		DoneDelayedSetup();
		ResetCanLoadFlags();  // From this point read all parameters from CDB.
	}
	
	Reset();
	AliHLTUInt32_t specification = 0;  // Contains the output data block spec bits.
	
	// Resize the rec hit arrays if we possibly will need more space.
	// To guarantee that they will not overflow we need to make sure each
	// array is at least as big as the number of input data blocks.
	if (fRecHitBlockArraySize < evtData.fBlockCnt)
	{
		// Release the old memory block and allocate more memory.
		if (fRecHitBlock[0] != NULL)
		{
			delete [] fRecHitBlock[0];
		}
		
		// Reset the number of records actually stored in the arrays.
		for (Int_t i = 0; i < 4; i++)
		{
			fRecHitBlockCount[i] = 0;
		}
		
		try
		{
			fRecHitBlock[0] = new AliRecHitBlockInfo[evtData.fBlockCnt*4];
		}
		catch (const std::bad_alloc&)
		{
			HLTError("Could not allocate more memory for the reconstructed hit arrays.");
			// Ok so now we need to clear all the pointers because we actually
			// deleted the memory.
			fRecHitBlockArraySize = 0;
			for (Int_t i = 0; i < 4; i++)
			{
				fRecHitBlock[i] = NULL;
			}
			if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
			return -ENOMEM;
		}
		// Only set the arrays' size once we have successfully allocated the memory for the arrays.
		fRecHitBlockArraySize = evtData.fBlockCnt;
		// Now we need to set the pointers fRecHitBlock[i] {i>0} relative to fRecHitBlock[0].
		for (Int_t i = 1; i < 4; i++)
		{
			fRecHitBlock[i] = fRecHitBlock[i-1] + fRecHitBlockArraySize;
		}
	}
	
	AliHLTMUONMansoTracksBlockWriter block(outputPtr, size);
	fBlock = &block;
	
	if (not block.InitCommonHeader())
	{
		Logging(kHLTLogError,
			"AliHLTMUONMansoTrackerFSMComponent::DoEvent",
			"Buffer overflow",
			"The buffer is only %d bytes in size. We need a minimum of %d bytes.",
			size, sizeof(AliHLTMUONMansoTracksBlockWriter::HeaderType)
		);
		if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
		size = 0; // Important to tell framework that nothing was generated.
		return -ENOBUFS;
	}

	// Loop over all input blocks in the event and add the ones that contain
	// reconstructed hits into the hit buffers. The blocks containing trigger
	// records are ignored for now and will be processed later.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
		HLTDebug("Handling block: %u, with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
			n, DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fPtr, blocks[n].fSize
		);
		
		if (blocks[n].fDataType == AliHLTMUONConstants::RecHitsBlockDataType())
		{
			specification |= blocks[n].fSpecification;
			
			AliHLTMUONRecHitsBlockReader inblock(blocks[n].fPtr, blocks[n].fSize);
			if (not BlockStructureOk(inblock))
			{
				if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
				continue;
			}
			
			if (inblock.Nentries() != 0)
				AddRecHits(blocks[n].fSpecification, inblock.GetArray(), inblock.Nentries());
			else
			{
				Logging(kHLTLogDebug,
					"AliHLTMUONMansoTrackerFSMComponent::DoEvent",
					"Block empty",
					"Received a reconstructed hits data block which contains no entries."
				);
			}
		}
		else if (blocks[n].fDataType != AliHLTMUONConstants::TriggerRecordsBlockDataType())
		{
			// Log a message indicating that we got a data block that we
			// do not know how to handle.
			if (fWarnForUnexpecedBlock)
				HLTWarning("Received a data block of a type we cannot handle: '%s', spec: 0x%X",
					DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fSpecification
				);
#ifdef __DEBUG
			else
				HLTDebug("Received a data block of a type we cannot handle: '%s', spec: 0x%X",
					DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fSpecification
				);
#endif
		}
	}
  
	// Again loop over all input blocks in the event, but this time look for
	// the trigger record blocks and process these.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
		if (blocks[n].fDataType != AliHLTMUONConstants::TriggerRecordsBlockDataType())
			continue;
		
		AliHLTMUONTriggerRecordsBlockReader inblock(blocks[n].fPtr, blocks[n].fSize);
		if (not BlockStructureOk(inblock))
		{
			if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
			continue;
		}
		
		DebugTrace("Processing a trigger block with "
			<< inblock.Nentries() << " entries."
		);
		
		specification |= blocks[n].fSpecification;
		
		for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
		{
			fTracker->FindTrack(inblock[i]);
			
			// Reset the tracker so that we do not double count tracks.
			fTracker->Reset();
		}
	}
	
	AliHLTComponentBlockData bd;
	FillBlockData(bd);
	bd.fPtr = outputPtr;
	bd.fOffset = 0;
	bd.fSize = block.BytesUsed();
	bd.fDataType = AliHLTMUONConstants::MansoTracksBlockDataType();
	bd.fSpecification = specification;
	outputBlocks.push_back(bd);
	AliHLTUInt32_t totalSize = block.BytesUsed();
	
	if (fTracker->MakeCandidates())
	{
		AliHLTMUONMansoCandidatesBlockWriter candidatesBlock(outputPtr+totalSize, size-totalSize);
		if (not candidatesBlock.InitCommonHeader())
		{
			HLTError("Buffer overflowed. There are only %d bytes left in the buffer,"
				" but we need a minimum of %d bytes.",
				size, sizeof(AliHLTMUONMansoCandidatesBlockWriter::HeaderType)
			);
			if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
			size = 0; // Important to tell framework that nothing was generated.
			return -ENOBUFS;
		}
		
		// Fill in the output block buffer.
		candidatesBlock.SetNumberOfEntries(fTracker->TrackCandidatesCount());
		for (AliHLTUInt32_t i = 0; i < fTracker->TrackCandidatesCount(); ++i)
		{
			candidatesBlock[i] = fTracker->TrackCandidates()[i];
		}
		
		fTracker->ZeroTrackCandidatesList();
		
		AliHLTComponentBlockData bdc;
		FillBlockData(bdc);
		bdc.fPtr = outputPtr;
		bdc.fOffset = totalSize;
		bdc.fSize = candidatesBlock.BytesUsed();
		bdc.fDataType = AliHLTMUONConstants::MansoCandidatesBlockDataType();
		bdc.fSpecification = specification;
		outputBlocks.push_back(bdc);
		totalSize += candidatesBlock.BytesUsed();
	}
	
	size = totalSize;
	return 0;
}


void AliHLTMUONMansoTrackerFSMComponent::Reset()
{
	///
	/// Reset the track count and reconstructed hit data block arrays.
	///
	
	DebugTrace("Resetting AliHLTMUONMansoTrackerFSMComponent.");

	//fTracker->Reset();  // Not necessary here because it is done after every FindTrack call.
	fTrackCount = 0;
	fBlock = NULL;  // Do not delete. Already done implicitly at the end of DoEvent.
	for (int i = 0; i < 4; i++)
	{
		fRecHitBlockCount[i] = 0;
	}
}


void AliHLTMUONMansoTrackerFSMComponent::FreeMemory()
{
	/// Deletes any objects and arrays allocated by this component and releases
	/// the memory used. This is called as a helper routine by the init and deinit
	/// methods. If some or all of the object pointers are already NULL then
	/// nothing is done for those. This method guarantees that all the relevant
	/// pointers will be NULL after returning from this method.

	if (fTracker != NULL)
	{
		delete fTracker;
		fTracker = NULL;
	}
	
	// Remember that only fRecHitBlock[0] stores the pointer to the allocated memory.
	// The other pointers are just reletive to this.
	if (fRecHitBlock[0] != NULL)
		delete [] fRecHitBlock[0];
	
	fRecHitBlockArraySize = 0;
	for (Int_t i = 0; i < 4; i++)
	{
		fRecHitBlockCount[i] = 0;
		fRecHitBlock[i] = NULL;
	}
}


void AliHLTMUONMansoTrackerFSMComponent::AddRecHits(
		AliHLTUInt32_t specification,
		const AliHLTMUONRecHitStruct* recHits,
		AliHLTUInt32_t count
	)
{
	///
	/// Adds a new reconstructed hit data block to the internal list of blocks
	/// for the tracker to process.
	/// These lists will later be used when the tracker requests them through
	/// the callback method 'RequestClusters'.
	///
	
	DebugTrace("AliHLTMUONMansoTrackerFSMComponent::AddRecHits called with specification = 0x"
		 << std::hex << specification << std::dec << " and count = "
		 << count << " rec hits."
	);
	
	AliHLTUInt8_t chamberMap[20] = {
		1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10
	};
	
	// Identify the chamber the rec hits came from using the specifications field.
	bool gotDataFromDDL[22];
	AliHLTMUONUtils::UnpackSpecBits(specification, gotDataFromDDL);
		
	AliHLTInt8_t chamber = -1;
	for (int i = 0; i < 20; i++)
	{
		if (not gotDataFromDDL[i]) continue;
		if (7 <= chamberMap[i] and chamberMap[i] <= 10)
		{
			if (chamber != -1 and chamber != chamberMap[i])
			{
				Logging(kHLTLogError,
					"AliHLTMUONMansoTrackerFSMComponent::AddRecHits",
					"Invalid block",
					"Received a data block with data from multiple chambers."
					  " This component cannot handle such a case."
				);
				return;
			}
			else
				chamber = chamberMap[i];
		}
		else
		{
			if (fWarnForUnexpecedBlock)
			{
				Logging(kHLTLogWarning,
					"AliHLTMUONMansoTrackerFSMComponent::AddRecHits",
					"Invalid chamber",
					"Received a data block with data from chamber %d"
					 " which is outside the expected range: [7..10].",
					chamberMap[i]
				);
			}
			return;
		}
	}
	
	// Make sure we got one chamber number.
	if (chamber < 7 or 10 < chamber)
	{
		Logging(kHLTLogError,
			"AliHLTMUONMansoTrackerFSMComponent::AddRecHits",
			"Invalid block",
			"Received a reconstructed hit data block with a null specification."
			 " Cannot know which chamber the data comes from."
		);
		return;
	}
	
	DebugTrace("Added " << count << " reconstructed hits from chamber "
		<< (int)chamber	<< " to the internal arrays."
	);
	
	assert( fRecHitBlockCount[chamber-7] < fRecHitBlockArraySize );
	AliRecHitBlockInfo info(count, recHits);
	fRecHitBlock[chamber-7][fRecHitBlockCount[chamber-7]] = info;
	fRecHitBlockCount[chamber-7]++;
}


void AliHLTMUONMansoTrackerFSMComponent::ResetCanLoadFlags()
{
	/// Resets all the fCanLoad* flags to true. This enables loading of all
	/// those CDB entries in the method ReadConfigFromCDB.
	
	fCanLoadZmiddle = true;
	fCanLoadBL = true;
	for (int i = 0; i < 4; i++)
	{
		fCanLoadA[i] = true;
		fCanLoadB[i] = true;
	}
	for (int i = 0; i < 6; i++)
	{
		fCanLoadZ[i] = true;
	}
}


bool AliHLTMUONMansoTrackerFSMComponent::AtLeastOneCanLoadFlagsIsSet() const
{
	/// Returns true if at least one fCanLoad* flag was true and false otherwise.

	if (fCanLoadZmiddle or fCanLoadBL) return true;
	for (int i = 0; i < 4; i++)
	{
		if (fCanLoadA[i]) return true;
		if (fCanLoadB[i]) return true;
	}
	for (int i = 0; i < 6; i++)
	{
		if (fCanLoadZ[i]) return true;
	}
	return false;
}


void AliHLTMUONMansoTrackerFSMComponent::RequestClusters(
		AliHLTMUONMansoTrackerFSM* tracker,
		AliHLTFloat32_t left, AliHLTFloat32_t right,
		AliHLTFloat32_t bottom, AliHLTFloat32_t top,
		AliHLTMUONChamberName chamber, const void* tag
	)
{
	///
	/// Inherited from AliHLTMUONMansoTrackerFSMCallback.
	/// This is the call back method used by the tracker algorithm to request
	/// clusters on a certain chamber.
	///

	DebugTrace("AliHLTMUONMansoTracker::RequestClusters(chamber = " << chamber << ")");
	void* ctag = const_cast<void*>(tag);
	int chNo = -1;
	AliHLTUInt32_t recHitsCount = 0;
	AliRecHitBlockInfo* recHitsBlock = NULL;
	switch (chamber)
	{
	case kChamber7:
		recHitsCount = fRecHitBlockCount[0];
		recHitsBlock = fRecHitBlock[0];
		chNo = 7;
		break;

	case kChamber8:
		recHitsCount = fRecHitBlockCount[1];
		recHitsBlock = fRecHitBlock[1];
		chNo = 8;
		break;

	case kChamber9:
		recHitsCount = fRecHitBlockCount[2];
		recHitsBlock = fRecHitBlock[2];
		chNo = 9;
		break;

	case kChamber10:
		recHitsCount = fRecHitBlockCount[3];
		recHitsBlock = fRecHitBlock[3];
		chNo = 10;
		break;

	default: return;
	}
	
	DebugTrace("Returning requested hits for chamber " << chNo << ":");
	for (AliHLTUInt32_t i = 0; i < recHitsCount; i++)
	for (AliHLTUInt32_t j = 0; j < recHitsBlock[i].Count(); j++)
	{
		const AliHLTMUONRecHitStruct* hit = &(recHitsBlock[i].Data()[j]);
		if (left < hit->fX and hit->fX < right and bottom < hit->fY and hit->fY < top)
			tracker->ReturnClusters(ctag, hit, 1);
	}
	DebugTrace("Done returning hits from chamber " << chNo << ".");
	tracker->EndOfClusters(ctag);
}


void AliHLTMUONMansoTrackerFSMComponent::EndOfClusterRequests(
		AliHLTMUONMansoTrackerFSM* /*tracker*/
	)
{
	///
	/// Inherited from AliHLTMUONMansoTrackerFSMCallback.
	/// Nothing special to do here.
	///
	
	DebugTrace("End of cluster requests.");
}


void AliHLTMUONMansoTrackerFSMComponent::FoundTrack(AliHLTMUONMansoTrackerFSM* tracker)
{
	///
	/// Inherited from AliHLTMUONMansoTrackerFSMCallback.
	/// This is the call back method used by the tracker algorithm to declare
	/// that a new track has been found.
	///
	
	DebugTrace("AliHLTMUONMansoTrackerFSMComponent::FoundTrack()");
	
	AliHLTMUONMansoTracksBlockWriter* block =
		reinterpret_cast<AliHLTMUONMansoTracksBlockWriter*>(fBlock);
	
	AliHLTMUONMansoTrackStruct newTrack;
	tracker->FillTrackData(newTrack);
	
	// The indicies of the duplicate tracks. If set to block->Nentries() then
	// this indicates the index is not used.
	AliHLTUInt32_t dup1 = block->Nentries();
	AliHLTUInt32_t dup2 = block->Nentries();
	
	// Check if there are any tracks that use the same hits as the one found.
	// If there are, then use the one that has the highest pT.
	// There will be at most 2 duplicate tracks.
	for (AliHLTUInt32_t i = 0; i < block->Nentries(); i++)
	{
		AliHLTMUONMansoTrackStruct& track = (*block)[i];
		bool hasNoDuplicates = true;
		for (AliHLTUInt32_t j = 0; j < 4; j++)
		{
			if (track.fHit[j] == AliHLTMUONConstants::NilRecHitStruct()) continue;
			if (newTrack.fHit[j] == AliHLTMUONConstants::NilRecHitStruct()) continue;
			if (track.fHit[j] == newTrack.fHit[j])
			{
				hasNoDuplicates = false;
				break;
			}
		}
		if (hasNoDuplicates) continue;
		
		if (dup1 == block->Nentries())
		{
			dup1 = i;
		}
		else if (dup2 == block->Nentries())
		{
			dup2 = i;
		}
		else
		{
			HLTError("Found more than 2 tracks with duplicate hits. This is completely unexpected. Something is seriously wrong!");
		}
	}
	
	if (dup1 != block->Nentries() and dup2 != block->Nentries())
	{
		// In this case we found 2 duplicate entries.
		// Figure out which one has the highest pT and keep only that one.
		AliHLTMUONMansoTrackStruct& track1 = (*block)[dup1];
		AliHLTMUONMansoTrackStruct& track2 = (*block)[dup2];
		double newPt = sqrt(newTrack.fPx * newTrack.fPx + newTrack.fPy * newTrack.fPy);
		double dupPt1 = sqrt(track1.fPx * track1.fPx + track1.fPy * track1.fPy);
		double dupPt2 = sqrt(track2.fPx * track2.fPx + track2.fPy * track2.fPy);
		
		if (newPt >= dupPt1 and newPt >= dupPt2)
		{
			// The new track must replace both existing tracks.
			track1 = newTrack;
			track2 = (*block)[block->Nentries()-1];
		}
		else if (dupPt1 >= newPt and dupPt1 >= dupPt2)
		{
			// track1 has the highest pT so ignore the new track and delete track2.
			track2 = (*block)[block->Nentries()-1];
		}
		else
		{
			// In this case track2 must have the highest pT so ignore the new
			// track and delete track1.
			track1 = (*block)[block->Nentries()-1];
		}
		
		// Decrement the number of entries because we deleted a track.
		assert(fTrackCount > 0);
		assert(block->Nentries() > 0);
		block->SetNumberOfEntries(block->Nentries()-1);
		fTrackCount--;
	}
	else if (dup1 != block->Nentries())
	{
		// Only one track with duplicate hits found.
		// See if the new track has higher pT. If it does then replace the
		// exisiting track, otherwise ignore the new track.
		AliHLTMUONMansoTrackStruct& track1 = (*block)[dup1];
		double newPt = sqrt(newTrack.fPx * newTrack.fPx + newTrack.fPy * newTrack.fPy);
		double dupPt1 = sqrt(track1.fPx * track1.fPx + track1.fPy * track1.fPy);
		if (newPt >= dupPt1)
		{
			track1 = newTrack;
		}
	}
	else
	{
		// No track found with duplicate hits so we can add the new track as it is.
		AliHLTMUONMansoTrackStruct* track = block->AddEntry();
		if (track == NULL)
		{
			Logging(kHLTLogError,
				"AliHLTMUONMansoTrackerFSMComponent::FoundTrack",
				"Buffer overflow",
				"We have overflowed the output buffer for Manso track data."
				" The output buffer size is only %d bytes.",
				block->BufferSize()
			);
			return;
		}
	
		fTrackCount++;
		*track = newTrack;
		DebugTrace("\tAdded new track: " << *track);
	}
}


void AliHLTMUONMansoTrackerFSMComponent::NoTrackFound(AliHLTMUONMansoTrackerFSM* /*tracker*/)
{
	///
	/// Inherited from AliHLTMUONMansoTrackerFSMCallback.
	/// Nothing special to do here.
	///
	
	DebugTrace("No track found.");
}

