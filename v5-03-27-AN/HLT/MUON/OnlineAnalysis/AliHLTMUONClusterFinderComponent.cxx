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

// $Id:  $

///
/// @file   AliHLTMUONClusterFinderComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 May 2007
/// @brief  Implementation of the offline algorithm cluster finding processing component.
///
/// The cluster finder component interfaces the offline MUON reconstruction
/// algorithms for cluster finding with the online HLT system.
///

#include "AliHLTMUONClusterFinderComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONDataBlockWriter.h"
#include "AliHLTLogging.h"
#include "AliHLTSystem.h"
#include "AliHLTDefinitions.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONClusterStoreV2.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONVClusterFinder.h"
#include "AliMUONSimpleClusterServer.h"
#include "AliMUONRecoParam.h"
#include "AliMUONReconstructor.h"
#include "AliMUONVCluster.h"
#include "AliMUONRawStreamTrackerHP.h"
#include "AliMpConstants.h"
#include "AliMpCDB.h"
#include "AliMpArea.h"
#include "AliRawReader.h"
#include "AliRawReaderMemory.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "TMap.h"
#include <cstdlib>
#include <cerrno>
#include <cassert>


ClassImp(AliHLTMUONClusterFinderComponent)


AliHLTMUONClusterFinderComponent::AliHLTMUONClusterFinderComponent() :
	AliHLTMUONProcessor(),
	fRawReader(NULL),
	fDigitMaker(NULL),
	fTransformer(NULL),
	fCalibrationData(NULL),
	fDigitCalibrator(NULL),
	fClusterFinder(NULL),
	fClusterServer(NULL),
	fRecoParam(NULL),
	fMakeClusterStore(true),
	fMakeRecHits(false)
{
	/// Default constructor.
}


AliHLTMUONClusterFinderComponent::~AliHLTMUONClusterFinderComponent()
{
	/// Default destructor.
	
	FreeObjects();
}

const char* AliHLTMUONClusterFinderComponent::GetComponentID()
{
	/// Inherited from AliHLTComponent. Returns the component ID.
	
	return AliHLTMUONConstants::ClusterFinderId();
}


void AliHLTMUONClusterFinderComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
	/// Inherited from AliHLTProcessor. Returns the list of expected input data types.
	
	list.clear();
	list.push_back( AliHLTMUONConstants::DDLRawDataType() );
}


AliHLTComponentDataType AliHLTMUONClusterFinderComponent::GetOutputDataType()
{
	/// Inherited from AliHLTComponent. Returns kAliHLTMultipleDataType
	/// refer to GetOutputDataTypes for all returned data types.
	
	return kAliHLTMultipleDataType;
}


int AliHLTMUONClusterFinderComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
	/// Inherited from AliHLTComponent. Returns the output data types.
	
	assert( list.empty() );
	list.push_back( AliHLTMUONConstants::ClusterStoreDataType() );
	list.push_back( AliHLTMUONConstants::RecHitsBlockDataType() );
	return list.size();
}


void AliHLTMUONClusterFinderComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	
	constBase = sizeof(AliMUONVClusterStore) + 1024*1024;
	inputMultiplier = 1;
}


AliHLTComponent* AliHLTMUONClusterFinderComponent::Spawn()
{
	/// Inherited from AliHLTComponent. Creates a new object instance.
	
	return new AliHLTMUONClusterFinderComponent;
}


int AliHLTMUONClusterFinderComponent::DoInit(int argc, const char** argv)
{
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.

	HLTInfo("Initialising dHLT cluster finder component.");

	// Inherit the parents functionality.
	int result = AliHLTMUONProcessor::DoInit(argc, argv);
	if (result != 0) return result;
	
	// Must make sure that all the offline reconstruction objects are released in case
	// they are still allocated for whatever reason.
	FreeObjects();
	
	// Initialise fields with default values then parse the command line.
	fMakeClusterStore = true;
	fMakeRecHits = false;
	bool tryRecover = false;
	
	for (int i = 0; i < argc; i++)
	{
		if (ArgumentAlreadyHandled(i, argv[i])) continue;
		
		if (strcmp( argv[i], "-tryrecover" ) == 0)
		{
			tryRecover = true;
			continue;
		}
		
		if (strcmp( argv[i], "-nostore" ) == 0)
		{
			fMakeClusterStore = false;
			continue;
		}
		
		if (strcmp( argv[i], "-makehits" ) == 0)
		{
			fMakeRecHits = true;
			continue;
		}
		
		HLTError("Unknown option '%s'", argv[i]);
		return -EINVAL;
	
	} // for loop

	try
	{
		fRawReader = new AliRawReaderMemory();
		fDigitMaker = new AliMUONDigitMaker(true, true, true);
		fTransformer = new AliMUONGeometryTransformer();
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the MUON offline reconstruction objects.");
		FreeObjects();
		return -ENOMEM;
	}
	
	if (not DelaySetup())
	{
		result = ReadConfigFromCDB();
		if (result != 0)
		{
			// Error messages already generated in ReadConfigFromCDB.
			FreeObjects(); // Make sure we cleanup to avoid partial initialisation.
			return result;
		}
	}

#ifndef HAVE_NOT_MUON_DIGITMAKER_GETRAWSTREAM	
	AliMUONRawStreamTrackerHP* stream = static_cast<AliMUONRawStreamTrackerHP*>(
			fDigitMaker->GetRawStreamTracker()
		);
	stream->TryRecover(tryRecover);
#endif //HAVE_NOT_MUON_DIGITMAKER_GETRAWSTREAM
	
	return 0;
}


int AliHLTMUONClusterFinderComponent::DoDeinit()
{
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	
	HLTInfo("Deinitialising dHLT cluster finder component.");
	FreeObjects();
	return 0;
}


int AliHLTMUONClusterFinderComponent::Reconfigure(
		const char* cdbEntry, const char* componentId
	)
{
	/// Inherited from AliHLTComponent. This method will reload CDB configuration
	/// entries for this component from the CDB.
	/// \param cdbEntry If set to NULL or starts with "MUON/" then the DDL store
	///      for MUON will be loaded which contains the mapping and also the
	///      calibration data will be updated.
	/// \param componentId  The name of the component in the current chain.
	
	bool startsWithMUON = TString(cdbEntry).Index("MUON/", 5, 0, TString::kExact) == 0;
	
	if (cdbEntry == NULL or startsWithMUON)
	{
		HLTInfo("Reading new configuration entries from CDB for component '%s'.", componentId);
		
		//TODO: add more granularity to the loading.
		return ReadConfigFromCDB();
	}
	
	return 0;
}


int AliHLTMUONClusterFinderComponent::ReadPreprocessorValues(const char* modules)
{
	/// Inherited from AliHLTComponent. 
	/// Updates the configuration of this component if either HLT or MUON have
	/// been specified in the 'modules' list.

	TString mods = modules;
	if (mods.Contains("ALL") or mods.Contains("MUON"))
	{
		return Reconfigure(NULL, GetComponentID());
	}
	if (mods.Contains("MUON"))
	{
		return Reconfigure("MUON/*", GetComponentID());
	}
	return 0;
}


int AliHLTMUONClusterFinderComponent::DoEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr,
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks
	)
{
	/// Inherited from AliHLTProcessor. Processes the new event data by
	/// applying the offline cluster finding algorithms.
	
	// Initialise the mapping and calibration from CDB if we were requested
	// to initialise only when the first event was received.
	if (DelaySetup())
	{
		int result = ReadConfigFromCDB();
		if (result != 0) return result;
		DoneDelayedSetup();
	}
	
	assert(fRawReader != NULL);
	assert(fDigitMaker != NULL);
	assert(fTransformer != NULL);
	assert(fCalibrationData != NULL);
	assert(fDigitCalibrator != NULL);
	assert(fClusterFinder != NULL);
	assert(fClusterServer != NULL);
	
	AliHLTUInt32_t specification = 0x0;  // Accumulated specification to use for output blocks.
	
	HLTDebug("Processing event %llu with %u input data blocks.",
		evtData.fEventID, evtData.fBlockCnt
	);

	AliMUONVDigitStore* digitStore = NULL;
	AliMUONVClusterStore* clusterStore = NULL;
	try
	{
		digitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2R");
		clusterStore = new AliMUONClusterStoreV2();
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the digit or cluster store object.");
		if (digitStore != NULL) delete digitStore;
		if (clusterStore != NULL) delete clusterStore;
		size = 0;
		return -ENOMEM;
	}

	const AliHLTComponentBlockData* block = NULL;
	const AliHLTComponentDataType& rawType = AliHLTMUONConstants::DDLRawDataType();
	for (block = GetFirstInputBlock(rawType); block != NULL; block = GetNextInputBlock())
	{
		HLTDebug("Handling block with fDataType = '%s', fSpecification = 0x%8.8X, fPtr = %p and fSize = %u bytes.",
			DataType2Text(block->fDataType).c_str(), block->fSpecification,
			block->fPtr, block->fSize
		);

		if (not AliHLTMUONUtils::IsTrackerDDL(block->fSpecification))
		{
			HLTError("Received raw data from a DDL that was not a tracker DDL."
				" The data block specification was: 0x%8.8X."
				" Will skip the data block.",
				block->fSpecification
			);
			continue;
		}

		specification |= block->fSpecification;

		fRawReader->SetMemory(reinterpret_cast<UChar_t*>(block->fPtr), ULong_t(block->fSize));
		fRawReader->SetEquipmentID(AliHLTMUONUtils::SpecToEquipId(block->fSpecification));
		fRawReader->Reset();
		fRawReader->NextEvent();
		fDigitMaker->Raw2Digits(fRawReader, digitStore, NULL);
#ifndef HAVE_NOT_MUON_DIGITMAKER_GETRAWSTREAM
		if (fDigitMaker->GetRawStreamTracker()->IsErrorMessage() and DumpDataOnError())
		{
			DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
		}
#endif //HAVE_NOT_MUON_DIGITMAKER_GETRAWSTREAM
	}

	fDigitCalibrator->Calibrate(*digitStore);
	TIter next(digitStore->CreateIterator());
	fClusterServer->UseDigits(next);
	AliMpArea area;
	for (Int_t i = 0; i < AliMpConstants::NofTrackingChambers(); ++i)
	{
		if (fRecoParam->UseChamber(i))
		{
			if ( ( i == 6 or i == 7 ) and fRecoParam->BypassSt4() ) continue;
			if ( ( i == 8 or i == 9 ) and fRecoParam->BypassSt5() ) continue;
			fClusterServer->Clusterize(i, *clusterStore, area);
		}
	}
	
	// Now write the clusters to output blocks.
	AliHLTUInt8_t* buffer = outputPtr;
	AliHLTUInt32_t bufferSize = size;

	if (fMakeClusterStore)
	{
		PushBack(clusterStore, AliHLTMUONConstants::ClusterStoreDataType(), specification);
	
		if (fMakeRecHits)
		{
			try
			{
				bufferSize = sizeof(AliHLTMUONRecHitsBlockWriter::HeaderType)
					+ sizeof(AliHLTMUONRecHitsBlockWriter::ElementType) * clusterStore->GetSize();
				buffer = new AliHLTUInt8_t[size];
			}
			catch (const std::bad_alloc&)
			{
				HLTError("Could not allocate more memory for the reconstructed hit buffer.");
				if (digitStore != NULL) delete digitStore;
				if (clusterStore != NULL) delete clusterStore;
				return -ENOMEM;
			}
		}
	}

	if (fMakeRecHits)
	{
		AliHLTMUONRecHitsBlockWriter outBlock(buffer, bufferSize);
		outBlock.InitCommonHeader();
		AliHLTUInt32_t i = 0;
		TIter next2(clusterStore->CreateIterator());
		AliMUONVCluster* cluster = NULL;
		while ( (cluster = static_cast<AliMUONVCluster*>(next2())) != NULL )
		{
			AliHLTMUONRecHitStruct& hit = outBlock[i++];
			hit.fFlags = AliHLTMUONUtils::PackRecHitFlags(
					AliHLTUInt8_t(cluster->GetChamberId()),
					AliHLTUInt16_t(cluster->GetDetElemId())
				);
			hit.fX = cluster->GetX();
			hit.fY = cluster->GetY();
			hit.fZ = cluster->GetZ();
		}
		outBlock.SetNumberOfEntries(i);

		PushBack(
			buffer, outBlock.BytesUsed(),
			AliHLTMUONConstants::RecHitsBlockDataType(), specification
		);

		if (fMakeClusterStore)
		{
			delete [] buffer;
		}
	}

	delete digitStore;
	delete clusterStore;

	return 0;
}


void AliHLTMUONClusterFinderComponent::FreeObjects()
{
	/// Deletes any allocated objects, if they are allocated, else nothing is
	/// done for objects not yet allocated.
	/// This is used as a helper method to make sure the corresponding pointers
	/// are NULL and we get back to a well defined state.

	if (fRawReader != NULL)
	{
		delete fRawReader;
		fRawReader = NULL;
	}
	if (fDigitMaker != NULL)
	{
		delete fDigitMaker;
		fDigitMaker = NULL;
	}
	if (fTransformer != NULL)
	{
		delete fTransformer;
		fTransformer = NULL;
	}
	if (fCalibrationData != NULL)
	{
		delete fCalibrationData;
		fCalibrationData = NULL;
	}
	if (fDigitCalibrator != NULL)
	{
		delete fDigitCalibrator;
		fDigitCalibrator = NULL;
	}
	if (fClusterServer != NULL)
	{
		delete fClusterServer;
		fClusterServer = NULL;
		fClusterFinder = NULL;  // The cluster server takes ownership.
	}

	// The following is just in case we created the cluster finder, but could
	// not yet create the cluster server.
	if (fClusterFinder != NULL)
	{
		delete fClusterFinder;
		fClusterFinder = NULL;
	}

	if (fRecoParam != NULL)
	{
		delete fRecoParam;
		fRecoParam = NULL;
	}
}


int AliHLTMUONClusterFinderComponent::ReadConfigFromCDB(
		bool loadParams, bool loadMapping, bool loadGeom, bool loadCalib
	)
{
	/// Loads the various configuration, calibration, mapping and geometry
	/// data from CDB.

	if (loadMapping)
	{
		HLTInfo("Loading mapping information from CDB.");
		int result = FetchMappingStores();
		if (result != 0) return result;
	}

	if (loadParams)
	{
		HLTInfo("Loading reconstruction parameters from CDB.");
		AliMUONRecoParam* recoParam = NULL;
		try
		{
			//TODO:
			HLTWarning("Have not yet implemented loading reco params from CDB! Using AliMUONRecoParam::GetLowFluxParam().");
			//result = LoadRecoParamsFromCDB(recoParam);
			//if (result != 0) return result;
			recoParam = AliMUONRecoParam::GetLowFluxParam();
		}
		catch (const std::bad_alloc&)
		{
			HLTError("Could not allocate more memory for the reconstruction parameter object.");
			return -ENOMEM;
		}
		if (fRecoParam != NULL) delete fRecoParam;
		fRecoParam = recoParam;
	}

	if (loadGeom)
	{
		HLTInfo("Loading geometry data from CDB.");
		// Only load geometry if not already loaded.
		if (AliGeomManager::GetGeometry() == NULL)
		{
			AliGeomManager::LoadGeometry();
		}
		assert(fTransformer != NULL);
		if (not fTransformer->LoadGeometryData())
		{
			HLTError("Could not load geometry data for transformation.");
			return -ENOENT;
		}
	}

	if (loadCalib)
	{
		assert(fRecoParam != NULL);

		HLTInfo("Loading calibration information from CDB.");
		AliMUONCalibrationData* calibData = NULL;
		try
		{
			assert(AliCDBManager::Instance() != NULL);
			Int_t runNumber = AliCDBManager::Instance()->GetRun();
			calibData = new AliMUONCalibrationData(runNumber);
		}
		catch (const std::bad_alloc&)
		{
			HLTError("Could not allocate more memory for the calibration data object.");
			return -ENOMEM;
		}
		
		if (not calibData->IsValid())
		{
			HLTError("Could not retrieve calibrations!");
			delete calibData;
			return -ENOENT;
		}

		// Check that we get all the calibrations we'll need.
		if (not calibData->Pedestals() or
		    not calibData->Gains() or
		    not calibData->HV() )
		{
			HLTError("Could not access all required calibration data.");
			delete calibData;
			return -ENOENT;
		}
		
		if (fCalibrationData != NULL) delete fCalibrationData;
		fCalibrationData = calibData;


		AliMUONDigitCalibrator* calibrator = NULL;
		AliMUONVClusterFinder* clusterFinder = NULL;
		AliMUONSimpleClusterServer* clusterServer = NULL;
		try
		{
			calibrator = new AliMUONDigitCalibrator(*fCalibrationData, fRecoParam);
		
			clusterFinder = AliMUONReconstructor::CreateClusterFinder(
					fRecoParam->GetClusteringMode()
				);
	
			clusterServer = new AliMUONSimpleClusterServer(clusterFinder, *fTransformer);
		}
		catch (const std::bad_alloc&)
		{
			HLTError("Could not allocate more memory for a offline reconstruction object.");
			if (calibrator != NULL) delete calibrator;
			if (clusterFinder != NULL) delete clusterFinder;
			if (clusterServer != NULL) delete clusterServer;
			return -ENOMEM;
		}
		
		if (fDigitCalibrator != NULL) delete fDigitCalibrator;
		fDigitCalibrator = calibrator;
		if (fClusterFinder != NULL) delete fClusterFinder;
		fClusterFinder = clusterFinder;
		if (fClusterServer != NULL) delete fClusterServer;
		fClusterServer = clusterServer;
	}

	return 0;
}

