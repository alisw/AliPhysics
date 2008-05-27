/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
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

/* $Id$ */

///
///
///  The HitRec Component is designed to deal the rawdata inputfiles to findout the 
///  the reconstructed hits. The output is send to the output block for further 
///  processing.
///
///  Author : Indranil Das ( indra.das@saha.ac.in || indra.ehep@gmail.com )
///

#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONHitReconstructorComponent.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONDataBlockWriter.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTLogging.h"
#include "AliHLTSystem.h"
#include "AliHLTDefinitions.h"
#include <cstdlib>
#include <cerrno>
#include <cassert>
#include <fstream>

//STEER 
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"

//MUON
#include "AliMUONGeometryTransformer.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVCalibParam.h"

//MUON/mapping
#include "AliMpCDB.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"

ClassImp(AliHLTMUONHitReconstructorComponent)


AliHLTMUONHitReconstructorComponent::AliHLTMUONHitReconstructorComponent() :
	AliHLTMUONProcessor(),
	fHitRec(NULL),
	fDDL(-1),
	fLutSize(0),
	fLut(NULL),
	fIdToEntry(),
	fWarnForUnexpecedBlock(false)
{
	///
	/// Default constructor.
	///
}


AliHLTMUONHitReconstructorComponent::~AliHLTMUONHitReconstructorComponent()
{
	///
	/// Default destructor.
	///
	
	if (fHitRec != NULL)
	{
		delete fHitRec;
	}
	if (fLut != NULL)
	{
		delete [] fLut;
	}
}

const char* AliHLTMUONHitReconstructorComponent::GetComponentID()
{
	///
	/// Inherited from AliHLTComponent. Returns the component ID.
	///
	
	return AliHLTMUONConstants::HitReconstructorId();
}


void AliHLTMUONHitReconstructorComponent::GetInputDataTypes( std::vector<AliHLTComponentDataType>& list)
{
	///
	/// Inherited from AliHLTProcessor. Returns the list of expected input data types.
	///
	
	list.clear();
	list.push_back( AliHLTMUONConstants::DDLRawDataType() );
}


AliHLTComponentDataType AliHLTMUONHitReconstructorComponent::GetOutputDataType()
{
	///
	/// Inherited from AliHLTComponent. Returns the output data type.
	///
	
	return AliHLTMUONConstants::RecHitsBlockDataType();
}


void AliHLTMUONHitReconstructorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
	///
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	///
	
	constBase = sizeof(AliHLTMUONRecHitsBlockWriter::HeaderType);
	inputMultiplier = 1;
}


AliHLTComponent* AliHLTMUONHitReconstructorComponent::Spawn()
{
	///
	/// Inherited from AliHLTComponent. Creates a new object instance.
	///
	
	return new AliHLTMUONHitReconstructorComponent;
}


int AliHLTMUONHitReconstructorComponent::DoInit(int argc, const char** argv)
{
	///
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	///

	HLTInfo("Initialising dHLT hit reconstruction component.");
	
	// Must make sure that fHitRec and fLut is deleted if it is still
	// allocated for whatever reason.
	FreeMemory();
	
	try
	{
		fHitRec = new AliHLTMUONHitReconstructor();
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the hit reconstructor component.");
		return -ENOMEM;
	}
	
	// Initialise fields with default values then parse the command line.
	fDDL = -1;
	fIdToEntry.clear();
	fWarnForUnexpecedBlock = false;
	
	const char* lutFileName = NULL;
	const char* cdbPath = NULL;
	Int_t run = -1;
	bool useCDB = false;
	
	for (int i = 0; i < argc; i++)
	{
		HLTDebug("argv[%d] == %s", i, argv[i]);
		
		if (strcmp( argv[i], "-ddl" ) == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The DDL number was not specified. Must be in the range [13..20].");
				FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			unsigned long num = strtoul( argv[i+1], &cpErr, 0 );
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to DDL a number.", argv[i+1] );
				FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
				return -EINVAL;
			}
			if (num < 13 or 20 < num)
			{
				HLTError("The DDL number must be in the range [13..20].");
				FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
				return -EINVAL;
			}
			fDDL = num - 1;  // convert to range [12..19]
			
			i++;
			continue;
		} // -ddl argument
		
		if (strcmp( argv[i], "-lut" ) == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The lookup table filename was not specified.");
				FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
				return -EINVAL;
			}
			lutFileName = argv[i+1];
			i++;
			continue;
		} // -lut argument
		
		if (strcmp( argv[i], "-cdb" ) == 0)
		{
			useCDB = true;
			continue;
		} // -cdb argument
		
		if (strcmp( argv[i], "-cdbpath" ) == 0)
		{
			if ( argc <= i+1 )
			{
				HLTError("The CDB path was not specified." );
				FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
				return -EINVAL;
			}
			cdbPath = argv[i+1];
			useCDB = true;
			i++;
			continue;
		} // -cdb argument
	
		if (strcmp( argv[i], "-run" ) == 0)
		{
			if ( argc <= i+1 )
			{
				HLTError("The RUN number was not specified." );
				FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			run = Int_t( strtoul(argv[i+1], &cpErr, 0) );
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a valid run number."
					" Expected an integer value.", argv[i+1]
				);
				FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
				return -EINVAL;
			}
			
			i++;
			continue;
		} // -run argument
		
		if (strcmp( argv[i], "-warn_on_unexpected_block" ) == 0)
		{
			fWarnForUnexpecedBlock = true;
			continue;
		}
	
		HLTError("Unknown option '%s'", argv[i]);
		FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
		return -EINVAL;
	
	} // for loop
	
	if (lutFileName == NULL) useCDB = true;
	
	if (fDDL == -1)
	{
		HLTWarning("DDL number not specified. Cannot check if incomming data is valid.");
	}
	
	int result = 0;
	if (useCDB)
	{
		HLTInfo("Loading lookup table information from CDB for DDL %d.", fDDL+1);
		if (fDDL == -1)
			HLTWarning("DDL number not specified. The lookup table loaded from CDB will be empty!");
		result = ReadCDB(cdbPath, run);
	}
	else
	{
		HLTInfo("Loading lookup table information from file %s.", lutFileName);
		result = ReadLookUpTable(lutFileName);
	}
	if (result != 0)
	{
		// Error messages already generated in ReadCDB or ReadLookUpTable.
		FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
		return result;
	}
	
	fHitRec->SetLookUpTable(fLut, &fIdToEntry);
	
	return 0;
}


int AliHLTMUONHitReconstructorComponent::DoDeinit()
{
	///
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	///
	
	HLTInfo("Deinitialising dHLT hit reconstruction component.");
	FreeMemory();
	return 0;
}


int AliHLTMUONHitReconstructorComponent::DoEvent(
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
	
	// Process an event
	unsigned long totalSize = 0; // Amount of memory currently consumed in bytes.

	HLTDebug("Processing event %llu with %u input data blocks.",
		evtData.fEventID, evtData.fBlockCnt
	);

	// Loop over all input blocks in the event
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
		HLTDebug("Handling block: %u, with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
			n, DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fPtr, blocks[n].fSize
		);

		if (blocks[n].fDataType != AliHLTMUONConstants::DDLRawDataType()
		    or not AliHLTMUONUtils::IsTrackerDDL(blocks[n].fSpecification)
		   )
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
			
			continue;
		}
		
		if (fDDL != -1)
		{
			if (AliHLTMUONUtils::SpecToDDLNumber(blocks[n].fSpecification) != fDDL)
			{
				HLTWarning("Received raw data from an unexpected DDL.");
			}
		}
		
		// Create a new output data block and initialise the header.
		AliHLTMUONRecHitsBlockWriter block(outputPtr+totalSize, size-totalSize);
		if (not block.InitCommonHeader())
		{
			HLTError("There is not enough space in the output buffer for the new data block."
				 " We require at least %u bytes, but have %u bytes left.",
				sizeof(AliHLTMUONRecHitsBlockWriter::HeaderType),
				block.BufferSize()
			);
			break;
		}
		
		AliHLTUInt32_t totalDDLSize = blocks[n].fSize / sizeof(AliHLTUInt32_t);
		AliHLTUInt32_t ddlRawDataSize = totalDDLSize - fHitRec->GetkDDLHeaderSize();
		AliHLTUInt32_t* buffer = reinterpret_cast<AliHLTUInt32_t*>(blocks[n].fPtr)
			+ fHitRec->GetkDDLHeaderSize();
		AliHLTUInt32_t nofHit = block.MaxNumberOfEntries();

#ifdef DEBUG
		HLTDebug("=========== Dumping DDL payload buffer ==========");
		for (AliHLTUInt32_t j = 0; j < totalDDLSize; j++)
			HLTDebug("buffer[%d] : %x",j,buffer[j]);
		HLTDebug("================== End of dump =================");
#endif // DEBUG

		if (not fHitRec->Run(buffer, ddlRawDataSize, block.GetArray(), nofHit))
		{
			HLTError("Error while processing the hit reconstruction algorithm.");
			size = totalSize; // Must tell the framework how much buffer space was used.
			return -EIO;
		}
		
		// nofHit should now contain the number of reconstructed hits actually found
		// and filled into the output data block, so we can set this number.
		assert( nofHit <= block.MaxNumberOfEntries() );
		block.SetNumberOfEntries(nofHit);
		
		HLTDebug("Number of reconstructed hits found is %d", nofHit);

		// Fill a block data structure for our output block.
		AliHLTComponentBlockData bd;
		FillBlockData(bd);
		bd.fPtr = outputPtr;
		// This block's start (offset) is after all other blocks written so far.
		bd.fOffset = totalSize;
		bd.fSize = block.BytesUsed();
		bd.fDataType = AliHLTMUONConstants::RecHitsBlockDataType();
		bd.fSpecification = blocks[n].fSpecification;
		outputBlocks.push_back(bd);

		// Increase the total amount of data written so far to our output memory
		totalSize += block.BytesUsed();
	}
	// Finally we set the total size of output memory we consumed.
	size = totalSize;

	return 0;
}


void AliHLTMUONHitReconstructorComponent::FreeMemory()
{
	/// Deletes any allocated objects if they are allocated else nothing is
	/// done for objects not yet allocated.
	/// This is used as a helper method to make sure the corresponding pointers
	/// are NULL and we get back to a well defined state.

	if (fHitRec != NULL)
	{
		delete fHitRec;
		fHitRec = NULL;
	}
	if (fLut != NULL)
	{
		delete [] fLut;
		fLut = NULL;
		fLutSize = 0;
	}
	
	fIdToEntry.clear();
}


int AliHLTMUONHitReconstructorComponent::ReadLookUpTable(const char* lutFileName)
{
	/// Read in the lookup table from a text file.
	/// Note that this method could leave fLut allocated which is cleaned up
	/// by DoInit with a call to FreeMemory().
	
	assert( fLut == NULL );
	assert( fLutSize == 0 );
	assert( fIdToEntry.empty() );
	
	std::ifstream file(lutFileName);
	if (not file.good())
	{
		HLTError("Could not open the LUT file %s", lutFileName);
		return -ENOENT;
	}
	
	// First count the number of lines of text in the LUT file before decoding.
	// This is not the most optimal. It would be better to read and decode at the
	// same time but we are not allowed to use STL and ROOT containers are too
	// heavy for this task. At least this is done only at the start of run.
	std::string str;
	AliHLTUInt32_t lineCount = 0;
	while (std::getline(file, str)) lineCount++;
	if (not file.eof())
	{
		HLTError("There was a problem reading the LUT file %s", lutFileName);
		return -EIO;
	}
	if (lineCount == 0)
	{
		HLTWarning("The LUT file %s was empty.", lutFileName);
	}
	
	// Add one extra LUT line for the first element which is used as a sentinel value.
	lineCount++;
	
	try
	{
		fLut = new AliHLTMUONHitRecoLutRow[lineCount];
		fLutSize = lineCount;
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the lookuptable.");
		return -ENOMEM;
	}
	
	// Initialise the sentinel value.
	fLut[0].fDetElemId = 0;
	fLut[0].fIX = 0;
	fLut[0].fIY = 0;
	fLut[0].fRealX = 0.0;
	fLut[0].fRealY = 0.0;
	fLut[0].fRealZ = 0.0;
	fLut[0].fHalfPadSize = 0.0;
	fLut[0].fPlane = -1;
	fLut[0].fPed = -1;
	fLut[0].fSigma = -1;
	fLut[0].fA0 = -1;
	fLut[0].fA1 = -1;
	fLut[0].fThres = -1;
	fLut[0].fSat = -1;
	
	// Clear the eof flag and start reading from the beginning of the file again.
	file.clear();
	file.seekg(0, std::ios::beg);
	if (not file.good())
	{
		HLTError("There was a problem seeking in the LUT file %s", lutFileName);
		return -EIO;
	}
	
	AliHLTInt32_t idManuChannel;
	for (AliHLTUInt32_t i = 1; i < fLutSize; i++)
	{
		if (std::getline(file, str).fail())
		{
			HLTError("There was a problem reading line %d of LUT file %s", i, lutFileName);
			return -EIO;
		}
		
		int result = sscanf(
			str.c_str(), "%d\t%d\t%d\t%d\t%e\t%e\t%e\t%e\t%d\t%e\t%e\t%e\t%e\t%d\t%d",
			&idManuChannel, &fLut[i].fDetElemId, &fLut[i].fIX,
			&fLut[i].fIY, &fLut[i].fRealX,
			&fLut[i].fRealY, &fLut[i].fRealZ,
			&fLut[i].fHalfPadSize, &fLut[i].fPlane,
			&fLut[i].fPed, &fLut[i].fSigma, &fLut[i].fA0,
			&fLut[i].fA1, &fLut[i].fThres, &fLut[i].fSat
		);
		
		if (result != 15)
		{
			HLTError("Line %d in LUT file %s does not contain 15 elements.", i, lutFileName);
			return -EIO;
		}
		
		fIdToEntry[idManuChannel] = i;
	}
	
	return 0;
}


int AliHLTMUONHitReconstructorComponent::ReadCDB(const char* cdbPath, Int_t run)
{
	/// Reads LUT from CDB.

	assert( fLut == NULL );
	assert( fLutSize == 0 );
	assert( fIdToEntry.empty() );
	
	std::vector<AliHLTMUONHitRecoLutRow> lutList;
	AliHLTMUONHitRecoLutRow lut;
	AliHLTUInt32_t iEntry = 0;
	
	int result = FetchMappingStores(cdbPath, run);
	// Error message already generated in FetchMappingStores.
	if (result != 0) return result;
	AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
	
	AliMpSegmentation* mpSegFactory = AliMpSegmentation::Instance();
	if (mpSegFactory == NULL)
	{
		HLTError("Could not find segmentation mapping (AliMpSegmentation) instance.");
		return -EIO;
	}
	
	AliGeomManager::LoadGeometry();
	AliMUONGeometryTransformer chamberGeometryTransformer;
	if (not chamberGeometryTransformer.LoadGeometryData())
	{
		HLTError("Failed to load geomerty data.");
		return -ENOENT;
	}
	
	AliMUONCalibrationData calibData(run);
	
	Int_t chamberId;
	
	for(Int_t iCh = 6; iCh < 10; iCh++)
	{
		chamberId = iCh;
		
		AliMpDEIterator it;
		for ( it.First(chamberId); ! it.IsDone(); it.Next() )
		{
			Int_t detElemId = it.CurrentDEId();
			int iDDL = ddlStore->GetDetElement(detElemId)->GetDdlId();
			if (iDDL != fDDL) continue;
		
			for (Int_t iCath = 0 ; iCath <= 1 ; iCath++)
			{
				AliMp::CathodType cath;
				
				if(iCath == 0)
					cath = AliMp::kCath0;
				else
					cath = AliMp::kCath1;
				
				const AliMpVSegmentation* seg = mpSegFactory->GetMpSegmentation(detElemId, cath);
				AliMp::PlaneType plane = seg->PlaneType();
				Int_t maxIX = seg->MaxPadIndexX();
				Int_t maxIY = seg->MaxPadIndexY();
				
				Int_t idManuChannel, manuId, channelId, buspatchId;
				AliHLTFloat32_t padSizeX, padSizeY;
				AliHLTFloat32_t halfPadSize;
				Double_t realX, realY, realZ;
				Double_t localX, localY, localZ;
				Float_t calibA0Coeff,calibA1Coeff,pedestal,sigma;
				Int_t thresold,saturation;
				
				// Pad Info of a slat to print in lookuptable
				for (Int_t iX = 0; iX<= maxIX ; iX++)
				for (Int_t iY = 0; iY<= maxIY ; iY++)
				{
					if (not seg->HasPad(AliMpIntPair(iX,iY))) continue;

					AliMpPad pad = seg->PadByIndices(AliMpIntPair(iX,iY), kFALSE);
					
					// Getting Manu id
					manuId = pad.GetLocation().GetFirst();
					manuId &= 0x7FF; // 11 bits 
			
					buspatchId = ddlStore->GetBusPatchId(detElemId,manuId);
					
					// Getting channel id
					channelId =  pad.GetLocation().GetSecond();
					channelId &= 0x3F; // 6 bits
					
					idManuChannel = buspatchId << 11;
					idManuChannel = (idManuChannel | manuId) << 6;
					idManuChannel |= channelId;
					
					localX = pad.Position().X();
					localY = pad.Position().Y();
					localZ = 0.0;
					
					chamberGeometryTransformer.Local2Global(
						detElemId,localX,localY,localZ,
						realX,realY,realZ
					);
					
					padSizeX = pad.Dimensions().X();
					padSizeY = pad.Dimensions().Y();
					
					calibA0Coeff = (calibData.Gains(detElemId, manuId))->ValueAsFloat(channelId, 0);
					calibA1Coeff = (calibData.Gains(detElemId, manuId))->ValueAsFloat(channelId, 1);
					thresold = (calibData.Gains(detElemId, manuId))->ValueAsInt(channelId, 2);
					saturation = (calibData.Gains(detElemId, manuId))->ValueAsInt(channelId, 4);
					
					pedestal = (calibData.Pedestals(detElemId, manuId))->ValueAsFloat(channelId, 0);
					sigma = (calibData.Pedestals(detElemId, manuId))->ValueAsFloat(channelId, 1);
					
					if (plane == 0)
						halfPadSize = padSizeX;
					else
						halfPadSize = padSizeY;
					
					fIdToEntry[idManuChannel] = iEntry+1;
			
					lut.fDetElemId = detElemId;
					lut.fIX = iX;
					lut.fIY = iY;
					lut.fRealX = realX;
					lut.fRealY = realY;
					lut.fRealZ = realZ;
					lut.fHalfPadSize = halfPadSize;
					lut.fPlane = plane;
					lut.fPed = pedestal;
					lut.fSigma = sigma;
					lut.fA0 = calibA0Coeff;
					lut.fA1 = calibA1Coeff;
					lut.fThres = thresold;
					lut.fSat = saturation;
					
					lutList.push_back(lut);
					iEntry++;
				} // iX, iY loop
			} // iCath loop
		} // detElemId loop
	} // ichamber loop

	try
	{
		// Use iEntry+1 since we add one extra LUT line for the first element
		// which is used as a sentinel value.
		fLut = new AliHLTMUONHitRecoLutRow[iEntry+1];
		fLutSize = iEntry+1;
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the lookuptable.");
		return -ENOMEM;
	}
	
	// Initialise the sentinel value.
	fLut[0].fDetElemId = 0;
	fLut[0].fIX = 0;
	fLut[0].fIY = 0;
	fLut[0].fRealX = 0.0;
	fLut[0].fRealY = 0.0;
	fLut[0].fRealZ = 0.0;
	fLut[0].fHalfPadSize = 0.0;
	fLut[0].fPlane = -1;
	fLut[0].fPed = -1;
	fLut[0].fSigma = -1;
	fLut[0].fA0 = -1;
	fLut[0].fA1 = -1;
	fLut[0].fThres = -1;
	fLut[0].fSat = -1;
	
	for (AliHLTUInt32_t i = 0; i < iEntry; i++)
		fLut[i+1] = lutList[i];
	
	return 0;
}


bool AliHLTMUONHitReconstructorComponent::GenerateLookupTable(
		AliHLTInt32_t ddl, const char* filename,
		const char* cdbPath, Int_t run
	)
{
	/// Generates a ASCII text file containing the lookup table (LUT) from
	/// the CDB, which can be used for the hit reconstructor component later.
	/// @param ddl  Must be the DDL for which to generate the DDL,
	///             in the range [13..20].
	/// @param filename  The name of the LUT file to generate.
	/// @param cdbPath  The CDB path to use.
	/// @param run  The run number to use for the CDB.
	/// @return  True if the generation of the LUT file succeeded.
	
	AliHLTMUONHitReconstructorComponent comp;
	
	if (ddl < 12 or 19 < ddl)
	{
		std::cerr << "ERROR: the DDL number must be in the range [12..19]." << std::endl;
		return false;
	}
	
	comp.fDDL = ddl;
	if (comp.ReadCDB(cdbPath, run) != 0) return false;
	
	char str[1024*4];
	std::fstream file(filename, std::ios::out);
	if (not file)
	{
		std::cerr << "ERROR: could not open file: " << filename << std::endl;
		return false;
	}
	
	assert( comp.fLut != NULL );
	
	for (IdManuChannelToEntry::iterator id = comp.fIdToEntry.begin();
	     id != comp.fIdToEntry.end();
	     id++
	    )
	{
		AliHLTInt32_t idManuChannel = id->first;
		AliHLTInt32_t row = id->second;
		
		assert( row < comp.fLutSize );
		
		sprintf(str, "%d\t%d\t%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\t%d\t%.15e\t%.15e\t%.15e\t%.15e\t%d\t%d",
			idManuChannel, comp.fLut[row].fDetElemId, comp.fLut[row].fIX,
			comp.fLut[row].fIY, comp.fLut[row].fRealX,
			comp.fLut[row].fRealY, comp.fLut[row].fRealZ,
			comp.fLut[row].fHalfPadSize, comp.fLut[row].fPlane,
			comp.fLut[row].fPed, comp.fLut[row].fSigma, comp.fLut[row].fA0,
			comp.fLut[row].fA1, comp.fLut[row].fThres, comp.fLut[row].fSat
		);
		
		file << str << endl;
		if (file.fail())
		{
			std::cerr << "ERROR: There was an I/O error when writing to the file: "
				<< filename << std::endl;
			return false;
		}
	}
	
	return true;
}
