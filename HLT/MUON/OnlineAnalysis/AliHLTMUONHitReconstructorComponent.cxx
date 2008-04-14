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
#include "AliHLTLogging.h"
#include "AliHLTSystem.h"
#include "AliHLTDefinitions.h"
#include <cstdlib>
#include <cerrno>
#include <cassert>

//STEER 
#include "AliCDBManager.h"
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
	AliHLTProcessor(),
	fHitRec(NULL),
	fDDL(0),
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
	
	try
	{
		fHitRec = new AliHLTMUONHitReconstructor();
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the hit reconstructor component.");
		return -ENOMEM;
	}
	
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
				return EINVAL;
			}
			
			char* cpErr = NULL;
			unsigned long num = strtoul( argv[i+1], &cpErr, 0 );
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to DDL Number ", argv[i+1] );
				return EINVAL;
			}
			if (num < 13 or 20 < num)
			{
				HLTError("The DDL number must be in the range [13..20].");
				return EINVAL;
			}
			fDDL = num - 1;
			
			i++;
			continue;
		} // -ddl argument
		
		if (strcmp( argv[i], "-lut" ) == 0)
		{
			if (argc <= i+1)
			{
				HLTError("The lookup table filename was not specified.");
				return EINVAL;
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
				return EINVAL;
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
				return EINVAL;
			}
			
			char* cpErr = NULL;
			run = Int_t( strtoul(argv[i+1], &cpErr, 0) );
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a valid run number."
					" Expected an integer value.", argv[i+1]
				);
				return EINVAL;
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
		return EINVAL;
	
	} // for loop
	
	if (lutFileName == NULL) useCDB = true;
	
	AliHLTMUONHitRecoLutRow* lookupTable;
	
	if (useCDB)
	{
		if(!ReadCDB(lookupTable,cdbPath,run))
		{
			HLTError("Failed to read cdb, cdb cannot be read, DoInit");
			
			if (fHitRec)
			{
				delete fHitRec;
				fHitRec = NULL;
			}
			if(lookupTable)
				delete []lookupTable;
			
			return ENOENT ; /* No such file or directory */
		}
	}
	else
	{
		AliHLTUInt32_t lutLine;
		if(!GetLutLine(lutFileName,fDDL,lutLine))
		{
			HLTError("Failed for lookuptable count the number of lines in lookuptable, DoInit");
			
			if(fHitRec)
				delete fHitRec;
			return EIO;
		}
	
		try
		{
			lookupTable = new AliHLTMUONHitRecoLutRow[lutLine];
		}
		catch(const std::bad_alloc&)
		{
			HLTError("Dynamic memory allocation failed for lookuptable, DoInit");
			
			if(fHitRec)
				delete fHitRec;
			
			return ENOMEM;
		}
	
	
		if(!ReadLookUpTable(lookupTable,lutFileName))
		{
			HLTError("Failed to read lut, lut cannot be read, DoInit");
			
			if(fHitRec)
				delete fHitRec;
			if(lookupTable)
				delete []lookupTable;
			
			return ENOENT ; /* No such file or directory */
		}
	}
	
	if(!fHitRec->SetIdManuChannelToEntry(fIdToEntry))
	{
		HLTError("Failed to set fIdToEntry mapping, DoInit");
		
		if(fHitRec)
		delete fHitRec;
		if(lookupTable)
		delete []lookupTable;
		
		fIdToEntry.clear();
		
		return ENOENT ; /* No such file or directory */
	}
	
	if(!fHitRec->LoadLookUpTable(lookupTable,fDDL))
	{
		HLTError("Cannot Laod hitrec lookuptable , DoInit");
		
		if(fHitRec)
		delete fHitRec;
		if(lookupTable)
		delete []lookupTable;
		
		fIdToEntry.clear();
		
		return ENOENT;
	}
	
	delete [] lookupTable;
	
	return 0;
}


int AliHLTMUONHitReconstructorComponent::DoDeinit()
{
	///
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	///
	
	HLTInfo("Deinitialising dHLT hit reconstruction component.");
	
	if (fHitRec != NULL)
	{
		delete fHitRec;
		fHitRec = NULL;
	}
	
	fIdToEntry.clear();
	
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
	for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
	{
#ifdef __DEBUG
		char id[kAliHLTComponentDataTypefIDsize+1];
		for (int i = 0; i < kAliHLTComponentDataTypefIDsize; i++)
			id[i] = blocks[n].fDataType.fID[i];
		id[kAliHLTComponentDataTypefIDsize] = '\0';
		char origin[kAliHLTComponentDataTypefOriginSize+1];
		for (int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++)
			origin[i] = blocks[n].fDataType.fOrigin[i];
		origin[kAliHLTComponentDataTypefOriginSize] = '\0';
#endif // __DEBUG
		HLTDebug("Handling block: %u, with fDataType.fID = '%s',"
			  " fDataType.fID = '%s', fPtr = %p and fSize = %u bytes.",
			n, static_cast<char*>(id), static_cast<char*>(origin),
			blocks[n].fPtr, blocks[n].fSize
		);

		if (blocks[n].fDataType != AliHLTMUONConstants::DDLRawDataType()
		    or not AliHLTMUONUtils::IsTrackerDDL(blocks[n].fSpecification)
		   )
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
				HLTWarning("Received a data block of a type we cannot handle: '%s' origin: '%s' spec: 0x%X",
					static_cast<char*>(id), static_cast<char*>(origin), blocks[n].fSpecification
				);
			else
				HLTDebug("Received a data block of a type we cannot handle: '%s' origin: '%s' spec: 0x%X",
					static_cast<char*>(id), static_cast<char*>(origin), blocks[n].fSpecification
				);
			
			continue;
		}
		
		bool ddl[22];
		AliHLTMUONUtils::UnpackSpecBits(blocks[n].fSpecification, ddl);
		if (not ddl[fDDL])
		{
			HLTWarning("Received raw data from an unexpected DDL.");
		}
		
		// Create a new output data block and initialise the header.
		AliHLTMUONRecHitsBlockWriter block(outputPtr+totalSize, size-totalSize);
		if (not block.InitCommonHeader())
		{
			HLTError("There is not enough space in the output buffer for the new data block.",
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
			HLTError("Error while processing of hit reconstruction algorithm.");
			size = totalSize; // Must tell the framework how much buffer space was used.
			return EIO;
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


bool AliHLTMUONHitReconstructorComponent::GetLutLine(
		const char* lutFileName, AliHLTInt32_t /*iDDL*/, AliHLTUInt32_t& iLine
	)
{
  // Reads LUT from CDB.
  // TODO: combine this with ReadLookUpTable().
  
  ifstream fin(lutFileName);

  if(!fin){
    HLTError("Failed to open file '%s' ",lutFileName);
    return false;
  }

  string s;
  iLine = 0;
  while(getline(fin,s)){
    iLine++;
  }

  fin.close();

  return true;
}


bool AliHLTMUONHitReconstructorComponent::ReadLookUpTable(
		AliHLTMUONHitRecoLutRow* lookupTable, const char* lutFileName
	)
{
	///
	/// Read in the lookup table from a text file.
	///
  
  if (fDDL < AliHLTMUONHitReconstructor::GetkDDLOffSet() ||
      fDDL >= AliHLTMUONHitReconstructor::GetkDDLOffSet() + AliHLTMUONHitReconstructor::GetkNofDDL())
  {
    HLTError("DDL number is out of range");
    return false;
  }
  
  AliHLTUInt32_t lutLine;
  if(!GetLutLine(lutFileName,fDDL,lutLine)){
    HLTError("Failed for lookuptable count the number of lines in lookuptable, DoInit");
    
    return false;
  }
  
  int idManuChannel;
  fIdToEntry.clear();
  
  FILE *fin = fopen(lutFileName,"r");
  if(fin == NULL){
    printf("Failed to open file %s\n",lutFileName);
    return false;
  }else{
    for(AliHLTUInt32_t i=0;i<lutLine;i++){
      fscanf(fin,"%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%d\t%d\n",
	     &idManuChannel,&lookupTable[i].fDetElemId,&lookupTable[i].fIX,
	     &lookupTable[i].fIY,&lookupTable[i].fRealX,
	     &lookupTable[i].fRealY,&lookupTable[i].fRealZ,
	     &lookupTable[i].fHalfPadSize,&lookupTable[i].fPlane,
	     &lookupTable[i].fPed,&lookupTable[i].fSigma,&lookupTable[i].fA0,
	     &lookupTable[i].fA1,&lookupTable[i].fThres,&lookupTable[i].fSat);
      
      fIdToEntry[idManuChannel] = i+1;

      
    }
  }

  fclose(fin);

  return true;
}


bool AliHLTMUONHitReconstructorComponent::ReadCDB(
		AliHLTMUONHitRecoLutRow*& lookupTable,
		const char* cdbPath, Int_t run
	)
{
  // Reads LUT from CDB.
  // TODO: merge this with CreateHitRecoLookupTables.C, make this static and use in the macro for example.

  vector<AliHLTMUONHitRecoLutRow> lutList;
  lutList.clear();
  AliHLTMUONHitRecoLutRow lut;
  int iEntry = 0;

  Bool_t warn = kTRUE;

  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if (cdbPath != NULL) cdbManager->SetDefaultStorage(cdbPath);
  if (run != -1) cdbManager->SetRun(run);

  if (! AliMpCDB::LoadDDLStore(warn)){
    HLTError("Failed to Load DDLStore specified for CDBPath '%s', and Run : '%d'",cdbPath,run);
    return false;
  }

  AliMpSegmentation *mpSegFactory = AliMpSegmentation::Instance();
  AliGeomManager::LoadGeometry();
  AliMUONGeometryTransformer* chamberGeometryTransformer = new AliMUONGeometryTransformer();
  if(! chamberGeometryTransformer->LoadGeometryData()){
    HLTError("Failed to Load Geomerty Data ");
    return false;
  }
  
  AliMUONCalibrationData calibData(run);

  int totMaxIX = -1;
  int totMaxIY = -1;
  Int_t chamberId;

  for(Int_t iCh = 6; iCh < 10; iCh++){ // max 4

    chamberId = iCh ;

    AliMpDEIterator it;
    for ( it.First(chamberId); ! it.IsDone(); it.Next() ) {
    
      Int_t detElemId = it.CurrentDEId();
      int iDDL = AliMpDDLStore::Instance()->GetDetElement(detElemId)->GetDdlId() + 1;
      if(iDDL == fDDL){

	for(Int_t iCath = 0 ; iCath <= 1 ; iCath++){
	
	  AliMp::CathodType cath;
	  
	  if(iCath == 0)
	    cath = AliMp::kCath0 ;
	  else
	  cath = AliMp::kCath1 ;
	  
	  const AliMpVSegmentation* seg = mpSegFactory->GetMpSegmentation(detElemId, cath);
	  AliMp::PlaneType plane = seg->PlaneType(); 
	  Int_t maxIX = seg->MaxPadIndexX();  
	  Int_t maxIY = seg->MaxPadIndexY(); 
	  if(maxIX > totMaxIX)
	    totMaxIX = maxIX;
	  if(maxIY > totMaxIY)
	    totMaxIY = maxIY;
	  
	  Int_t idManuChannel, manuId, channelId, buspatchId;
	  float padSizeX, padSizeY;
	  float halfPadSize ;
	  Double_t realX, realY, realZ;
	  Double_t localX, localY, localZ;
	  Float_t calibA0Coeff,calibA1Coeff,pedestal,sigma;
	  Int_t thresold,saturation;
	  
	  // 	cout<<"Running for detElemId :"<<detElemId<<", and plane : "<<plane<<endl;
	  //Pad Info of a segment to print in lookuptable
	  for(Int_t iX = 0; iX<= maxIX ; iX++){
	    for(Int_t iY = 0; iY<= maxIY ; iY++){
	      if(seg->HasPad(AliMpIntPair(iX,iY))){
		AliMpPad pad = seg->PadByIndices(AliMpIntPair(iX,iY),kFALSE);
		
		// Getting Manu id
		manuId = pad.GetLocation().GetFirst();
		manuId &= 0x7FF; // 11 bits 

		buspatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
		
		// Getting channel id
		channelId =  pad.GetLocation().GetSecond();
		channelId &= 0x3F; // 6 bits
		
		idManuChannel &= 0x0;
		idManuChannel = (idManuChannel|buspatchId)<<11;  
		idManuChannel = (idManuChannel|manuId)<<6 ;
		idManuChannel |= channelId ;
		
		localX = pad.Position().X();
		localY = pad.Position().Y();
		localZ = 0.0;
		
		chamberGeometryTransformer->Local2Global(detElemId,localX,localY,localZ,
 						       realX,realY,realZ);
		
		padSizeX = pad.Dimensions().X();
		padSizeY = pad.Dimensions().Y();
		
		calibA0Coeff = (calibData.Gains(detElemId,manuId))->ValueAsFloat(channelId,0) ;
		calibA1Coeff = (calibData.Gains(detElemId,manuId))->ValueAsFloat(channelId,1) ;
		thresold = (calibData.Gains(detElemId,manuId))->ValueAsInt(channelId,2) ;
		saturation = (calibData.Gains(detElemId,manuId))->ValueAsInt(channelId,4) ;
		
		pedestal = (calibData.Pedestals(detElemId,manuId))->ValueAsFloat(channelId,0);
		sigma = (calibData.Pedestals(detElemId,manuId))->ValueAsFloat(channelId,1);
		
		if(plane==0)
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
	      }// HasPad Condn
	    }// iY loop
	  }// iX loop
	
	}// iPlane
      }// iDDL
    } // detElemId loop

  }// ichamber loop

  AliHLTMUONHitRecoLutRow *temp;

  try{
    temp = new AliHLTMUONHitRecoLutRow[iEntry];
  }
  catch(const std::bad_alloc&){
    HLTError("Dynamic memory allocation failed for temp");

    return false;
  }
  
  for(int iterm = 0; iterm < iEntry ;iterm++)
    temp[iterm] = lutList.at(iterm);
  
  lookupTable = temp;

  return true;
}
