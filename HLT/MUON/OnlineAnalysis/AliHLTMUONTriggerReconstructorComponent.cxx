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

/** @file   AliHLTMUONTriggerReconstructorComponent.cxx
    @author Indranil Das
    @date   
    @brief  Implementation of the trigger DDL reconstructor component. */

#include "AliHLTMUONTriggerReconstructorComponent.h"
#include "AliHLTMUONTriggerReconstructor.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONDataBlockWriter.h"
#include <cstdlib>
#include <cerrno>
#include <cassert>

namespace
{
	// This is a global object used for automatic component registration,
	// do not use this for calculation.
	AliHLTMUONTriggerReconstructorComponent gAliHLTMUONTriggerReconstructorComponent;

} // end of namespace


ClassImp(AliHLTMUONTriggerReconstructorComponent)
    
    
AliHLTMUONTriggerReconstructorComponent::AliHLTMUONTriggerReconstructorComponent() :
	fTrigRec(NULL),
	fDDLDir(""),
	fDDL(0),
	fWarnForUnexpecedBlock(false)
{
}


AliHLTMUONTriggerReconstructorComponent::~AliHLTMUONTriggerReconstructorComponent()
{
}


const char* AliHLTMUONTriggerReconstructorComponent::GetComponentID()
{
	return AliHLTMUONConstants::TriggerReconstructorId();
}


void AliHLTMUONTriggerReconstructorComponent::GetInputDataTypes( std::vector<AliHLTComponentDataType>& list)
{
	list.clear();
	list.push_back( AliHLTMUONConstants::TriggerDDLRawDataType() );
}


AliHLTComponentDataType AliHLTMUONTriggerReconstructorComponent::GetOutputDataType()
{
	return AliHLTMUONConstants::TriggerRecordsBlockDataType();
}


void AliHLTMUONTriggerReconstructorComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	constBase = sizeof(AliHLTMUONTriggerRecordsBlockWriter::HeaderType);
	inputMultiplier = 1;
}


AliHLTComponent* AliHLTMUONTriggerReconstructorComponent::Spawn()
{
	return new AliHLTMUONTriggerReconstructorComponent;
}


int AliHLTMUONTriggerReconstructorComponent::DoInit(int argc, const char** argv)
{
  // perform initialization. We check whether our relative output size is
  // specified in the arguments.
  
  HLTInfo("Initialising DHLT Trigger Record Component");

  fWarnForUnexpecedBlock = false;
  fTrigRec = new AliHLTMUONTriggerReconstructor();
      
  // this is just to get rid of the warning "unused parameter"
  if (argc==0 && argv==NULL) {
    HLTError("Arguments missing, no arguments" );
  }

  char lutFileName[500],reglocFileName[500];

  int i = 0;
  char* cpErr;
  while ( i < argc )
    {
      HLTDebug("argv[%d] == %s", i, argv[i] );
      
      if ( !strcmp( argv[i], "lut" ) ) {
	if ( argc <= i+1 ) {
	  HLTError("LookupTable filename not specified" );
	  return EINVAL; /* Invalid argument */ 
	}
	
	sprintf(lutFileName,"%s",argv[i+1]);
	
	i += 2;
	continue;
      }// lut argument
      
      if ( !strcmp( argv[i], "ddl" ) ) {
	if ( argc <= i+1 ) {
	  HLTError("DDL number not specified" );
	  return EINVAL;  /* Invalid argument */
	}

	fDDL = strtoul( argv[i+1], &cpErr, 0 );
	if ( *cpErr )
	  {
	    HLTError("Cannot convert '%s' to DDL Number ", argv[i+1] );
	    return EINVAL;
	  }
	//fDDL = atoi(argv[i+1]);
	
	i += 2;
	continue;
      }// ddl argument
	  
      if ( !strcmp( argv[i], "rawdir" ) ) {
	if ( argc <= i+1 ) {
	  HLTError("DDL directory not specified" );
	  return EINVAL;  /* Invalid argument */
	}

	fDDLDir = argv[i+1] ;
	i += 2;
	continue;
      }// ddl directory argument

      if ( !strcmp( argv[i], "reglocmap" ) ) {
	if ( argc <= i+1 ) {
	  HLTError("Regional to Local Card mapping  filename not specified" );
	  return EINVAL; /* Invalid argument */
	}

	sprintf(reglocFileName,"%s",argv[i+1]);

	i += 2;
	continue;
      }// regtolocalmap argument
	  
      if ( !strcmp( argv[i], "-warn_on_unexpected_block" ) ) {
        fWarnForUnexpecedBlock = true;
	i++;
	continue;
      }

      HLTError("Unknown option '%s'", argv[i] );
      return EINVAL;
	  
    }//while loop

    int lutline = fTrigRec->GetLutLine();
    AliHLTMUONHitReconstructor::DHLTLut* lookupTable = new AliHLTMUONHitReconstructor::DHLTLut[lutline];
    if(!ReadLookUpTable(lookupTable,lutFileName)){
      HLTError("Failed to read lut, lut cannot be read");
      return ENOENT ; /* No such file or directory */
    }else{
      
      fTrigRec->LoadLookUpTable(lookupTable,fDDL+AliHLTMUONTriggerReconstructor::GetkDDLOffSet());

      AliHLTMUONTriggerReconstructor::RegToLoc regToLocMap[128]; // 16(locCard)*8(regCard)
      if(!ReadRegToLocMap(regToLocMap,reglocFileName)){
	HLTError("Failed to read RegToLocMap file");
	return ENOENT ; /* No such file or directory */
      }

      if(!(fTrigRec->SetRegToLocCardMap(regToLocMap))){
	HLTError("Failed to assign RegToLocMap to TrigRec Class due to memory problem");
	return ENOMEM ; /*cannot allocate memory*/
      }
      
    }// reading lut

    delete []lookupTable;

    HLTInfo("Initialisation of DHLT Trigger Record Component is done");

    return 0;
}


int AliHLTMUONTriggerReconstructorComponent::DoDeinit()
{
	HLTInfo("Deinitialising DHLT Trigger Record Component");

	if(fTrigRec)
	{
		delete fTrigRec;
		fTrigRec = NULL;
	}
	return 0;
}


int AliHLTMUONTriggerReconstructorComponent::DoEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks, 
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		std::vector<AliHLTComponentBlockData>& outputBlocks
	)
{
	// Process an event
	unsigned long totalSize = 0; // Amount of memory currently consumed in bytes.

	HLTDebug("Processing event %llu with %u input data blocks.",
		evtData.fEventID, evtData.fBlockCnt
	);
	
	// Loop over all input blocks in the event and run the trigger DDL
	// reconstruction algorithm on the raw data.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
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

		if (blocks[n].fDataType != AliHLTMUONConstants::TriggerDDLRawDataType())
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
				HLTWarning("Received a data block of a type we can not handle: %s origin %s",
					static_cast<char*>(id), static_cast<char*>(origin)
				);
			else
				HLTDebug("Received a data block of a type we can not handle: %s origin %s",
					static_cast<char*>(id), static_cast<char*>(origin)
				);
			
			continue;
		}
		
		// Create a new output data block and initialise the header.
		AliHLTMUONTriggerRecordsBlockWriter block(outputPtr+totalSize, size-totalSize);
		if (not block.InitCommonHeader())
		{
			HLTError("There is not enough space in the output buffer for the new data block.",
				 " We require at least %u bytes, but have %u bytes left.",
				sizeof(AliHLTMUONTriggerRecordsBlockWriter::HeaderType),
				block.BufferSize()
			);
			break;
		}

		AliHLTUInt32_t totalDDLSize = blocks[n].fSize / sizeof(AliHLTUInt32_t);
		AliHLTUInt32_t ddlRawDataSize = totalDDLSize - fTrigRec->GetkDDLHeaderSize();
		AliHLTUInt32_t* buffer = reinterpret_cast<AliHLTUInt32_t*>(blocks[n].fPtr)
			+ fTrigRec->GetkDDLHeaderSize();
		AliHLTUInt32_t nofTrigRec = block.MaxNumberOfEntries();

		if (not fTrigRec->Run(buffer, ddlRawDataSize, block.GetArray(), nofTrigRec))
		{
			HLTError("Error while processing of trigger DDL reconstruction algorithm.");
			size = totalSize; // Must tell the framework how much buffer space was used.
			return EIO;
		}
		
		// nofTrigRec should now contain the number of triggers actually found
		// and filled into the output data block, so we can set this number.
		assert( nofTrigRec <= block.MaxNumberOfEntries() );
		block.SetNumberOfEntries(nofTrigRec);
		
		HLTDebug("Number of trigger records found is %d", nofTrigRec);
		
		// Fill a block data structure for our output block.
		AliHLTComponentBlockData bd;
		FillBlockData(bd);
		bd.fPtr = outputPtr;
		// This block's start (offset) is after all other blocks written so far.
		bd.fOffset = totalSize;
		bd.fSize = block.BytesUsed();
		bd.fDataType = AliHLTMUONConstants::TriggerRecordsBlockDataType();
		bd.fSpecification = blocks[n].fSpecification;
		outputBlocks.push_back(bd);
		
		HLTDebug("Created a new output data block at fPtr = %p,"
			  " with fOffset = %u (0x%.X) and fSize = %u bytes.", 
			bd.fPtr, bd.fOffset, bd.fOffset, bd.fSize
		);
		
		// Increase the total amount of data written so far to our output memory.
		totalSize += block.BytesUsed();
	}
	
	// Finally we set the total size of output memory we consumed.
	size = totalSize;
	return 0;
}


bool AliHLTMUONTriggerReconstructorComponent::ReadLookUpTable(AliHLTMUONHitReconstructor::DHLTLut* lookupTable, const char* lutpath)
{
  if (fDDL < 0 || fDDL >= 2){
    HLTError("DDL number is out of range");
    return false;
  }
  
  int lutLine = fTrigRec->GetLutLine();
  
  FILE* fin = fopen(lutpath, "r");
  if (fin == NULL){
    HLTError("Failed to open file: %s",lutpath);
    return false;
  }
  
  for(int i=0;i<lutLine;i++){
    fscanf(
	   fin,
	   "%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n",
	   &lookupTable[i].fIdManuChannel,
	   &lookupTable[i].fIX,
	   &lookupTable[i].fIY,
	   &lookupTable[i].fRealX,
	   &lookupTable[i].fRealY,
	   &lookupTable[i].fRealZ,
	   &lookupTable[i].fPcbZone,
	   &lookupTable[i].fPlane
	   );
  }
  
  fclose(fin);
  return true;
}


bool AliHLTMUONTriggerReconstructorComponent::ReadRegToLocMap(AliHLTMUONTriggerReconstructor::RegToLoc* regToLocMap,const char* reglocFileName)
{
  int iTrigDDL,iReg,iLoc,locId,switchWord,detElemId[4];
  int index;

  memset(regToLocMap,-1,128*sizeof(AliHLTMUONTriggerReconstructor::RegToLoc));

  char s[100];
  ifstream fin(reglocFileName);
  
  if(!fin){
    HLTError("Failed to open file %s",reglocFileName);
    return false;
  }

  while(fin.getline(s,100)){
    sscanf(s,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
	   &iTrigDDL,&iReg,&iLoc,&locId,&switchWord,&detElemId[0],&detElemId[1],&detElemId[2],&detElemId[3]);
    if(iTrigDDL==fDDL){
      index = iReg*16 + iLoc;
      regToLocMap[index].fTrigDDL = iTrigDDL ; 
      regToLocMap[index].fRegId = iReg ;
      regToLocMap[index].fLoc = iLoc ;
      regToLocMap[index].fLocId = locId ;  
      regToLocMap[index].fSwitch = switchWord ;
      for(int idet = 0; idet<4; idet++)
	regToLocMap[index].fDetElemId[idet] = detElemId[idet] ;
    }// if matches with fDDL
  }//file loop
  
  fin.close();
  return true;
}
