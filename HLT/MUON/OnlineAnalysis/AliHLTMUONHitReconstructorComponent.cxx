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

///*
//
//  The HitRec Component is designed to deal the rawdata inputfiles to findout the 
//  the reconstructed hits. The output is send to the output block for further 
//  processing.
//
//  Author : Indranil Das ( indra.das@saha.ac.in || indra.ehep@gmail.com )
// 
//*/

#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONHitReconstructorComponent.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONDataBlockWriter.h"
#include "AliHLTLogging.h"
#include "AliHLTSystem.h"
#include "AliHLTDefinitions.h"
#include <cstdlib>
#include <cerrno>
#include <cassert>

namespace
{
	// The global object used for automatic component registration, 
	// Note DO NOT use this component for calculation!
	AliHLTMUONHitReconstructorComponent gAliHLTMUONHitReconstructorComponent;
}

ClassImp(AliHLTMUONHitReconstructorComponent)


AliHLTMUONHitReconstructorComponent::AliHLTMUONHitReconstructorComponent() :
	fHitRec(NULL),
	fDDLDir(""),
	fDDL(0),
	fReaderType(false),
	fWarnForUnexpecedBlock(false)
{
}


AliHLTMUONHitReconstructorComponent::~AliHLTMUONHitReconstructorComponent()
{
}

const char* AliHLTMUONHitReconstructorComponent::GetComponentID()
{
	return AliHLTMUONConstants::HitReconstructorId();
}


void AliHLTMUONHitReconstructorComponent::GetInputDataTypes( std::vector<AliHLTComponentDataType>& list)
{
	list.clear();
	list.push_back( AliHLTMUONConstants::TrackingDDLRawDataType() );
}


AliHLTComponentDataType AliHLTMUONHitReconstructorComponent::GetOutputDataType()
{
	return AliHLTMUONConstants::RecHitsBlockDataType();
}


void AliHLTMUONHitReconstructorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
	constBase = sizeof(AliHLTMUONRecHitsBlockWriter::HeaderType);
	inputMultiplier = 1;
}


// Spawn function, return new instance of this class
AliHLTComponent* AliHLTMUONHitReconstructorComponent::Spawn()
{
	return new AliHLTMUONHitReconstructorComponent;
}


int AliHLTMUONHitReconstructorComponent::DoInit(int argc, const char** argv)
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
     
  HLTInfo("Initialising DHLT HitReconstruction Component");

  fHitRec = new AliHLTMUONHitReconstructor();
  fWarnForUnexpecedBlock = false;

  // this is to get rid of the warning "unused parameter"
  if (argc==0 && argv==NULL) {
    HLTError("Arguments missing", " no arguments" );
  }

  char lutFileName[500],buspatchFileName[500];

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
	  
	  
      if ( !strcmp( argv[i], "buspatchmap" ) ) {
	if ( argc <= i+1 ) {
	  HLTError("Buspatch filename not specified" );
	  return EINVAL; /* Invalid argument */
	}
	    
	sprintf(buspatchFileName,"%s",argv[i+1]);
	
	i += 2;
	continue;
      }// buspatch argument

      if ( !strcmp( argv[i], "rawreader" ) ) {
	fReaderType = true; // true when using rawreader for standalone it is set to false.
	i += 1;
	continue;
      }
	  
      if ( !strcmp( argv[i], "-warn_on_unexpected_block" ) ) {
        fWarnForUnexpecedBlock = true;
	i++;
	continue;
      }

      HLTError("Unknown option '%s'", argv[i] );
      return EINVAL;
      
    }//while loop

  int lutline = fHitRec->GetLutLine(fDDL);
  AliHLTMUONHitReconstructor::DHLTLut* lookupTable = new AliHLTMUONHitReconstructor::DHLTLut[lutline];
  if(!ReadLookUpTable(lookupTable,lutFileName)){
    HLTError("Failed to read lut, lut cannot be read, DoInit");
    return ENOENT ; /* No such file or directory */
  }else{
    
    BusToDetElem busToDetElem;
    BusToDDL busToDDL;
    if(!ReadBusPatchToDetElemFile(busToDetElem,busToDDL,buspatchFileName)){
      HLTError("Failed to read buspatchmap, buspatchmap cannot be read, DoInit");
      return ENOENT ; /* No such file or directory */
    }
    
    fHitRec->SetBusToDetMap(busToDetElem);
    fHitRec->SetBusToDDLMap(busToDDL);
    fHitRec->LoadLookUpTable(lookupTable,fDDL);
    
  }// reading lut

  delete []lookupTable;
  
  HLTInfo("Initialisation of DHLT HitReconstruction Component is done");
  
  return 0;
}


int AliHLTMUONHitReconstructorComponent::DoDeinit()
{
  if(fHitRec)
    delete fHitRec;
  
  HLTInfo(" Deinitialising DHLT HitReconstruction Component");
  
  return 0;
}


int AliHLTMUONHitReconstructorComponent::DoEvent(
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

		if (blocks[n].fDataType != AliHLTMUONConstants::TrackingDDLRawDataType())
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
				HLTWarning("Received a data block of a type we cannot handle: '%s' origin: '%s'",
					static_cast<char*>(id), static_cast<char*>(origin)
				);
			else
				HLTDebug("Received a data block of a type we cannot handle: '%s' origin: '%s'",
					static_cast<char*>(id), static_cast<char*>(origin)
				);
			
			continue;
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

		HLTDebug("=========== Dumping DDL payload buffer ==========");
		for (AliHLTUInt32_t j = 0; j < totalDDLSize; j++)
			HLTDebug("buffer[%d] : %x",j,buffer[j]);
		HLTDebug("================== End of dump =================");

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


bool AliHLTMUONHitReconstructorComponent::ReadLookUpTable(AliHLTMUONHitReconstructor::DHLTLut* lookupTable, const char* lutpath)
{
  if (fDDL < AliHLTMUONHitReconstructor::GetkDDLOffSet() ||
      fDDL >= AliHLTMUONHitReconstructor::GetkDDLOffSet() + AliHLTMUONHitReconstructor::GetkNofDDL()){
    HLTError("DDL number is out of range");
    return false;
  }
  
  int lutLine = fHitRec->GetLutLine(fDDL);
  
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


bool AliHLTMUONHitReconstructorComponent::ReadBusPatchToDetElemFile(BusToDetElem& busToDetElem, BusToDDL& busToDDL, const char* buspatchmappath)
{
  char getLine[80];
  char temp;
  int detElem, minBusPatch, maxBusPatch, ddl;
  
  FILE* fin = fopen(buspatchmappath, "r");
  if (fin == NULL){
    HLTError("Failed to open file: %s",buspatchmappath);
    return false;
  }
  
  while (feof(fin)==0){
    fgets(getLine,80,fin);
    sscanf(getLine, "%d\t%d %c %d\t%d", &detElem, &minBusPatch, &temp, &maxBusPatch,&ddl);
    if (detElem >= 700 && detElem <= 1025){
      
      for(int i = minBusPatch; i <= maxBusPatch; i++){
	busToDetElem[i] = detElem;
	busToDDL[i] = ddl;
      }//for loop
    } // detElem condn
  } // while loop for file
  
  fclose(fin);
  return true;
}
