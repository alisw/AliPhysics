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
    @brief  A processing component for the dHLT TrigRec. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTSystem.h"
#include "AliHLTMUONTriggerReconstructorComponent.h"
#include "AliHLTDefinitions.h"
#include <stdlib.h>
#include <errno.h>

// this is a global object used for automatic component registration, do not use this
AliHLTMUONTriggerReconstructorComponent gAliHLTMUONTriggerReconstructorComponent;

ClassImp(AliHLTMUONTriggerReconstructorComponent)
    
AliHLTMUONTriggerReconstructorComponent::AliHLTMUONTriggerReconstructorComponent()
  :
  fOutputPercentage(100), // By default we copy to the output exactly what we got as input
  fDDLDir(""),
  fDDL(0),
  fTrigRec(NULL)
    {
    }

AliHLTMUONTriggerReconstructorComponent::~AliHLTMUONTriggerReconstructorComponent()
    {
    }

const char* AliHLTMUONTriggerReconstructorComponent::GetComponentID()
    {
    return "MUONTrigRec"; // The ID of this component
    }

void AliHLTMUONTriggerReconstructorComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
    {
      /* in order to be backward compatible we have to keep the old code, at
       * least for a while. Remember to use the new const kAliHLTVoidDataType
       * if you are using a more recent AliRoot version (from Jan 07)
       list.push_back(kAliHLTAnyDataType); // We do not have any requirements for our input data type(s).
      */

      list.clear();
      list.push_back( AliHLTMUONConstants::TriggerDDLRawDataType() );
    }

AliHLTComponentDataType AliHLTMUONTriggerReconstructorComponent::GetOutputDataType()
    {
      /* in order to be backward compatible we have to keep the old code, at
       * least for a while. Remember to use the new const kAliHLTVoidDataType
       * if you are using a more recent AliRoot version (from Jan 07)
      return kAliHLTVoidDataType;
      */
      return AliHLTMUONConstants::TriggerRecordsBlockDataType();
    }

void AliHLTMUONTriggerReconstructorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
    constBase = 0;
    inputMultiplier = ((double)fOutputPercentage)/100.0;
    }



// Spawn function, return new instance of this class
AliHLTComponent* AliHLTMUONTriggerReconstructorComponent::Spawn()
    {
    return new AliHLTMUONTriggerReconstructorComponent;
    }

int AliHLTMUONTriggerReconstructorComponent::DoInit( int argc, const char** argv )
{
    // perform initialization. We check whether our relative output size is specified in the arguments.
      
      fTrigRec = new AliHLTMUONTriggerReconstructor();
      
      HLTInfo("dHLT trigrec");
      if (argc==0 && argv==NULL) {
	Logging( kHLTLogError, "AliHLTMUONTriggerReconstructorComponent::DoInit", "Arguments missing", " no arguments" );
	// this is just to get rid of the warning "unused parameter"
      }

      //Int_t i = 0;
      char lutFileName[500], ddlDir[500];

      fOutputPercentage = 100;
      int i = 0;
      char* cpErr;
      while ( i < argc )
	{
	  Logging( kHLTLogDebug, "HLT::MUONTrigRec::DoInit", "Arguments", "argv[%d] == %s", i, argv[i] );
	  if ( !strcmp( argv[i], "output_percentage" ) )
	    {
	      if ( i+1>=argc )
		{
		  Logging(kHLTLogError, "HLT::MUONTrigRec::DoInit", "Missing Argument", "Missing output_percentage parameter");
		  return ENOTSUP;
		}
	      Logging( kHLTLogDebug, "HLT::MUONTrigRec::DoInit", "Arguments", "argv[%d+1] == %s", i, argv[i+1] );
	      fOutputPercentage = strtoul( argv[i+1], &cpErr, 0 );
	      if ( *cpErr )
		{
		  Logging(kHLTLogError, "HLT::MUONTrigRec::DoInit", "Wrong Argument", "Cannot convert output_percentage parameter '%s'", argv[i+1] );
		  return EINVAL;
		}
	      Logging( kHLTLogInfo, "HLT::MUONTrigRec::DoInit", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
	      i += 2;
	    continue;
	    }
	

	  if ( !strcmp( argv[i], "lut" ) ) {
	    if ( argc <= i+1 ) {
	      Logging( kHLTLogError, "AliHLTMUONTriggerReconstructorComponent::DoInit", "Missing LookupTable filename", "LookupTable filename not specified" );
	      return EINVAL; /* Invalid argument */ 
	    }
	    
	    sprintf(lutFileName,"%s",argv[i+1]);
	    
	    i += 2;
	    continue;
	  }// lut argument
	  
	  
	  if ( !strcmp( argv[i], "ddl" ) ) {
	    if ( argc <= i+1 ) {
	      Logging( kHLTLogError, "AliHLTMUONTriggerReconstructorComponent::DoInit", "Missing DDL argument", "DDL number not specified" );
	      HLTError("AliHLTMUONTriggerReconstructorComponent::DoInit : DDL number is not specified ");
	      return EINVAL;  /* Invalid argument */
	    }
	    
	    fDDL = atoi(argv[i+1]);
	    
	    i += 2;
	    continue;
	  }// ddl argument
	  

	  if ( !strcmp( argv[i], "rawdir" ) ) {
	    if ( argc <= i+1 ) {
	      Logging( kHLTLogError, "AliHLTMUONTriggerReconstructorComponent::DoInit", "Missing DDL directory", "DDL directory not specified" );
	      HLTError("AliHLTMUONTriggerReconstructorComponent::DoInit : DDL directory is not specified ");
	      return EINVAL;  /* Invalid argument */
	    }
	    
	    fDDLDir = argv[i+1] ;
	    
	    i += 2;
	    continue;
	  }// ddl directory argument

	  Logging(kHLTLogError, "HLT::MUONTrigRec::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
	  return EINVAL;
	  
	}//while loop

    int lutline = fTrigRec->GetLutLine();
    AliHLTMUONHitReconstructor::DHLTLut* lookupTable = new AliHLTMUONHitReconstructor::DHLTLut[lutline];
    if(!ReadLookUpTable(lookupTable,lutFileName)){
      Logging(kHLTLogInfo, "AliHLTMUONTriggerReconstructorComponent::DoInit", "Failed to read lut", "lut cannot be read, DoInit");
      return ENOENT ; /* No such file or directory */
    }else{
      
      fTrigRec->LoadLookUpTable(lookupTable,fDDL);
      
    }// reading lut

    delete []lookupTable;

    return 0;
}

int AliHLTMUONTriggerReconstructorComponent::DoDeinit()
    {
      if(fTrigRec)
	delete fTrigRec;
      HLTInfo("dHLT trigrec");
      return 0;
  
    return 0;
    }

int AliHLTMUONTriggerReconstructorComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
				      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
    {
      // Process an event
    unsigned long totalSize = 0;
    Logging( kHLTLogInfo, "HLT::MUONTrigRec::DoEvent", "Output percentage set", "Output percentage set to %lu %% and totalSize %lu", fOutputPercentage,totalSize );

    // Loop over all input blocks in the event
    for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
      {

	if ( totalSize > size )
	    break;

	int totalDDLSize = blocks[n].fSize/4;
	int  ddlRawDataSize = totalDDLSize - AliHLTMUONTriggerReconstructor::GetkDDLHeaderSize();
	//cout<<"ddlRawDataSize :"<<ddlRawDataSize<<endl;
	int *buffer = (int *)((int *)blocks[n].fPtr + AliHLTMUONTriggerReconstructor::GetkDDLHeaderSize()) ;
	
	AliHLTMUONTriggerRecordStruct trigRecord[300];
	int nofTrigRec = 300;
	
 	if(! (fTrigRec->Run(buffer,&ddlRawDataSize,trigRecord,&nofTrigRec))){
	  HLTError("ERROR In Processing of TrigRec Algo ");
	  return EIO;
	}

// 	if(! (fTrigRec->Run((int)evtData.fEventID,fDDL,trigRecord,&nofTrigRec))){
// 	  HLTError("ERROR In Processing of TrigRec Algo ");
// 	  return EIO;
//	}
	
	unsigned long mySize = sizeof(AliHLTMUONTriggerRecordStruct)*nofTrigRec;
    
//	cout<<"nofHit "<<nofHit<<endl;
// 	for(int i=0;i<nofHit;i++)
// 	  cout<<"\t 0 : recHit["<<i<<"].fX :"<<recHit[i].fX
// 	      <<"  recHit["<<i<<"].fY :"<<recHit[i].fY
// 	      <<"  recHit["<<i<<"].fZ :"<<recHit[i].fZ
// 	      <<"  recHit["<<i<<"].fDetElemId :"<<recHit[i].fDetElemId
// 	      <<endl;

	//unsigned long mySize = (blocks[n].fSize * fOutputPercentage) / 100;

	Logging( kHLTLogInfo, "HLT::MUONTrigRec::DoEvent", "mySize set (1)", "mySize == %lu B - blocks[%lu].fSize == %lu - fOutputPercentage == %lu", 
		 mySize, n, blocks[n].fSize, fOutputPercentage );

	// Check how much space we have left and adapt this output block's size accordingly.
	if ( totalSize + mySize > size )
	    mySize = size-totalSize;

	Logging( kHLTLogInfo, "HLT::MUONTrigRec::DoEvent", "mySize set (2)", "mySize == %lu B - totalSize == %lu - size == %lu", 
		 mySize, totalSize, size );

	if ( mySize<=0 )
	    continue; // No room left to write a further block.

	// Now copy the input block
	unsigned long copied = 0;
	// First copy all full multiples of the input block

	// And the copy the remaining fragment of the block
	Logging( kHLTLogInfo, "1 : HLT::MUONTrigRec::DoEvent", "Copying", "Copying %lu B - Copied: %lu B - totalSize: %lu B", 
		 mySize-copied, copied, totalSize );
	//memcpy( outputPtr+totalSize+copied, blocks[n].fPtr, mySize-copied );
	memcpy( outputPtr+totalSize+copied, &trigRecord[0], mySize);
	Logging( kHLTLogInfo, "HLT::MUONTrigRec::DoEvent", "Copied", "Copied: %lu B - totalSize: %lu B", 
		 copied, totalSize );
	// Fill a block data structure for our output block.
	AliHLTComponentBlockData ob;
	// Let the structure be filled with the default values.
	// This takes care of setting the shared memory and data type values to default values,
	// so that they can be filled in by the calling code.
	FillBlockData( ob );
	// This block's start (offset) is after all other blocks written so far
	ob.fOffset = totalSize;
	// the size of this block's data.
	ob.fSize = mySize;
	// The specification of the data is copied from the input block.
	ob.fSpecification = blocks[n].fSpecification;
	// The data type is set automatically to the component's specified output data type.
	// Place this block into the list of output blocks
	outputBlocks.push_back( ob );
	// Increase the total amount of data written so far to our output memory
	totalSize += mySize;
	}
    // Finally we set the total size of output memory we consumed.
    size = totalSize;

    return 0;
    }

bool AliHLTMUONTriggerReconstructorComponent::ReadLookUpTable(AliHLTMUONHitReconstructor::DHLTLut* lookupTable, const char* lutpath)
{
  if (fDDL < AliHLTMUONTriggerReconstructor::GetkDDLOffSet() ||
      fDDL >= AliHLTMUONTriggerReconstructor::GetkDDLOffSet() + AliHLTMUONTriggerReconstructor::GetkNofDDL()){
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
	   &lookupTable[i].fPlane,
	   &lookupTable[i].fPcbZone
	   );
  }
  
  fclose(fin);
  return true;
}

// implement this as well.

