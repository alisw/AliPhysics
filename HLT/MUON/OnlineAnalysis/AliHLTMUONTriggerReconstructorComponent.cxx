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

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTSystem.h"
#include "AliHLTMUONTriggerReconstructorComponent.h"
#include "AliHLTMUONTriggerReconstructor.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTMUONConstants.h"

#include <stdlib.h>
#include <errno.h>

namespace
{
	// This is a global object used for automatic component registration,
	// do not use this for calculation.
	AliHLTMUONTriggerReconstructorComponent gAliHLTMUONTriggerReconstructorComponent;
} // end of namespace


ClassImp(AliHLTMUONTriggerReconstructorComponent)
    
    
AliHLTMUONTriggerReconstructorComponent::AliHLTMUONTriggerReconstructorComponent()
  :
  fTrigRec(NULL),
  fDDLDir(""),
  fDDL(0)
{
}


AliHLTMUONTriggerReconstructorComponent::~AliHLTMUONTriggerReconstructorComponent()
{
}


const char* AliHLTMUONTriggerReconstructorComponent::GetComponentID()
{
  return "MUONTrigRec"; // The ID of this component
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


void AliHLTMUONTriggerReconstructorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  constBase = 0;
  inputMultiplier = 1;
}


// Spawn function, return new instance of this class
AliHLTComponent* AliHLTMUONTriggerReconstructorComponent::Spawn()
{
  return new AliHLTMUONTriggerReconstructorComponent;
}


int AliHLTMUONTriggerReconstructorComponent::DoInit(int argc, const char** argv)
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  
  HLTInfo("Initialising DHLT Trigger Record Component");

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
  if(fTrigRec)
    delete fTrigRec;
  
  HLTInfo(" Deinitialising DHLT Trigger Record Component");
  
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
  unsigned long totalSize = 0;
  HLTDebug("Output percentage set to %lu and totalSize %lu",fOutputPercentage,totalSize );
    
  HLTDebug("Event : %d has : %lu  blocks",(int)evtData.fEventID,evtData.fBlockCnt);
  // Loop over all input blocks in the event
  for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
    {

      HLTDebug("block : %d, block rawData : %p, block.fSize (bytes) : %d, blocks.fDataType.fID : %s, blocks.fDataType.fOrigin  : %s, required type : %s\n",
	       n,blocks[n].fPtr,blocks[n].fSize,(char *)(blocks[n].fDataType.fID),
	      (char *)(blocks[n].fDataType.fOrigin,(char *)(AliHLTMUONConstants::TriggerDDLRawDataType().fID)));
      

      if(strncmp((char *)(blocks[n].fDataType.fID),(char *)(AliHLTMUONConstants::TriggerDDLRawDataType().fID),kAliHLTComponentDataTypefIDsize)) continue;

      if ( totalSize > size )
	break;
      
      int totalDDLSize = blocks[n].fSize/sizeof(int);
      int  ddlRawDataSize = totalDDLSize - fTrigRec->GetkDDLHeaderSize();

      int *buffer = (int *)((int *)blocks[n].fPtr + fTrigRec->GetkDDLHeaderSize()) ;
      
      AliHLTMUONTriggerRecordStruct trigRecord[300];
      int nofTrigRec = 300;
	
      if(! (fTrigRec->Run(buffer,&ddlRawDataSize,&trigRecord[0],&nofTrigRec))){
	HLTError("ERROR In Processing of TrigRec Algo ");
	return EIO;
      }
      
      // 	if(! (fTrigRec->Run((int)evtData.fEventID,fDDL,trigRecord,&nofTrigRec))){
      // 	  HLTError("ERROR In Processing of TrigRec Algo ");
      // 	  return EIO;
      //	}
	
      unsigned long mySize = sizeof(AliHLTMUONTriggerRecordStruct)*nofTrigRec;
    
      HLTDebug("Number record found is %d",nofTrigRec);
//       for(int ihit=0;ihit<nofTrigRec;ihit++)
// 	cout<<"\tdetelem : "<<trigRecord[ihit].fId
// 	    <<"\t"<<trigRecord[ihit].fHit[0].fX
// 	    <<"\t"<<trigRecord[ihit].fHit[0].fY
// 	    <<"\t"<<trigRecord[ihit].fHit[0].fZ
// 	    <<endl;

	// Check how much space we have left and adapt this output block's size accordingly.
	if ( totalSize + mySize > size )
	    mySize = size-totalSize;

	Logging( kHLTLogDebug, "AliHLTMUONTriggerReconstructor::DoEvent", "mySize set (2)", "mySize == %lu B - totalSize == %lu - size == %lu", 
		 mySize, totalSize, size );

	if ( mySize<=0 )
	    continue; // No room left to write a further block.

	// Now copy the input block
	unsigned long copied = 0;
	// First copy all full multiples of the input block

	// And the copy the remaining fragment of the block
	Logging( kHLTLogDebug, "AliHLTMUONTriggerReconstructor::DoEvent", "Copying", "Copying %lu B - Copied: %lu B - totalSize: %lu B", 
		 mySize-copied, copied, totalSize );
	//memcpy( outputPtr+totalSize+copied, blocks[n].fPtr, mySize-copied );
	memcpy( outputPtr+totalSize+copied, &trigRecord[0], mySize);
	Logging( kHLTLogDebug, "AliHLTMUONTriggerReconstructor::DoEvent", "Copied", "Copied: %lu B - totalSize: %lu B", 
		 copied, totalSize );
	
	// Fill a block data structure for our output block.
	AliHLTComponentBlockData bd;
	FillBlockData(bd);
	bd.fPtr = outputPtr;
	// This block's start (offset) is after all other blocks written so far.
	bd.fOffset = totalSize;
	bd.fSize = mySize;
	bd.fDataType = AliHLTMUONConstants::TriggerRecordsBlockDataType();
	bd.fSpecification = blocks[n].fSpecification;
	outputBlocks.push_back(bd);
	
	// Increase the total amount of data written so far to our output memory
	totalSize += mySize;
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
