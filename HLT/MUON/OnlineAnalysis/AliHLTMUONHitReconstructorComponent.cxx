/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
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

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTMUONHitReconstructorComponent.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTLogging.h"
#include "AliHLTSystem.h"
#include "AliHLTDefinitions.h"
#include <stdlib.h>
#include <errno.h>

namespace
{
	// The global object used for automatic component registration, 
	// Note DO NOT use this component for calculation!
	AliHLTMUONHitReconstructorComponent gAliHLTMUONHitReconstructorComponent;
}

ClassImp(AliHLTMUONHitReconstructorComponent)


AliHLTMUONHitReconstructorComponent::AliHLTMUONHitReconstructorComponent()
  : 
  fHitRec(NULL)
{
  // see header file for class documentation
}

AliHLTMUONHitReconstructorComponent::~AliHLTMUONHitReconstructorComponent()
{
  // see header file for class documentation

}

int AliHLTMUONHitReconstructorComponent::DoInit( int argc, const char** argv ){
  // see header file for class documentation

  fHitRec = new AliHLTMUONHitReconstructor();

  HLTInfo("dHLT hitrec");
  if (argc==0 && argv==NULL) {
    Logging( kHLTLogError, "AliHLTMUONHitReconstructorComponent::DoInit", "Arguments missing", " no arguments" );
    // this is just to get rid of the warning "unused parameter"
  }

  Int_t i = 0;
  char lutFileName[500],buspatchFileName[500];
  int iDDL = -1;
  while(i<argc){
    if ( !strcmp( argv[i], "lut" ) ) {
      if ( argc <= i+1 ) {
	Logging( kHLTLogError, "AliHLTMUONHitReconstructorComponent::DoInit", "Missing LookupTable filename", "LookupTable filename not specified" );
	return EINVAL; /* Invalid argument */ 
      }
      
      sprintf(lutFileName,"%s",argv[i+1]);

      i += 2;
      continue;
    }// lut argument

    if ( !strcmp( argv[i], "ddl" ) ) {
      if ( argc <= i+1 ) {
	Logging( kHLTLogError, "AliHLTMUONHitReconstructorComponent::DoInit", "Missing DDL argument", "DDL number not specified" );
	HLTError("AliHLTMUONHitReconstructorComponent::DoInit : DDL number is not specified ");
	return EINVAL;  /* Invalid argument */
      }

      iDDL = atoi(argv[i+1]);

      i += 2;
      continue;
    }// ddl argument

    if ( !strcmp( argv[i], "buspatchmap" ) ) {
      if ( argc <= i+1 ) {
	Logging( kHLTLogError, "AliHLTMUONHitReconstructorComponent::DoInit", "Missing buspatch filename", "buspatch filename not specified" );
	return EINVAL; /* Invalid argument */
      }
      
      sprintf(buspatchFileName,"%s",argv[i+1]);

      i += 2;
      continue;
    }// buspatch argument



  }// end of while loop

  int lutline = fHitRec->GetLutLine(iDDL);
  AliHLTMUONHitReconstructor::DHLTLut* lookupTable = new AliHLTMUONHitReconstructor::DHLTLut[lutline];
  if(!ReadLookUpTable(lookupTable,lutFileName,iDDL)){
    Logging(kHLTLogInfo, "AliHLTMUONHitReconstructorComponent::DoInit", "Failed to read lut", "lut cannot be read, DoInit");
    return ENOENT ; /* No such file or directory */
  }else{
    for(int i = 0;i<lutline; i++){
      //#ifdef MY_DEBUG
//       printf("%d\t%d\t%d\t%f\t%f\t%f\t%d\t\n",
// 	     lookupTable[i].fIdManuChannel,
// 	     lookupTable[i].fIX,
// 	     lookupTable[i].fIY,
// 	     lookupTable[i].fRealX,
// 	     lookupTable[i].fRealY,
// 	     lookupTable[i].fRealZ,
// 	     lookupTable[i].fPcbZone,
// 	     lookupTable[i].fPlane
// 	     );
      //#endif
    }

    BusToDetElem busToDetElem;
    if(!ReadBusPatchToDetElemFile(busToDetElem,buspatchFileName)){
      Logging(kHLTLogInfo, "AliHLTMUONHitReconstructorComponent::DoInit", "Failed to read buspatchmap", "buspatchmap cannot be read, DoInit");
      return ENOENT ; /* No such file or directory */
    }

    fHitRec->SetBusToDetMap(busToDetElem);
    fHitRec->LoadLookUpTable(lookupTable,iDDL);

  }// reading lut

  return 0;
}

int AliHLTMUONHitReconstructorComponent::DoDeinit(){
  // see header file for class documentation
  if(fHitRec)
    delete fHitRec;
  HLTInfo("dHLT hitrec");
  return 0;
}

int AliHLTMUONHitReconstructorComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
				      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ) {
  // see header file for class documentation
  HLTInfo("dHLT hitrec");

//   if (evtData.fStructSize==0 && blocks==NULL && trigData.fStructSize==0 &&
//       outputPtr==0 && size==0)
//   {
//     outputBlocks.clear();
//     // this is just to get rid of the warning "unused parameter"
//  }

  unsigned long totalSize = 0;
  unsigned long mySize;
  //cout<<"Block Count : "<<evtData.fBlockCnt<<endl;
  for(UInt_t i=0;i<evtData.fBlockCnt;i++){

    cout<<"0: totalsize : "<<totalSize<<"\tkAliHLTBlockAlignment :"<<kAliHLTBlockAlignment<<"\t size :"<<size<<endl;

    // Align the beginning of this  block to the required value.
    // 	  if ( totalSize % kAliHLTBlockAlignment ){
    // 	    totalSize += kAliHLTBlockAlignment-(totalSize % kAliHLTBlockAlignment);
    // 	  }

    if ( totalSize > size )
      break;
    
    cout<<"1: totalsize : "<<totalSize<<"\tkAliHLTBlockAlignment :"<<kAliHLTBlockAlignment<<endl;
    // Determine the size we should use for the output for this block (the input block's size times the relative output size)
    

    int totalDDLSize = blocks[i].fSize/4;
    int ddlRawDataSize = totalDDLSize - AliHLTMUONHitReconstructor::GetkDDLHeaderSize();
    int ddlHeader[8];
    memcpy((char *) & ddlHeader,blocks[i].fPtr,(size_t)4*(AliHLTMUONHitReconstructor::GetkDDLHeaderSize()));

//     for(int j=0;j<8;j++)
//       HLTDebug("ddlHeader[%d] : %d\n",j,ddlHeader[j]);

    int* buffer = new int[ddlRawDataSize];
    memcpy((int*)buffer,((int*)blocks[i].fPtr + AliHLTMUONHitReconstructor::GetkDDLHeaderSize()),(sizeof(int)*ddlRawDataSize));

//     for(int j=0;j<ddlRawDataSize;j++)
//       HLTDebug("buffer[%d] : %x\n",j,buffer[j]);

    
    AliHLTMUONRecHitStruct recHit[300];
    int nofHit = 300;
   
    if(! (fHitRec->Run(buffer,&ddlRawDataSize,recHit,&nofHit))){
      cerr <<"AliHLTMUONHitReconstructorComponent::DoEvent : ERROR In Processing of HitRec Algo "<< endl;
      return EIO;
    }

    mySize = sizeof(AliHLTMUONRecHitStruct)*nofHit;
    
    HLTInfo("mySize set (1) mySize == %lu B - blocks[%lu].fSize == %lu", 
	     mySize, i, blocks[i].fSize);
    
    // Check how much space we have left and adapt this output block's size accordingly.
    if ( totalSize + mySize > size )
      mySize = size-totalSize;
    
    Logging( kHLTLogInfo, "HLT::Dummy::DoEvent", "mySize set (2)", "mySize == %lu B - totalSize == %lu - size == %lu", 
		 mySize, totalSize, size );

    if ( mySize<=0 )
      continue; // No room left to write a further block.

    // Now copy the input block
    unsigned long copied = 0;
    // First copy all full multiples of the input block
    while ( copied+blocks[i].fSize <= mySize )
      {
	Logging( kHLTLogInfo, "0 : HLT::Dummy::DoEvent", "Copying", "Copying %lu B - Copied: %lu B - totalSize: %lu B", 
		 blocks[i].fSize, copied, totalSize );
	memcpy( outputPtr+totalSize+copied, blocks[i].fPtr, blocks[i].fSize );
	copied += blocks[i].fSize;
      }
    // And the copy the remaining fragment of the block
    Logging( kHLTLogInfo, "1 : HLT::Dummy::DoEvent", "Copying", "Copying %lu B - Copied: %lu B - totalSize: %lu B", 
	     mySize-copied, copied, totalSize );
    memcpy( outputPtr+totalSize+copied, blocks[i].fPtr, mySize-copied );
    Logging( kHLTLogInfo, "HLT::Dummy::DoEvent", "Copied", "Copied: %lu B - totalSize: %lu B", 
	     copied, totalSize );
    

    
    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset = totalSize;;
    bd.fSize = mySize;
    //bd.fPtr = (AliHLTUInt8_t*)(&recHit[0]);
    bd.fSpecification = blocks[i].fSpecification;
    outputBlocks.push_back( bd );
    
    totalSize += mySize;

    //    for(int j=0; j<nofHit; j++){
//       printf("%d\t\t%d\t\t%f\t%f\t%f\n",
// 	     i,outPtr->fRecPoint[j].fDetElemId,outPtr->fRecPoint[j].fX,
// 	     outPtr->fRecPoint[j].fY,outPtr->fRecPoint[j].fZ);
//       printf("1 : %d\t\t%d\t\t%f\t%f\t%f\n",
// 	     i,(((AliHLTMUONRecHitStruct*)bd.fPtr) + j)->fDetElemId,
// 	     (((AliHLTMUONRecHitStruct*)bd.fPtr) + j)->fX,
// 	     (((AliHLTMUONRecHitStruct*)bd.fPtr) + j)->fY,
// 	     (((AliHLTMUONRecHitStruct*)bd.fPtr) + j)->fZ);
//    }// nof RecHits


    if(totalSize > size){
      cout<<"size : "<<size<<"\t  totalSize: "<<totalSize<<"\t mySize :"<<mySize<<endl; 
      Logging( kHLTLogError, "AliHLTMUONHitReconstructorComponent::DoEvent", "Size exceeds the quota", "Size Exceeds the maximum quota" );
      return EMSGSIZE;
    }

  }// block loop

  size = totalSize;
  
  return 0;
}


bool AliHLTMUONHitReconstructorComponent::ReadLookUpTable(
		AliHLTMUONHitReconstructor::DHLTLut* lookupTable, const char* lutpath, int iDDL
	)
{
// Reads in a lookup table for the hit reconstruction algorithm from file.

	if (iDDL < AliHLTMUONHitReconstructor::GetkDDLOffSet() ||
		iDDL >= AliHLTMUONHitReconstructor::GetkDDLOffSet() + AliHLTMUONHitReconstructor::GetkNofDDL())
	{
// 		Logging(kHLTLogError, "AliHLTMUONHitReconstructorComponent::ReadLookUpTable", "Invalid DDL")
// 			<< "DDL number is out of range (must be " << AliHLTLog::kDec << HLTMUONHitRec::fgkDDLOffSet
// 			<< " <= iDDL < " << AliHLTLog::kDec
// 			<< HLTMUONHitRec::fgkDDLOffSet + HLTMUONHitRec::fgkNofDDL
// 			<< ENDLOG;
		return false;
	}

	int lutLine = fHitRec->GetLutLine(iDDL);
// 	cout<<"LutLine :"<<lutLine<<endl;
	FILE* fin = fopen(lutpath, "r");
	if (fin == NULL)
	{
// 		LOG(kHLTLogError, "AliHLTMUONHitReconstructorComponent::ReadLookUpTable", "I/O error")
// 			<< "Failed to open file: " << lutpath << ENDLOG;
		return false;
  	}

#	ifdef DEBUG
// 	Logging(kHLTLogError, "AliHLTMUONHitReconstructorComponent::ReadLookUpTable", "Trace")
// 		<< "Reading LUT file: " << lutpath << ENDLOG;
#	endif

	for(int i=0;i<lutLine;i++)
	{
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

// implement this as well.

bool AliHLTMUONHitReconstructorComponent::ReadBusPatchToDetElemFile(
		BusToDetElem& busToDetElem, const char* buspatchmappath
	)
{
// Loads the bus patch to detector element ID map from an ASCII file. 
	
	char getLine[80];
	char temp;
	int detElem, minBusPatch, maxBusPatch;

	FILE* fin = fopen(buspatchmappath, "r");
	if (fin == NULL)
	{
// 		Logging(kHLTLogError, "AliHLTMUONHitReconstructorComponent::ReadBusPatchToDetElemFile", "I/O error")
// 			<< "Failed to open file: " << buspatchmappath << ENDLOG;
		return false;
  	}

#	ifdef DEBUG
// 	Logging(kHLTLogError, "AliHLTMUONHitReconstructorComponent::ReadBusPatchToDetElemFile", "Trace")
// 		<< "Reading bus patch mapping file: " << buspatchmappath << ENDLOG;
#	endif

	while (feof(fin)==0)
	{
		fgets(getLine,80,fin);
		sscanf(getLine, "%d\t%d %c %d", &detElem, &minBusPatch, &temp, &maxBusPatch);
		if (detElem >= 700 && detElem <= 1025)
		{
			for(int i = minBusPatch; i <= maxBusPatch; i++)
				busToDetElem[i] = detElem;
		} // detElem condn
	} // while loop for file

	fclose(fin);
	return true;
}
