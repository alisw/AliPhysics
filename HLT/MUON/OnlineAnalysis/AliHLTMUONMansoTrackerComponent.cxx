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

/**
 *  @file   AliHLTMUONMansoTrackerComponent.cxx
 *  @author Artur Szostak <artursz@iafrica.com>
 *  @date   
 *  @brief  Implementation of AliHLTMUONMansoTrackerComponent class.
 */

#include "AliHLTMUONMansoTrackerComponent.h"
#include "Util/AliHLTMUONRecPoint.h"
#include <stdlib.h>
#include <errno.h>

namespace
{
	// this is a global object used for automatic component registration, do not use this
	AliHLTMUONMansoTrackerComponent gAliHLTMUONMansoTrackerComponent;
}


ClassImp(AliHLTMUONMansoTrackerComponent);


AliHLTMUONMansoTrackerComponent::AliHLTMUONMansoTrackerComponent()
  :
    fOutputPercentage(100) // By default we copy to the output exactly what we got as input
    {
    }

AliHLTMUONMansoTrackerComponent::~AliHLTMUONMansoTrackerComponent()
    {
    }

const char* AliHLTMUONMansoTrackerComponent::GetComponentID()
    {
    return "tracking"; // The ID of this component
    }

void AliHLTMUONMansoTrackerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
    {
      /* in order to be backward compatible we have to keep the old code, at
       * least for a while. Remember to use the new const kAliHLTVoidDataType
       * if you are using a more recent AliRoot version (from Jan 07)
       list.push_back(kAliHLTAnyDataType); // We do not have any requirements for our input data type(s).
      */

      AliHLTComponentDataType dt = 
	{ sizeof(AliHLTComponentDataType),
	  {'*','*','*','*','*','*','*','\0'},
	  {'*','*','*','\0'}};
       list.push_back(dt);
    }

AliHLTComponentDataType AliHLTMUONMansoTrackerComponent::GetOutputDataType()
    {
      /* in order to be backward compatible we have to keep the old code, at
       * least for a while. Remember to use the new const kAliHLTVoidDataType
       * if you are using a more recent AliRoot version (from Jan 07)
      return kAliHLTVoidDataType;
      */
      AliHLTComponentDataType dt = 
	{ sizeof(AliHLTComponentDataType),
	  {'\0','\0','\0','0','\0','\0','\0','\0'},
	  {'\0','\0','\0','\0'}};
      return dt;
    }

void AliHLTMUONMansoTrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
    constBase = 0;
    inputMultiplier = ((double)fOutputPercentage)/100.0;
    }



// Spawn function, return new instance of this class
AliHLTComponent* AliHLTMUONMansoTrackerComponent::Spawn()
    {
    return new AliHLTMUONMansoTrackerComponent;
    }

int AliHLTMUONMansoTrackerComponent::DoInit( int argc, const char** argv )
    {
    // perform initialization. We check whether our relative output size is specified in the arguments.
    fOutputPercentage = 100;
						
    fTracker = new AliHLTMUONMansoTracker();

    Logging(kHLTLogInfo, "dHLT", "Tracking", "hitrec, DoInit");
    if (argc==0 && argv==NULL) {
      // this is just to get rid of the warning "unused parameter"
    }

    int i = 0;
    char* cpErr;
    while ( i < argc )
      {
	Logging( kHLTLogDebug, "AliHLTMUONMansoTrackerComponent::DoInit", "Arguments", "argv[%d] == %s", i, argv[i] );
	if ( !strcmp( argv[i], "output_percentage" ) )
	    {
	    if ( i+1>=argc )
		{
		Logging(kHLTLogError, "AliHLTMUONMansoTrackerComponent::DoInit", "Missing Argument", "Missing output_percentage parameter");
		return ENOTSUP;
		}
	    Logging( kHLTLogDebug, "AliHLTMUONMansoTrackerComponent::DoInit", "Arguments", "argv[%d+1] == %s", i, argv[i+1] );
	    fOutputPercentage = strtoul( argv[i+1], &cpErr, 0 );
	    if ( *cpErr )
		{
		  Logging(kHLTLogError, "AliHLTMUONMansoTrackerComponent::DoInit", "Wrong Argument", "Cannot convert output_percentage parameter '%s'", argv[i+1] );
		  return EINVAL;
		}
	    Logging( kHLTLogInfo, "AliHLTMUONMansoTrackerComponent::DoInit", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
	    i += 2;
	    continue;
	    }
	Logging(kHLTLogError, "AliHLTMUONMansoTrackerComponent::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
	return EINVAL;
      }// while loop
    return 0;
    }

int AliHLTMUONMansoTrackerComponent::DoDeinit()
{
  if(fTracker)
    delete fTracker;
  return 0;
}

int AliHLTMUONMansoTrackerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
				      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  Logging( kHLTLogInfo, "AliHLTMUONMansoTrackerComponent::DoEvent", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
  // Process an event
  unsigned long totalSize = 0;
  AliHLTUInt32_t maxTracksPerEvent = 1024;

  fTracker->Reset();
  fTracker->SetTrackOutputBuffer((AliHLTMUONTrackPoints*)outputPtr,maxTracksPerEvent);
  
//  cout<<"BlockSize :"<<evtData.fBlockCnt<<" size :"<<outputBlocks.size()<<endl;
  // Loop over all input blocks in the event
  for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
    {
      // Align the beginning of this  block to the required value.
//       if ( totalSize % kAliHLTBlockAlignment )
// 	totalSize += kAliHLTBlockAlignment-(totalSize % kAliHLTBlockAlignment);

      if ( totalSize > size )
	break;
      // Determine the size we should use for the output for this block (the input block's size times the relative output size)
      unsigned long mySize = (blocks[n].fSize * fOutputPercentage) / 100;
      Logging( kHLTLogInfo, "AliHLTMUONMansoTrackerComponent::DoEvent", "mySize set (1)", "mySize == %lu B - blocks[%lu].fSize == %lu - fOutputPercentage == %lu - totalSize == %lu", 
	       mySize, n, blocks[n].fSize, fOutputPercentage, totalSize );

      Int_t noftrigData,nofhitData;
      if(n== (evtData.fBlockCnt - 1)){// for trigger
	UInt_t trigDataSize;
	memcpy(&trigDataSize,blocks[n].fPtr,sizeof(UInt_t));
//	cout<<"trigDataSize :"<<trigDataSize<<endl;
	UInt_t trigOffSet;
	memcpy(&trigOffSet,(UInt_t*)blocks[n].fPtr + 1,sizeof(UInt_t));
//	cout<<"trigOffSet :"<<trigOffSet<<endl;
	noftrigData = (trigDataSize -  sizeof(UInt_t))*sizeof(AliHLTUInt32_t)/sizeof(AliHLTMUONCoreTriggerRecord) ;
	fTrigData = new AliHLTMUONCoreTriggerRecord[noftrigData];

	for(Int_t i=0;i<noftrigData;i++){
	  AliHLTMUONCoreTriggerRecord record;
	  memcpy(&record,(UInt_t*)blocks[n].fPtr + 2 + i*(sizeof(AliHLTMUONCoreTriggerRecord))/4,sizeof(AliHLTMUONCoreTriggerRecord));
	  fTrigData[i] = record;

//	  cout<<" Sign : "<<fTrigData[i].fSign
//	      <<" Pt : "<<fTrigData[i].fPt
//	      <<"\t X1:"<<fTrigData[i].fStation1impact.X()<<" Y1 :"<<fTrigData[i].fStation1impact.Y()
//	      <<"\t X2:"<<fTrigData[i].fStation2impact.X()<<" Y2 :"<<fTrigData[i].fStation2impact.Y()
//	      << endl;
	}// for
	
	fTracker->FindTracks(fTrigData,noftrigData);
//	cout<<"Nof tracks found :"<<fTracker->TrackCount()<<endl;
//	cout<<"Z7 : "<<AliHLTMUONCoreMansoTracker::GetZ7()<<endl;

      }else{ // for hitrec
	UInt_t hitDataSize;
	memcpy(&hitDataSize,blocks[n].fPtr,sizeof(UInt_t));
//	cout<<"hitDataSize :"<<hitDataSize<<endl;
	nofhitData = hitDataSize*sizeof(AliHLTUInt32_t)/sizeof(AliHLTMUONCorePoint) ;
	
	Int_t chamber = n + 6;
	AliHLTMUONRecPoint *recHit = new AliHLTMUONRecPoint[nofhitData];
	for(Int_t i=0;i<nofhitData;i++){
	  AliHLTMUONCorePoint point;
	  memcpy(&point,(UInt_t*)blocks[n].fPtr + 1 + i*(sizeof(AliHLTMUONCorePoint))/4,sizeof(AliHLTMUONCorePoint));
//	  cout <<"chamber :"<<chamber<<"\tX : "<<point.X()<<"\t Y : "<<point.Y()<< endl;
	  recHit[i].fX = point.X();
	  recHit[i].fY = point.Y();
	}// for
	
	fTracker->AddRecHits(chamber,recHit,nofhitData);

      }// hit or trig condn
      

      //for(Int_t itrig = 0 ; itrig < nofTrigData ; itrig++){
      
      //}


	// Check how much space we have left and adapt this output block's size accordingly.
// 	if ( totalSize + mySize > size )
// 	    mySize = size-totalSize;
// 	Logging( kHLTLogInfo, "AliHLTMUONMansoTrackerComponent::DoEvent", "mySize set (2)", "mySize == %lu B - totalSize == %lu - size == %lu", 
// 		 mySize, totalSize, size );
// 	if ( mySize<=0 )
// 	    continue; // No room left to write a further block.
// 	// Now copy the input block
// 	unsigned long copied = 0;
// 	// First copy all full multiples of the input block
// 	while ( copied+blocks[n].fSize <= mySize )
// 	    {
// 	    Logging( kHLTLogInfo, "AliHLTMUONMansoTrackerComponent::DoEvent", "Copying", "Copying %lu B - Copied: %lu B - totalSize: %lu B", 
// 		     blocks[n].fSize, copied, totalSize );
// 	    memcpy( outputPtr+totalSize+copied, blocks[n].fPtr, blocks[n].fSize );
// 	    copied += blocks[n].fSize;
// 	    }
// 	// And the copy the remaining fragment of the block
// 	Logging( kHLTLogInfo, "AliHLTMUONMansoTrackerComponent::DoEvent", "Copying", "Copying %lu B - Copied: %lu B - totalSize: %lu B", 
// 		 mySize-copied, copied, totalSize );
// 	memcpy( outputPtr+totalSize+copied, blocks[n].fPtr, mySize-copied );
// 	Logging( kHLTLogInfo, "AliHLTMUONMansoTrackerComponent::DoEvent", "Copied", "Copied: %lu B - totalSize: %lu B", 
// 		 copied, totalSize );
// 	// Fill a block data structure for our output block.
// 	AliHLTComponentBlockData ob;
// 	// Let the structure be filled with the default values.
// 	// This takes care of setting the shared memory and data type values to default values,
// 	// so that they can be filled in by the calling code.
// 	FillBlockData( ob );
// 	// This block's start (offset) is after all other blocks written so far
// 	ob.fOffset = totalSize;
// 	// the size of this block's data.
// 	ob.fSize = mySize;
// 	// The specification of the data is copied from the input block.
// 	ob.fSpecification = blocks[n].fSpecification;
// 	// The data type is set automatically to the component's specified output data type.
// 	// Place this block into the list of output blocks
// 	outputBlocks.push_back( ob );
// 	// Increase the total amount of data written so far to our output memory
// 	totalSize += mySize;
    }
    // Finally we set the total size of output memory we consumed.
    size = totalSize;
    return 0;
}
