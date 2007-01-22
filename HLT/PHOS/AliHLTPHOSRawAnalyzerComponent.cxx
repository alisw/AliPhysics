
/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE HLT Project.                    *
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

#include "AliHLTPHOSRawAnalyzerComponent.h"
#include <iostream>
#include "stdio.h"

#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"

//#include "AliHLTTPCDefinitions.h"

const AliHLTComponentDataType AliHLTPHOSRawAnalyzerComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
const AliHLTComponentDataType AliHLTPHOSRawAnalyzerComponent::outputDataType=kAliHLTVoidDataType;


//AliHLTPHOSRawAnalyzerComponent gAliHLTPHOSRawAnalyzerComponent;
//ClassImp(AliHLTPHOSRawAnalyzerComponent) 
AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent():AliHLTProcessor(),  eventCount(0), fPHOSRawStream(), fRawMemoryReader(0)
{
  //  fRawMemoryReader = NULL;
} 

AliHLTPHOSRawAnalyzerComponent::~AliHLTPHOSRawAnalyzerComponent()
{
  if(fRawMemoryReader != 0)
    {
      delete fRawMemoryReader;
    }
    if(fPHOSRawStream != 0)
    {
      delete fPHOSRawStream;
    }

}



AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & ) : AliHLTProcessor(),  eventCount(0), fPHOSRawStream(),fRawMemoryReader(0)
{
  //  fRawMemoryReader = NULL;
}


int 
AliHLTPHOSRawAnalyzerComponent::Deinit()
{
  return 0;
}

int 
AliHLTPHOSRawAnalyzerComponent::DoDeinit()
{
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSRawAnalyzerComponen DoDeinit");

  if(fRawMemoryReader !=0)
    {
      delete fRawMemoryReader;
    }
    
  if(fPHOSRawStream != 0)
    {
      delete fPHOSRawStream;
    }
  return 0;

}

const char* 
AliHLTPHOSRawAnalyzerComponent::GetComponentID()
{
  return "AliPhosTestRaw";
}

void
AliHLTPHOSRawAnalyzerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=inputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSRawAnalyzerComponent::GetOutputDataType()
{
  return outputDataType;
}

void
AliHLTPHOSRawAnalyzerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  constBase = 0;inputMultiplier = 0;
}


int AliHLTPHOSRawAnalyzerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  //  int tmpCnt = 0;
  cout << "processing Event" << endl;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;
    
  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      cout <<"Inside for loop ndx =" << ndx << endl;
      iter = blocks+ndx;
      
       if (eventCount == 0)
	{
	  continue;
	}
 
      
      if ( iter->fDataType == AliHLTPHOSDefinitions::gkDDLPackedRawDataType ) cout << "data type is :gkDDLPackedRawDataType " <<endl;
      if ( iter->fDataType == AliHLTPHOSDefinitions::gkPackedRawDataType ) cout << "data type is : gkPackedRawDataType" <<endl;
      if ( iter->fDataType == AliHLTPHOSDefinitions::gkUnpackedRawDataType) cout << "data type is gkUnpackedRawDataType" <<endl;
      if ( iter->fDataType == AliHLTPHOSDefinitions::gkClustersDataType) cout << "data type is ::gkClustersDataType" <<endl;
      if ( iter->fDataType == AliHLTPHOSDefinitions::gkVertexDataType ) cout << "data type is ::gkVertexDataType " <<endl;
      if ( iter->fDataType == AliHLTPHOSDefinitions::gkTrackSegmentsDataType) cout << "data type is :::gkTrackSegmentsDataType" <<endl;

      if ( iter->fDataType != AliHLTPHOSDefinitions::gkDDLPackedRawDataType )
	{
	  cout << "Warning: data type = is nOT gkDDLPackedRawDataType " << endl;
	  continue;
	}

     cout <<"PHOSHLT DoEvent: processing event:"<< eventCount << endl;
     cout <<"Struct size = " << evtData.fStructSize << endl;
     cout <<"Event ID = " <<evtData.fEventID << endl;
     cout <<"Block count = " <<evtData.fBlockCnt << endl;
     cout <<"Block size = " << blocks->fSize << endl;
     cout <<"printing out start od data block" << endl;
     //     cout << "content of data pointer =" << tmpDtaPtr << endl;  

     //     UChar_t *tmpDtaPtr  = (UChar_t *)blocks->fPtr;
     UInt_t *tmpDtaPtr  = (UInt_t *)blocks->fPtr;

    
     for(unsigned int i = 0; i < 15; i++)
       {
	 // getc();
	 //     sleep(10);
	 //    printf("\ntype return to continue; ");
	 //     getc(stdin);
	 //    printf("\nThanks; read in \n");
	 cout << "entry:" <<i <<" =  " << *tmpDtaPtr << endl;
	 tmpDtaPtr ++;
	 //	 cout << "entry:" <<i <<" =  " << *blocks << endl;
	 //	 blocks ++;
       } 


     fRawMemoryReader->SetMemory( reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize );
     bool readValue = true;
     
 
     cout << endl << "eventCount =" << eventCount << endl;
    fRawMemoryReader->DumpData();

     readValue = fPHOSRawStream->Next();
     //    fPHOSRawStream->Next();  

     while( readValue )
       {
	 cout <<"reading value" << endl;
	 readValue = fPHOSRawStream->Next();
       }
     cout << "end of for lop" << endl;
   }

  eventCount++; 

 return 0;
}



int
AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv )
{
  cout << "AliHLTPHOSRawAnalyzerComponent::DoInit Creating new  AliRawReaderMemory()" << endl; 
  fRawMemoryReader = new AliRawReaderMemory();
  fPHOSRawStream = new  AliCaloRawStream(fRawMemoryReader,"PHOS");
  cout <<"AliHLTPHOSRawAnalyzerComponent::DoIni  DONE!" << endl;
  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}



/*
int 
AliHLTPHOSRawAnalyzerComponent::DoEvent(const AliHLTComponentEventData& evtDtaPtr, const AliHLTComponentBlockData* dtaPtr, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&)
{
  const Bool_t skipBadEvent=kFALSE;
  const AliHLTComponentBlockData* iter = NULL;

  fRawMemoryReader->SetMemory( reinterpret_cast<UChar_t*>( dtaPtr->fPtr ), dtaPtr->fSize );

  //  fPHOSRawStream = new  AliCaloRawStream(fRawMemoryReader,"PHOS");

  fRawMemoryReader->DumpData();

  Logging(kHLTLogInfo, "HLT", "Sample", "PhosHLTRawAnalyzerComonent, DoEvent");

  //  UChar_t *tmpDtaPtr  = (UChar_t *)dtaPtr->fPtr;
  //  UInt32_t *tmpDtaPtr  = (UInt32_t *)dtaPtr->fPtr;
  //  unsigned long *tmpDtaPtr  = (unsigned long *)dtaPtr->fPtr;

  
  cout <<"PHOSHLT DoEvent: processing event:"<< eventCount << endl;
  cout <<"Struct size = " <<evtDtaPtr.fStructSize << endl;
  cout <<"Event ID = " <<evtDtaPtr.fEventID << endl;
  cout <<"Block count = " <<evtDtaPtr.fBlockCnt << endl;
  cout <<"Block size = " << dtaPtr->fSize << endl;
  cout <<"printing out start od data block" << endl;
  cout << "content of data pointer =" << tmpDtaPtr << endl;
  

  // UChar_t *tmpDtaPtr  = (UChar_t *)blokcs->fPtr;

  //  AliCaloRawStream in(fRawMemoryReader,"PHOS");
  //  for(unsigned int i = 0; i < tmpSize; i++)
  for(unsigned int i = 0; i < 15; i++)
    {
      // getc();
      //     sleep(10);
      //    printf("\ntype return to continue; ");
      //     getc(stdin);
      //    printf("\nThanks; read in \n");
      //    cout << "entry:" <<i <<" =  " << *tmpDtaPtr << endl;
      //   tmpDtaPtr ++;
   }

  
  
  cout << "Entering while loop" << endl;
  while ( fPHOSRawStream->Next() ) 
    { 
      cout << "Inside while loop" << endl;
      if (fPHOSRawStream->IsNewHWAddress()) cout << endl << " New HW address: " << endl;
      cout << "time bin=" << fPHOSRawStream->GetTime() << " of "<< fPHOSRawStream->GetTimeLength()
	   << ", amp=" <<  fPHOSRawStream->GetSignal() <<" gain="<< fPHOSRawStream->IsLowGain()
	   << " HW address="<<fPHOSRawStream->GetHWAddress()
	   << " channel=("<<fPHOSRawStream->GetModule() <<","<<fPHOSRawStream->GetColumn()<<","<<fPHOSRawStream->GetRow()<<")"<< endl;
      
      
      if (skipBadEvent && (fPHOSRawStream->GetTime() < 0 || fPHOSRawStream->GetTimeLength() <=0) ) {
	cout << "Wrong time bin or time length. Skip this event" << endl;
	break;
      }
    }
  cout << "Finnsihed while loop" << endl << endl;

  eventCount++;
  return 0;
}
*/
