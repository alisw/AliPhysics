
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

const AliHLTComponentDataType AliHLTPHOSRawAnalyzerComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
const AliHLTComponentDataType AliHLTPHOSRawAnalyzerComponent::outputDataType=kAliHLTVoidDataType;


//AliHLTPHOSRawAnalyzerComponent gAliHLTPHOSRawAnalyzerComponent;
//ClassImp(AliHLTPHOSRawAnalyzerComponent) 
AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent():AliHLTProcessor(), eventCount(0)
{
  //  ddlDecoderPtr = new AliRawReaderMemory();
  //  fRawMemoryReader = NULL;
} 

AliHLTPHOSRawAnalyzerComponent::~AliHLTPHOSRawAnalyzerComponent()
{
  // fRawMemoryReader = NULL;
}


AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & ) : AliHLTProcessor(), eventCount(0)
{
  //  ddlDecoderPtr = new AliRawReaderMemory();
  //  fRawMemoryReader = NULL;
}


//AliHLTComponent*
//AliHLTPHOSRawAnalyzerComponent::Spawn()
//{
//  return new AliHLTPHOSRawAnalyzerComponent;
//}



int 
AliHLTPHOSRawAnalyzerComponent::Deinit()
{
  return 0;
}

int 
AliHLTPHOSRawAnalyzerComponent::DoDeinit()
{
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSRawAnalyzerComponen DoDeinit");

  //  if ( fRawMemoryReader )
  //   delete fRawMemoryReader;
  //  fRawMemoryReader = NULL;
 
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

int 
AliHLTPHOSRawAnalyzerComponent::DoEvent(const AliHLTComponentEventData& evtDtaPtr, const AliHLTComponentBlockData* dtaPtr, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&)
{
 
  //  AliRawReaderMemory *testDDLDecoderPtr = new AliRawReaderMemory() ;
 //  virtual Bool_t   SetMemory( UChar_t* memory, ULong_t size );

  //  ddlDecoderPtr->SetMemory( reinterpret_cast<UChar_t*>( dtaPtr->fPtr ), dtaPtr->fSize );

  //  analyzerPtr->Evaluate(0, 100);

  Logging(kHLTLogInfo, "HLT", "Sample", "PhosHLTRawAnalyzerComonent, DoEvent");

  //  AliHLTUInt32_t *tmpDtaPtr  = (AliHLTUInt32_t *)dtaPtr->fPtr;


  //  UChar_t *dataPtr = address;
  //ULong_t address,

  UChar_t *tmpDtaPtr  = (UChar_t *)dtaPtr->fPtr;
  // UChar_t *tmpDtaPtr  = (ULong_t *)dtaPtr->fPtr;


  AliHLTUInt32_t tmpSize =  dtaPtr->fSize;
  AliHLTUInt32_t tmpSize32 = dtaPtr->fSize/4;

  cout <<"PHOSHLT DoEvent: processing event:"<< eventCount << endl;
  cout <<"Struct size = " <<evtDtaPtr.fStructSize << endl;
  cout <<"Event ID = " <<evtDtaPtr.fEventID << endl;
  cout <<"Block count = " <<evtDtaPtr.fBlockCnt << endl;
  cout <<"Block size = " << dtaPtr->fSize << endl;
  cout <<"printing out start od data block" << endl;

  
  //  tmpDtaPtr = dtaPtr->fPtr;
  //  AliHLTUInt32_t *tmpDtaPtr  = (AliHLTUInt32_t *)dtaPtr->fPtr;

  
  //  testDDLDecoderPtr->SetMemory(tmpDtaPtr, tmpSize);
  //  cout << "Dumoing data from the AliRawReader" << endl;
  // ddlDecoderPtr->DumpData(30);
  // cout << "finnished dumping data from the AliRawReader" << endl;
  
  //  unsigned int tmpSize =   (dtaPtr->fSize)/4;
 
  cout << "content of data pointer =" << tmpDtaPtr << endl;
  

  //  for(unsigned int i = 0; i < tmpSize; i++)
  for(unsigned int i = 0; i < 10; i++)
    {
      // getc();
      //     sleep(10);
      printf("\ntype return to continue; ");
      getc(stdin);
      printf("\nThanks; read in \n");
      cout << "entry:" <<i <<" =  " << *tmpDtaPtr << endl;
      tmpDtaPtr ++;
   }
 

  eventCount++;
  return 0;
}

//void
//AliHLTPHOSRawAnalyzerComponent::Evaluate(int start, int length)
//{
//  cout <<" AliHLTPHOSRawAnalyzerComponent::Evaluate" << endl ;
//}


int
AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv )
{
  cout << "AliHLTPHOSRawAnalyzerComponent::DoInit Creating new  AliRawReaderMemory()" << endl; 
  fRawMemoryReader = new AliRawReaderMemory;



  /*
  if(fRawMemoryReader)
    {
      return EINPROGRESS;
    }
  int i = 0;
  //     const char* tableFileBaseDir = NULL;
  while ( i < argc )
    {
      Logging(kHLTLogError, "HLT::TPCRawDataUnpacker::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
      return EINVAL;
    }
  
  //  fRawMemoryReader = new AliRawReaderMemory;
  */





  cout <<"AliHLTPHOSRawAnalyzerComponent::DoIni  DONE!" << endl;
  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}
