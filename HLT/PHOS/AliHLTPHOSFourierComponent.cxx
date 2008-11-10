// $Id$

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliHLTPHOSFourierComponent.h"
#include "AliHLTPHOSFourier.h"
#include "AliHLTPHOSValidCellDataStruct.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSSharedMemoryInterface.h"

#include "AliHLTPHOSRcuFFTDataStruct.h"

AliHLTPHOSFourierComponent gAliHLTPHOSFourierComponent;

AliHLTPHOSFourierComponent::AliHLTPHOSFourierComponent(): AliHLTPHOSRcuProcessor(),fFourierPtr(0), fShmPtr(0),fOutPtr(0)
{
  fFourierPtr = new AliHLTPHOSFourier();
  fShmPtr = new AliHLTPHOSSharedMemoryInterface();
}


AliHLTPHOSFourierComponent::~AliHLTPHOSFourierComponent()
{

}


int
AliHLTPHOSFourierComponent::Deinit()
{
  Logging(kHLTLogInfo, "HLT", "PHOS", ",Deinitializing AliHLTPHOSFourierComponent");
  return 0;
}


const char* 
AliHLTPHOSFourierComponent::GetComponentID()
{
  return "PhosFourier";
  //     "PhosFourier"
}


AliHLTComponent*
AliHLTPHOSFourierComponent::Spawn()

{
  return new AliHLTPHOSFourierComponent;
}


void 
AliHLTPHOSFourierComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkCellEnergyDataType);

  /*
  const AliHLTComponentDataType* pType=fgkInputDataTypes;

  while (pType->fID!=0) 
    {
      list.push_back(*pType);
      pType++;
    }
  */
}


AliHLTComponentDataType 
AliHLTPHOSFourierComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkFourierTransform;
  // return AliHLTPHOSDefinitions::fgkCellEnergyDataType;
  // return kAliHLTMultipleDataType;
}


/*
int 
AliHLTPHOSFourierComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  tgtList.clear();
  tgtList.push_back(AliHLTPHOSDefinitions::fgkFourierTransform);
  //  tgtList.push_back(kAliHLTDataTypeHwAddr16);
  return tgtList.size();
}
*/



void 
AliHLTPHOSFourierComponent::GetOutputDataSize(unsigned long& /*constBase*/, double& inputMultiplier)
{
 //  constBase =1;
//   inputMultiplier = 1;
  // constBase = 30;
  // inputMultiplier = 1.2;
  inputMultiplier = 0.1;
}
 
int 
AliHLTPHOSFourierComponent::DoInit(int /*argc*/, const char** /*argv*/)
{
  cout << "AliHLTPHOSFourierComponent::DoInit !!!!!!!!!! " << endl;
  return 0;
}

/*
int 
AliHLTPHOSFourierComponent::DoDeinit()
{

}
*/

int 
AliHLTPHOSFourierComponent::DoEvent(const AliHLTComponentEventData& evtData,
				    const AliHLTComponentBlockData* blocks, 
				    AliHLTComponentTriggerData& /*trigData*/,
				    AliHLTUInt8_t* outputPtr, 
				    AliHLTUInt32_t& size,
				    AliHLTComponentBlockDataList& outputBlocks )
{
  fPhosEventCount ++;
  AliHLTPHOSValidCellDataStruct *currentChannel =0;
  UInt_t offset            = 0; 
  UInt_t mysize            = 0;
  UInt_t tSize             = 0;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;
  UInt_t specification = 0;
 
  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
  
      if(iter->fDataType != AliHLTPHOSDefinitions::fgkCellEnergyDataType)
	{
	  
	  continue;
	}
      
      specification = specification|iter->fSpecification;	 
      AliHLTPHOSRcuCellEnergyDataStruct *cellDataPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)( iter->fPtr);
      AliHLTPHOSRcuFFTDataStruct *fourierOutPtr =  (AliHLTPHOSRcuFFTDataStruct*)outputPtr;
      fShmPtr->SetMemory(cellDataPtr);
      currentChannel = fShmPtr->NextChannel();
      
      while(currentChannel != 0)
	{
	  Int_t *data;
	  Int_t nsamples = 0;
	  data= fShmPtr->GetRawData(nsamples);
	  fFourierPtr->ProcessFourier(data, nsamples, currentChannel->fZ, currentChannel->fX, currentChannel->fGain  ); 
	  currentChannel = fShmPtr->NextChannel();
	}
          
      *fourierOutPtr =  fFourierPtr->GetPSD();
  
      for(int i=0; i < 500; i++)
	{
	  if(i%16 == 0)
	    {
	      printf("\n");
	    }

	  cout << fourierOutPtr->fGlobalAccumulatedPSD[1][i] <<  "\t";
	}
      
      mysize += sizeof(AliHLTPHOSRcuFFTDataStruct);
      AliHLTComponentBlockData bd;
      bd.fOffset = offset;
      bd.fSize = mysize;
      bd.fDataType = AliHLTPHOSDefinitions::fgkFourierTransform;
      // bd.fSpecification = 0xFFFFFFFF;
      //     bd.fSpecification = specification;
      bd.fSpecification = 1;
      outputBlocks.push_back( bd );
      



      cout <<"size left is "<<  size  <<"  FourierComponenet: offset ="<< offset << "mysize =" << mysize << "specification =" << specification <<endl;
      tSize += mysize;
   
      outputPtr += mysize;

      if( tSize > size )
	{
	  cout <<"HLT::AliHLTFourierComponent::DoEvent Too much data Data written over allowed buffer. Amount written:"<< tSize<<"allowed amount"<< size << endl;
	  Logging( kHLTLogFatal, "HLT::AliHLTFourierComponent::DoEvent", "Too much data",
		   "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu." , tSize, size );

	  return EMSGSIZE;
	}
    
    }

  return 0;

}
