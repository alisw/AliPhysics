
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Per Thomas Hille, Oystein Djuvsland                   *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSRawAnalyzerComponentv3.h"
#include "AliHLTPHOSRawAnalyzer.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSUtilities.h"
#include "AliHLTPHOSMapper.h"
#include "AliLog.h"


AliHLTPHOSRawAnalyzerComponentv3::AliHLTPHOSRawAnalyzerComponentv3() :
   AliHLTCaloRawAnalyzerComponentv3("PHOS")
   ,fCurrentSpec(-1)
{
   // See header file for class documentation
   InitMapping(0x1); //using 0x1 to avoid error message
}

AliHLTPHOSRawAnalyzerComponentv3::~AliHLTPHOSRawAnalyzerComponentv3()
{
  //comment
}


void
AliHLTPHOSRawAnalyzerComponentv3::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //comment
  list.clear();
  list.push_back( AliHLTPHOSDefinitions::fgkDDLPackedRawDataType);
}

AliHLTComponentDataType 
AliHLTPHOSRawAnalyzerComponentv3::GetOutputDataType()
{
  //comment
  return AliHLTPHOSDefinitions::fgkChannelDataType;
}

void
AliHLTPHOSRawAnalyzerComponentv3::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  //comment
  constBase = sizeof(AliHLTCaloChannelDataHeaderStruct);
  inputMultiplier = 1.5;
}

int 
AliHLTPHOSRawAnalyzerComponentv3::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& /*trigData*/, 
					 AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  //comment


   if(!IsDataEvent())
   {
      size = 0;
      return 0;
   }


  Int_t blockSize          = 0;
  UInt_t totSize           = 0;

  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if ( iter->fDataType != AliHLTPHOSDefinitions::fgkDDLPackedRawDataType  )
	{
	  HLTDebug("Data block is not of type fgkDDLPackedRawDataType");
	  continue; 
	}
      if(iter->fSpecification != fCurrentSpec)
      {
	 fCurrentSpec = iter->fSpecification;
	 InitMapping(iter->fSpecification);
      }
      blockSize = DoIt(iter, outputPtr, size, totSize); // Processing the block

      if(blockSize == -1) // If the processing returns -1 we are out of buffer and return an error msg.
	{
	  return -ENOBUFS;
	}

      totSize += blockSize; //Keeping track of the used size
      // HLTDebug("Output data size: %d - Input data size: %d", totSize, iter->fSize);

      //Filling the block data
      AliHLTComponentBlockData bdChannelData;
      FillBlockData( bdChannelData );
      bdChannelData.fOffset = 0; //FIXME
      bdChannelData.fSize = blockSize;
      bdChannelData.fDataType = AliHLTPHOSDefinitions::fgkChannelDataType;
      bdChannelData.fSpecification = iter->fSpecification;
      outputBlocks.push_back(bdChannelData);

      outputPtr += blockSize; //Updating position of the output buffer
    }

  fCaloEventCount++; 
  size = totSize; //telling the framework how much buffer space we have used.
  
  return 0;
}//end DoEvent


void AliHLTPHOSRawAnalyzerComponentv3::InitMapping ( const int specification )
{
   // See header file for class documentation
   fMapperPtr = new AliHLTPHOSMapper;
   fMapperPtr->InitDDLSpecificationMapping();
   fMapperPtr->InitAltroMapping(specification);

}

