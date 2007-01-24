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

//#include "AliHLTPHOSModuleMergerComponent.h"

#include "AliHLTPHOSModuleMergerComponent.h"
#include <iostream>
#include "stdio.h"

#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include <cstdlib>
//#include "TH2.h"

//#include "AliHLTTPCDefinitions.h"

const AliHLTComponentDataType  AliHLTPHOSModuleMergerComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
const AliHLTComponentDataType  AliHLTPHOSModuleMergerComponent::outputDataType=kAliHLTVoidDataType;


AliHLTPHOSModuleMergerComponent gAliHLTPHOSModuleMergerComponent;
//ClassImp( AliHLTPHOSModuleMergerComponent) 
AliHLTPHOSModuleMergerComponent:: AliHLTPHOSModuleMergerComponent():AliHLTProcessor(),  fEventCount(0),  fEquippmentID(0)
{

} 

AliHLTPHOSModuleMergerComponent::~ AliHLTPHOSModuleMergerComponent()
{

}

AliHLTPHOSModuleMergerComponent:: AliHLTPHOSModuleMergerComponent(const  AliHLTPHOSModuleMergerComponent & ) : AliHLTProcessor(),  fEventCount(0),  fEquippmentID(0)
{

}


int 
AliHLTPHOSModuleMergerComponent::Deinit()
{
  return 0;
}

int 
AliHLTPHOSModuleMergerComponent::DoDeinit()
{
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSModuleMerger DoDeinit");
  return 0;

}

const char* 
AliHLTPHOSModuleMergerComponent::GetComponentID()
{
  return "ModuleMerger";
}

void
 AliHLTPHOSModuleMergerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=inputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSModuleMergerComponent::GetOutputDataType()
{
  // return AliHLTPHOSDefinitions::gkUnpackedRawDataType;
  return AliHLTPHOSDefinitions::gkCellEnergyDataType;
  //  return outputDataType;

 }

void
AliHLTPHOSModuleMergerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  constBase = 30;
  inputMultiplier = 1;
}


int  AliHLTPHOSModuleMergerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  unsigned long ndx;
  const AliHLTComponentBlockData* iter = NULL;   
 
  // cout << "Inside AliHLTPHOSModuleMergerComponent::DoEvent" << endl;
  //  AliHLTUInt32_t *tmp = 
  //  cout << " AliHLTPHOSModuleMergerComponen: the size of ouputblocks is " << outputBlocks.size()  <<endl;
  //  cout << " AliHLTPHOSModuleMergerComponen: the size of inputblock is "  << blocks->fSize  <<endl;
  //  cout << " AliHLTPHOSModuleMergerComponen: the size "  << size  <<endl;
  
  //  cout << "AliHLTPHOSModuleMergerComponen: evtData  fStructSize  =  "<< evtData.fStructSize<<endl;
  //  cout << "AliHLTPHOSModuleMergerComponen: evtData  EventID      =  "<< evtData.fEventID<<endl;
  //  cout << "AliHLTPHOSModuleMergerComponen: evtData  block count  =  "<< evtData.fBlockCnt<<endl;
  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      //   iter = ndx; 

      AliHLTUInt32_t *tmpPtr = reinterpret_cast<AliHLTUInt32_t*>( iter->fPtr); 
 
      //
      //  if( iter->fDataType != AliHLTPHOSDefinitions::gkDDLPackedRawDataType )
      //  	{
      // 	  cout << "Warning: data type = is nOT gkDDLPackedRawDataType " << endl;
      //	  continue;
      //	}


      if (iter->fDataType == kAliHLTVoidDataType){ cout << "ModuleMerger: datatype is kAliHLTVoidDataType :" << endl;} 
      if (iter->fDataType ==  kAliHLTAnyDataType){ cout << "ModuleMerger: datatype is kAliHLTAnyDataType :" << endl;} 
      if (iter->fDataType == AliHLTPHOSDefinitions::gkDDLPackedRawDataType){ cout << "ModuleMerger: datatType is : AliHLTPHOSDefinitions::gkDDLPackedRawDataType :" << endl;} 
      if (iter->fDataType == AliHLTPHOSDefinitions::gkCellEnergyDataType){ cout << "ModuleMerger: datatype isAliHLTPHOSDefinitions::gkCellEnergyDataType" << endl;} 
   
      //  if (iter->fDataType == AliHLTPHOSDefinitions::gkPackedRawDataType){ cout << "ModuleMerger: datatype is: AliHLTPHOSDefinitions::gkPackedRawDataType" << endl;} 
      //  if (iter->fDataType == AliHLTPHOSDefinitions::gkUnpackedRawDataType){ cout << "ModuleMerger: datatype is: AliHLTPHOSDefinitions::gkUnpackedRawDataType" << endl;} 
      //  if (iter->fDataType == AliHLTPHOSDefinitions::gkClustersDataType){ cout << "ModuleMerger: datatype is: AliHLTPHOSDefinitions::gkClustersDataType" << endl;} 
      //  if (iter->fDataType == AliHLTPHOSDefinitions::gkVertexDataType){ cout << "ModuleMerger: datatype is: AliHLTPHOSDefinitions::gkVertexDataTypeg" << endl;} 
      //  if (iter->fDataType == AliHLTPHOSDefinitions::gkTrackSegmentsDataType){ cout << "ModuleMerger: datatype is: AliHLTPHOSDefinitions::gkTrackSegmentsDataTypeg" << endl;} 

    }



  cout << "blocks.fSize = " << blocks->fSize << endl ; 

  //  AliHLTUInt32_t* tmpPtr = reinterpret_cast<AliHLTUInt32_t*>(blocks->fPtr);

  AliHLTUInt32_t* tmpPtr = (AliHLTUInt32_t*)(blocks->fPtr);
  //  *tmpPtr = 100;



  cout << "ModuleMerge*tmpPtr =" << *tmpPtr << endl;
  cout << "ModuleMerge tmpPtr =" << tmpPtr << endl;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      //     AliHLTUInt16_t* tmpPtr = reinterpret_cast<AliHLTUInt16_t*>(blocks->fPtr);   
      iter = blocks+ndx;
      //      AliHLTUInt16_t *tmp =  (AliHLTUInt16_t *)(iter->fPtr);
      //    cout <<"AliHLTPHOSModuleMergerComponent::Equippment ID *tmp= " <<  *tmp  << endl;
      //      cout <<"AliHLTPHOSModuleMergerComponent::Equippment ID  blocks->fPtr[0] = " << tmpPtr[0] << endl;
    }
  
  
  fEventCount++; 
  return 0;
}//end DoEvent



int
AliHLTPHOSModuleMergerComponent::DoInit( int argc, const char** argv )
{
  Reset();

  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}

void
AliHLTPHOSModuleMergerComponent::DumpData()
{
  for(int mod = 0; mod <5; mod ++)
    {
      printf("\n ***********  MODULE %d ************\n", mod);
      for(int row = 0; row < 64; row ++)
	{
	  for(int col = 0; col < 56; col ++)
	    {
	      if( fMaxValues[mod][row][col][0] != 0)
		{ 
		  cout << fMaxValues[mod][row][col][0] << "\t";
		}
	    }
	} 
    }
}


void
AliHLTPHOSModuleMergerComponent::Reset()
{
  for(int mod = 0; mod <5; mod ++)
    {
      for(int row = 0; row < 64; row ++)
	{
	  for(int col = 0; col < 56; col ++)
	    {
	      for(int gain = 0; gain <2; gain ++ )
		{
		  fMaxValues[mod][row][col][gain] = 0;
		}
	    }
	}
    }

  for(int i = 0 ; i< 1008; i++)
    {
      fTmpChannelData[i] = 0;
    }
} // end Reset

void
AliHLTPHOSModuleMergerComponent::ResetDataPtr()
{
  for(int i = 0 ; i< 1008; i++)
    {
      fTmpChannelData[i] = 0;
    }
}


void 
AliHLTPHOSModuleMergerComponent::SetEquippmentId(int id)
{
  fEquippmentID = id;
}

int 
AliHLTPHOSModuleMergerComponent::GetEquippmentId()
{
  return  fEquippmentID;
}


AliHLTComponent*
AliHLTPHOSModuleMergerComponent::Spawn()
{
  return new AliHLTPHOSModuleMergerComponent;
}


