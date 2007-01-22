
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
#include <cstdlib>
//#include "TH2.h"

//#include "AliHLTTPCDefinitions.h"

const AliHLTComponentDataType AliHLTPHOSRawAnalyzerComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
const AliHLTComponentDataType AliHLTPHOSRawAnalyzerComponent::outputDataType=kAliHLTVoidDataType;


//AliHLTPHOSRawAnalyzerComponent gAliHLTPHOSRawAnalyzerComponent;
//ClassImp(AliHLTPHOSRawAnalyzerComponent) 
AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent():AliHLTProcessor(),  fEventCount(0),  fEquippmentId(0), fPHOSRawStream(), fRawMemoryReader(0)
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



AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & ) : AliHLTProcessor(),  fEventCount(0),  fEquippmentId(0), fPHOSRawStream(),fRawMemoryReader(0)
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
  constBase = 30;
  inputMultiplier = 0.1;
}


int AliHLTPHOSRawAnalyzerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{

  Int_t tmpMod  = 0;
  Int_t tmpRow  = 0;
  Int_t tmpCol  = 0;
  Int_t tmpGain = 0;
  Int_t processedChannels = 0;

  //  AliHLTUInt8_t *tmpOut = &((AliHLTUInt8_t)(GetEquippmentId()));

  AliHLTUInt8_t *tmpOut =( AliHLTUInt8_t*)&fEquippmentId;

  outputPtr[0] = tmpOut[0];
  outputPtr[1] = tmpOut[1];

  AliHLTUInt16_t  tmpID = *((AliHLTUInt16_t*)(outputPtr));

  cout << "OutputBuffer setting equippment ID to " << tmpID << endl;

  //  Int_t tmpMax  = 0;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;
  Reset();

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;

      if ( iter->fDataType != AliHLTPHOSDefinitions::gkDDLPackedRawDataType )
	{
	  cout << "Warning: data type = is nOT gkDDLPackedRawDataType " << endl;
	  continue;
	}

     fRawMemoryReader->SetMemory( reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize );
     fRawMemoryReader->DumpData();
     fRawMemoryReader->RewindEvents();
     analyzerPtr->SetData(fTmpChannelData);
 
      while(fPHOSRawStream->Next())
	{
	  if (fPHOSRawStream->IsNewHWAddress())
	    {
	      if(processedChannels > 0)
		{
		  analyzerPtr->SetData(fTmpChannelData);
		  analyzerPtr->Evaluate(0, 1008);
		  fMaxValues[tmpMod][tmpRow][tmpCol][tmpGain] = analyzerPtr->GetEnergy();
		  ResetDataPtr(); 
		}

	      tmpMod  =  fPHOSRawStream->GetModule();
	      tmpRow  =  fPHOSRawStream->GetRow();
	      tmpCol  =  fPHOSRawStream->GetColumn();
	      tmpGain =  fPHOSRawStream->IsLowGain(); 
	      processedChannels ++;
	    }

	  fTmpChannelData[fPHOSRawStream->GetTime()] =  fPHOSRawStream->GetSignal();
	}
    }
  DumpData();
  fEventCount++; 
  return 0;
}//end DoEvent



int
AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv )
{
  cout << "DOINIT argc =" << argc << endl;
  cout << "DOINIT argv[0] =" << argv[0] << endl;
  cout << "DOINIT argv[1] =" << argv[1] << endl;
  cout << "DOINIT argv[2] =" << argv[2] << endl;
  cout << "DOINIT argv[3] =" << argv[3] << endl;
  cout << "DOINIT argv[4] =" << argv[4] << endl;
  cout << "DOINIT argv[5] =" << argv[5] << endl;
  cout << "DOINIT argv[6] =" << argv[6] << endl;
 
  int equippmentId = atoi(argv[6]);
  cout << "The equipment ID was set to " <<equippmentId << endl;
  
 //fRawMemoryReader->SetEquipmentID(1806); 

  Reset();
  cout << "AliHLTPHOSRawAnalyzerComponent::DoInit Creating new  AliRawReaderMemory()" << endl; 
  //legoPlotPtr   = new TH2S("Lego plot 1","Phi0 20Gev, High gain", 56*5, 0, 56*5, 64, 0, 64);
  fRawMemoryReader = new AliRawReaderMemory();
  fPHOSRawStream = new  AliCaloRawStream(fRawMemoryReader,"PHOS");
  fRawMemoryReader->SetEquipmentID(equippmentId); 

  SetEquippmentId(equippmentId);

  cout <<"AliHLTPHOSRawAnalyzerComponent::DoIni  DONE!" << endl;
  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}

void
AliHLTPHOSRawAnalyzerComponent::DumpData()
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
AliHLTPHOSRawAnalyzerComponent::Reset()
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
AliHLTPHOSRawAnalyzerComponent::ResetDataPtr()
{
  for(int i = 0 ; i< 1008; i++)
    {
      fTmpChannelData[i] = 0;
    }
}


void 
AliHLTPHOSRawAnalyzerComponent::SetEquippmentId(int id)
{
  fEquippmentId = id;
}

int 
AliHLTPHOSRawAnalyzerComponent::GetEquippmentId()
{
  return  fEquippmentId;
}
