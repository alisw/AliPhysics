/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Boris Polichtchouk & Per Thomas Hille for the ALICE           *
 * offline/HLT Project. Contributors are mentioned in the code where      *
 * appropriate.                                                           *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <iostream>
#include "stdio.h"
#include <cstdlib>
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
//#include "AliHLTPHOSModuleCellAccumulatedEnergyDataStruct.h"
#include "AliHLTPHOSRcuHistogramProducer.h"
#include "AliHLTPHOSRcuHistogramProducerComponent.h"



const AliHLTComponentDataType  AliHLTPHOSRcuHistogramProducerComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
const AliHLTComponentDataType  AliHLTPHOSRcuHistogramProducerComponent::outputDataType=kAliHLTVoidDataType;
AliHLTPHOSRcuHistogramProducerComponent gAliHLTPHOSRcuHistogramProducerComponent;

//AliHLTPHOSHistogramProducerComponent gAliHLTPHOSHistogramProducerComponent;
/*************************************************************************
* Class AliHLTPHOSRcuHistogramProducerComponent accumulating histograms  *
* with amplitudes per PHOS channel                                       *
* It is intended to run at the HLT farm                                  *
* and it fills the histograms with amplitudes per channel.               * 
* Usage example see in PHOS/macros/Shuttle/AliPHOSCalibHistoProducer.C   *
**************************************************************************/
AliHLTPHOSRcuHistogramProducerComponent:: AliHLTPHOSRcuHistogramProducerComponent():AliHLTProcessor(), fEventCount(0), fRcuHistoProducerPtr(0)
{
  //  Reset();
} 



AliHLTPHOSRcuHistogramProducerComponent::~ AliHLTPHOSRcuHistogramProducerComponent()
{

}


AliHLTPHOSRcuHistogramProducerComponent::AliHLTPHOSRcuHistogramProducerComponent(const  AliHLTPHOSRcuHistogramProducerComponent & ) : AliHLTProcessor(), fEventCount(0), fRcuHistoProducerPtr(0)
{

}


int 
AliHLTPHOSRcuHistogramProducerComponent::Deinit()
{
  cout << "AliHLTPHOSRcuHistogramProducerComponent::Deinit()" << endl;
  fRcuHistoProducerPtr->WriteEnergyHistograms();
  return 0;
}


int 
AliHLTPHOSRcuHistogramProducerComponent::DoDeinit()
{
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSRcuHistogramProducer DoDeinit");
  return 0;
}


const char* 
AliHLTPHOSRcuHistogramProducerComponent::GetComponentID()
{
  return "RcuHistogramProducer";
}


void
 AliHLTPHOSRcuHistogramProducerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=inputDataTypes;
  while (pType->fID!=0) 
    {
      list.push_back(*pType);
      pType++;
    }
}


AliHLTComponentDataType 
AliHLTPHOSRcuHistogramProducerComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::gkCellEnergyDataType;
}


void
AliHLTPHOSRcuHistogramProducerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  constBase = 30;
  inputMultiplier = 1;
}


int  AliHLTPHOSRcuHistogramProducerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  unsigned long ndx       = 0;
  UInt_t offset           = 0; 
  UInt_t mysize           = 0;
  UInt_t tSize            = 0;
  const AliHLTComponentBlockData* iter = NULL;   
  AliHLTPHOSRcuCellEnergyDataStruct *cellDataPtr;
  AliHLTUInt8_t* outBPtr;
 
  // outBPtr = outputPtr;
  // fOutPtr =  (AliHLTPHOSRcuCellAccumulatedEnergyDataStruct*)outBPtr;
  int tmpCnt;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      cellDataPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)( iter->fPtr);
      tmpCnt =  cellDataPtr->fCnt;

      for(int i= 0; i <= tmpCnt; i ++)
	{
	  fRcuHistoProducerPtr->FillEnergy(cellDataPtr->fValidData[i].fX,
					   cellDataPtr->fValidData[i].fZ, 
					   cellDataPtr->fValidData[i].fGain, 
					   cellDataPtr->fValidData[i].fEnergy);
	}
    }
  
  outBPtr = outputPtr;
  fOutPtr =  (AliHLTPHOSRcuCellAccumulatedEnergyDataStruct*)outBPtr;
  const AliHLTPHOSRcuCellAccumulatedEnergyDataStruct  &innPtr = fRcuHistoProducerPtr->GetCellAccumulatedEnergies();

  fOutPtr->fModuleID = fModuleID;
  fOutPtr->fRcuX     = fRcuX;
  fOutPtr->fRcuZ     = fRcuZ;


  for(int x=0; x < N_XCOLUMNS_RCU; x ++)
    {
      for(int z=0; z < N_XCOLUMNS_RCU; z ++)
	{
	  for(int gain =0;  gain < N_GAINS; gain ++)
	    {
	      fOutPtr->fAccumulatedEnergies[x][z][gain] = innPtr.fAccumulatedEnergies[x][z][gain];
	      fOutPtr->fHits[x][z][gain] = innPtr.fHits[x][z][gain];
	    }
	}
    }


  //pushing data to shared output memory
  mysize += sizeof(AliHLTPHOSRcuCellAccumulatedEnergyDataStruct);
  AliHLTComponentBlockData bd;
  FillBlockData( bd );
  bd.fOffset = offset;
  bd.fSize = mysize;
  bd.fDataType = AliHLTPHOSDefinitions::gkCellAccumulatedEnergyDataType;
  bd.fSpecification = 0xFFFFFFFF;
  outputBlocks.push_back( bd );
  tSize += mysize;
  outBPtr += mysize;

  if( tSize > size )
    {
      Logging( kHLTLogFatal, "HLT::AliHLTRcuHistogramProducerComponent::DoEvent", "Too much data",
	       "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
	       , tSize, size );
      return EMSGSIZE;
    }

  fEventCount++; 
  return 0;
}//end DoEvent


int
AliHLTPHOSRcuHistogramProducerComponent::DoInit( int argc, const char** argv )
{
  int iResult=0;
  TString argument="";
  //  fRcuHistoProducerPtr = new AliHLTPHOSRcuHistogramProducer();
  AliHLTUInt8_t tmpRcuX;
  AliHLTUInt8_t tmpRcuZ; 
  AliHLTUInt8_t tmpModuleID;
  Bool_t isSetEquippmentID = kFALSE;

  for(int i=0; i<argc && iResult>=0; i++) 
    {
      argument=argv[i];
      
      if (argument.IsNull()) 
	{
	  continue;
	}
      if (argument.CompareTo("-equipmentID") == 0) 
	{
	  if(i+1 <= argc)
	    {
	      fEquippmentID = atoi(argv[i+1]);
	      isSetEquippmentID = kTRUE;
	    }
	  else
	    {
	       iResult= -1;
	       Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuHistogramProducerComponent::DoInt( int argc, const char** argv )", "Missing argument",
			"The argument -equippmentID expects a number");
	       return  iResult;   
	    }
	}

      int rcuIndex =  (fEquippmentID - 1792)%N_RCUS_PER_MODULE;
      //     fModuleID = (fEquippmentID  -1792 -rcuIndex)/N_RCUS_PER_MODULE;
      tmpModuleID = ((fEquippmentID -1792 -rcuIndex)/N_RCUS_PER_MODULE);
      SetModuleID(tmpModuleID);

      if(rcuIndex == 0)
	{
	  tmpRcuX = 0; 
	  tmpRcuZ = 0;
	}
      
      if(rcuIndex == 1)
	{
	 tmpRcuX = 0; 
	 tmpRcuZ = 1;
	}
      
      if(rcuIndex == 2)
	{
	  tmpRcuX = 1; 
	  tmpRcuZ = 0;
	}
      
      if(rcuIndex == 3)
	{
	  tmpRcuX = 1; 
	  tmpRcuZ = 1;
	}

      SetRcuX(tmpRcuX);
      SetRcuZ(tmpRcuZ); 
      cout <<"********InitInfo************"<< endl;
      cout <<"AliHLTPHOSRcuHistogramProducerComponent::SetCoordinate"<< endl;
      cout <<"Equpippment ID =\t"<< fEquippmentID <<endl;
      cout <<"Module ID =\t"<<  (int) tmpModuleID<<endl;
      cout <<"RCUX =\t\t" << (int)tmpRcuX << endl;
      cout <<"RCUZ =\t\t" << (int)tmpRcuZ << endl;
     }

  if(isSetEquippmentID == kFALSE)
    {
       Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuHistogramProducerComponent::DoInt( int argc, const char** argv )", "Missing argument",
			"The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>");
       iResult = -2; 
    }

  
  fRcuHistoProducerPtr = new AliHLTPHOSRcuHistogramProducer( tmpModuleID, tmpRcuX, tmpRcuZ);


  return  iResult;
}


void 
AliHLTPHOSRcuHistogramProducerComponent::SetRcuX(AliHLTUInt8_t X)
{
  fRcuX = X;
}


void 
AliHLTPHOSRcuHistogramProducerComponent::SetRcuZ(AliHLTUInt8_t Z)
{
  fRcuZ = Z;
}


void 
AliHLTPHOSRcuHistogramProducerComponent::SetModuleID(AliHLTUInt8_t moduleID)
{
  fModuleID = moduleID;
}


void 
AliHLTPHOSRcuHistogramProducerComponent::SetEquippmentId(int id)
{
  fEquippmentID = id;
  fRcuHistoProducerPtr->SetEquippmentId(id);
}


int 
AliHLTPHOSRcuHistogramProducerComponent::GetEquippmentId()
{
  return  fEquippmentID;
}


AliHLTComponent*
AliHLTPHOSRcuHistogramProducerComponent::Spawn()
{
  return new AliHLTPHOSRcuHistogramProducerComponent;
}


