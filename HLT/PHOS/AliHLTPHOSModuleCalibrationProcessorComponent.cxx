// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliPHOSRcuDA1.h"
#include "AliHLTPHOSSharedMemoryInterface.h"
#include "AliHLTPHOSModuleCalibrationProcessorComponent.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"

//#include <iostream>

/** @file   AliHLTPHOSModuleCalibrationProcessorComponent.cxx
    @author Oystein Djuvsland
    @date   
    @brief  A module calibration component for PHOS HLT, using the PHOS DA's
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

using namespace PhosHLTConst;

AliHLTPHOSModuleCalibrationProcessorComponent gAliHLTPHOSModuleCalibrationProcessorComponent;

AliHLTPHOSModuleCalibrationProcessorComponent::AliHLTPHOSModuleCalibrationProcessorComponent() :
  AliHLTCalibrationProcessor(),
  fCnt(0),
  fPHOSDAPtr(0),
  fShmPtr(0) 
{
  //See header file for documentation
}


AliHLTPHOSModuleCalibrationProcessorComponent::~AliHLTPHOSModuleCalibrationProcessorComponent() 
{
  //See header file for documentation
  if(fShmPtr)
    {
      delete fShmPtr;
      fShmPtr = 0;
    }
  if(fPHOSDAPtr)
    {
      delete fPHOSDAPtr;
      fPHOSDAPtr = 0;
    }
}

void AliHLTPHOSModuleCalibrationProcessorComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{
  //See header file for documentation
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkCellEnergyDataType);
}


AliHLTComponentDataType AliHLTPHOSModuleCalibrationProcessorComponent::GetOutputDataType()
{
  //See header file for documentation
  return AliHLTPHOSDefinitions::fgkEmcCalibDataType;
}
                                     
void AliHLTPHOSModuleCalibrationProcessorComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //See header file for documentation
  constBase = 0;
  inputMultiplier = 2;
}

AliHLTComponent* 
AliHLTPHOSModuleCalibrationProcessorComponent::Spawn()
{
  //See header file for documentation
  return new AliHLTPHOSModuleCalibrationProcessorComponent();
}

const char* 
AliHLTPHOSModuleCalibrationProcessorComponent::GetComponentID()
{
  //See header file for documentation
  return "PhosModuleCalibrationProcessor";  
}

Int_t 
AliHLTPHOSModuleCalibrationProcessorComponent::ScanArgument( Int_t /*argc*/, const char** /*argv*/)
{
  //See header file for documentation
  return 0;
}

Int_t AliHLTPHOSModuleCalibrationProcessorComponent::InitCalibration()
{  
  //See header file for documentation
  fShmPtr = new AliHLTPHOSSharedMemoryInterface();
  return 0;
}


Int_t AliHLTPHOSModuleCalibrationProcessorComponent::DeinitCalibration()
{
   //See header file for documentation

  cout << "Int_t AliHLTPHOSModuleCalibrationProcessorComponent::DeinitCalibration()"  << endl;
  //AliHLTPHOSModuleCalibrationProcessorComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) {
  
  AliHLTComponentEventData dummyEvtData;
  AliHLTComponentTriggerData dummyTrgData;

  ShipDataToFXS(dummyEvtData, dummyTrgData); 

  //See header file for documentation
  if(fShmPtr)
    {
      delete fShmPtr;
      fShmPtr = 0;
    }
  if(fPHOSDAPtr)
    {
      delete fPHOSDAPtr;
      fPHOSDAPtr = 0;
    }
  return 0;
}

Int_t AliHLTPHOSModuleCalibrationProcessorComponent::ProcessCalibration(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData)
{
  
  //  fCnt ++;

  fCnt ++;

  //See header file for documentation
  const  AliHLTComponentEventData eDta  = evtData;
  AliHLTComponentTriggerData  tDta =  trigData;


  AliHLTPHOSValidCellDataStruct *currentChannel =0;

  UInt_t specification = 0;
  const AliHLTComponentBlockData* iter = 0;
  iter = GetFirstInputBlock( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC);
  AliHLTPHOSRcuCellEnergyDataStruct* cellDataPtr = 0;
  Int_t xOffset = 0;
  Int_t zOffset = 0;
  Int_t module = -1;

  Float_t energyArray[NXCOLUMNSMOD][NZROWSMOD][NGAINS];
  Float_t timeArray[NXCOLUMNSMOD][NZROWSMOD][NGAINS];

  while(iter != 0)
    {

      specification = specification|iter->fSpecification;

      cellDataPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)( iter->fPtr);

      // Only necessary until the specification gives the module number
      module = cellDataPtr->fModuleID;
      if(!fPHOSDAPtr)
	{
	  //	  fPHOSDAPtr = new AliPHOSDA1(module,0);
	  fPHOSDAPtr = new AliPHOSRcuDA1(module,0);
	}
      
      xOffset = cellDataPtr->fRcuX*NXCOLUMNSRCU;
      zOffset = cellDataPtr->fRcuZ*NZROWSRCU;
      
      fShmPtr->SetMemory(cellDataPtr);
      currentChannel = fShmPtr->NextChannel();
      
      while(currentChannel != 0)
	{
      	  Int_t tmpZ =  currentChannel->fZ;
	  Int_t tmpX =  currentChannel->fX;
	  Int_t tmpGain =  currentChannel->fGain;
	  
	  energyArray[tmpX+xOffset][tmpZ+zOffset][tmpGain] = currentChannel->fEnergy;
	  timeArray[tmpX+xOffset][tmpZ+zOffset][tmpGain] = currentChannel->fTime;
	}
//       for(Int_t x = 0; x < NXCOLUMNSRCU; x++)
// 	{
// 	  for(Int_t z = 0; z < NZROWSRCU; z++)
// 	    {
// 	      for(Int_t gain = 0; gain < NGAINS; gain++)
// 		{
// 		  energyArray[x+xOffset][z+zOffset][gain] = cellDataPtr->fValidData[x][z][gain].fEnergy;
// 		  timeArray[x+xOffset][z+zOffset][gain] = cellDataPtr->fValidData[x][z][gain].fTime;
// 		}
// 	    }
// 	}
      iter = GetNextInputBlock(); 
    }
  
  if(fPHOSDAPtr)
    {
      // fPHOSDAPtr->FillHistograms(energyArray, timeArray);
    }

  ResetArrays(energyArray, timeArray);

  //PushBack((TObject*) fCalibDataPtr,  AliHLTPHOSDefinitions::fgkEmcCalibDataType, specification);
 
  return 0; 
}

  
Int_t 
AliHLTPHOSModuleCalibrationProcessorComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) {

  cout << "AliHLTPHOSModuleCalibrationProcessorComponent::ShipDataToFXS TP0" << endl;
 
  //://:ShipDataToFXS(

  // see header file for class documentation
 
  // ** PushBack data to FXS ...
  cout << "AliHLTPHOSModuleCalibrationProcessorComponent::ShipDataToFXS TP0.1" << endl;
  Char_t filename[200];

  for(int i=0; i < 200; i++)
    {
      filename[i] = 0;
    }

  cout << "AliHLTPHOSModuleCalibrationProcessorComponent::ShipDataToFXS TP0.2" << endl;
 
  fPHOSDAPtr->GetModule();

  sprintf(filename, "PHOS_Module%d_Calib_HLT.root",  fPHOSDAPtr->GetModule() );
  cout << "AliHLTPHOSModuleCalibrationProcessorComponent::ShipDataToFXS TP0.3" << endl;
  cout << "AliHLTPHOSModuleCalibrationProcessorComponent::ShipDataToFXS TP1" << endl;

  for(Int_t x = 0; x < NXCOLUMNSRCU; x++)
    {
      for(Int_t z = 0; z < NZROWSRCU; z++)
	{
	  cout << "AliHLTPHOSModuleCalibrationProcessorComponent::ShipDataToFXS TP2" << endl;
	  PushToFXS( (TObject*) fPHOSDAPtr->GetHgLgRatioHistogram(x, z), "PHOS", filename);
	  for(Int_t gain = 0; gain < NGAINS; gain++)
	    {
	      cout << "AliHLTPHOSModuleCalibrationProcessorComponent::ShipDataToFXS TP3" << endl; 
	      PushToFXS( (TObject*) fPHOSDAPtr->GetTimeEnergyHistogram(x, z, gain), "PHOS", filename);
	      cout << "AliHLTPHOSModuleCalibrationProcessorComponent::ShipDataToFXS TP4" << endl;
	    }
	}
    }
  return 0;
}  

void
AliHLTPHOSModuleCalibrationProcessorComponent::ResetArrays(Float_t e[NXCOLUMNSMOD][NZROWSMOD][NGAINS], Float_t t[NXCOLUMNSMOD][NZROWSMOD][NGAINS])
{
  for(Int_t x = 0; x < NXCOLUMNSRCU; x++)
    {
      for(Int_t z = 0; z < NZROWSRCU; z++)
	{
	  for(Int_t gain = 0; gain < NGAINS; gain++)
	    {
	      e[x][z][gain] = 0;
	      t[x][z][gain] = 0;
	    }
	}
    }
}
