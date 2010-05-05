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
#include "AliHLTPHOSSharedMemoryInterfacev2.h"
#include "AliHLTPHOSRcuDAComponent.h"
#include "AliHLTPHOSDefinitions.h"

// #include "AliHLTPHOSConstant.h"
//#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "TObjArray.h"
#include "AliHLTPHOSUtilities.h"
#include "AliHLTReadoutList.h"

//#include <iostream>

/** @file   AliHLTPHOSRcuDAComponent.cxx
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

AliHLTPHOSRcuDAComponent gAliHLTPHOSRcuDAComponent;

AliHLTPHOSRcuDAComponent::AliHLTPHOSRcuDAComponent() : 
  //AliHLTPHOSRcuProperties(),
  AliHLTCalibrationProcessor(),
  fPhosEventCount(0),
  fPHOSDAPtr(0),
  fShmPtr(0) 
  //    fTest(-2)
{
  fShmPtr = new AliHLTPHOSSharedMemoryInterfacev2();
}



AliHLTPHOSRcuDAComponent::~AliHLTPHOSRcuDAComponent() 
{
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


void AliHLTPHOSRcuDAComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkCellEnergyDataType);
}


AliHLTComponentDataType AliHLTPHOSRcuDAComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkEmcCalibDataType;
}
  
                                   
void AliHLTPHOSRcuDAComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  constBase = 0;
  inputMultiplier = 2;
}


AliHLTComponent* 
AliHLTPHOSRcuDAComponent::Spawn()
{
  return new AliHLTPHOSRcuDAComponent();
}

const char* 
AliHLTPHOSRcuDAComponent::GetComponentID()
{
  return "PhosRcuDAProcessor";  
}


Int_t
AliHLTPHOSRcuDAComponent::ScanArgument( Int_t /*argc*/, const char** /*argv*/)
{
  //CRAP PTH
  //  AliHLTPHOSUtilities::ScanArguments(argc, argv);
 

 return 0;
}


Int_t AliHLTPHOSRcuDAComponent::InitCalibration()
{  
  return 0;
}


Int_t AliHLTPHOSRcuDAComponent::DeinitCalibration()
{
  AliHLTComponentEventData dummyEvtData;
  AliHLTComponentTriggerData dummyTrgData;
  ShipDataToFXS(dummyEvtData, dummyTrgData); 


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



Int_t AliHLTPHOSRcuDAComponent::ProcessCalibration(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  /*
  fPhosEventCount ++;
  const  AliHLTComponentEventData eDta  = evtData;
  AliHLTComponentTriggerData  tDta =  trigData;

  UInt_t specification = 0;
  const AliHLTComponentBlockData* iter = 0;
  iter = GetFirstInputBlock( AliHLTPHOSDefinitions::fgkCellEnergyDataType | kAliHLTDataOriginPHOS);

  AliHLTPHOSRcuCellEnergyDataStruct* cellDataPtr = 0;
  AliHLTPHOSValidCellDataStruct *currentChannel =0;
  Int_t xOffset = 0;
  Int_t zOffset = 0;
  Int_t module = -1;

  Float_t energyArray[NXCOLUMNSMOD][NZROWSMOD][NGAINS];
  Float_t timeArray[NXCOLUMNSMOD][NZROWSMOD][NGAINS];
  ResetArrays(energyArray, timeArray);

  while(iter != 0)
    {
      specification = specification|iter->fSpecification;
      cellDataPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)( iter->fPtr);
      module = cellDataPtr->fModuleID;
      xOffset = cellDataPtr->fRcuX*NXCOLUMNSRCU;
      zOffset = cellDataPtr->fRcuZ*NZROWSRCU;

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
  
  fPHOSDAPtr->FillHistograms(energyArray, timeArray);

  ResetArrays(energyArray, timeArray);
  */
  return 0; 
}

  
Int_t 
AliHLTPHOSRcuDAComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) 
{
  Char_t filename[200];

  for(int i=0; i < 200; i++)
    {
      filename[i] = 0;
    }
  const TObjArray *calibPtr = fPHOSDAPtr->GetHistoContainer();
  sprintf(filename, "/home/perthi/hlt/rundir/test/outdata/%s.root", fPHOSDAPtr->GetName() );
  TFile *outFile =  new TFile(filename, "recreate");
  calibPtr->Write(); 
  outFile->Close();
  static AliHLTReadoutList rdList(AliHLTReadoutList::kPHOS);
  PushToFXS( (TObject*)fPHOSDAPtr->GetHistoContainer(), "PHOS",  filename, rdList.Buffer());
  cout << "Finnished pushing data to HLT FXS" << endl;
  return 0;
}  



/*
void
AliHLTPHOSRcuDAComponent::ResetArrays(Float_t e[NXCOLUMNSMOD][NZROWSMOD][NGAINS], Float_t t[NXCOLUMNSMOD][NZROWSMOD][NGAINS])
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
*/
