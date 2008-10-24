/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 s*                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 /** 
 * @file   AliHLTPHOSClusterizer.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Digit maker for PHOS HLT  
 */
  

    

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSDigitMaker.h"
#include "AliHLTPHOSDigit.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSBaseline.h"
#include "AliHLTPHOSMapper.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH2F.h"

#include "AliHLTPHOSValidCellDataStruct.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSChannelDataStruct.h"
#include "AliHLTPHOSChannelDataHeaderStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "AliHLTPHOSSharedMemoryInterfacev2.h" // added by PTH


ClassImp(AliHLTPHOSDigitMaker);

using namespace PhosHLTConst;

AliHLTPHOSDigitMaker::AliHLTPHOSDigitMaker() :
  AliHLTPHOSBase(),
  fCellDataPtr(0),
  fDigitContainerStructPtr(0),
  fDigitArrayPtr(0),
  fDigitStructPtr(0),
  fShmPtr(0),
  fDigitCount(0),
  fOrdered(true)
{
  // See header file for documentation

  fShmPtr = new AliHLTPHOSSharedMemoryInterfacev2();

  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(int z = 0; z < N_ZROWS_MOD; z++)
	{
	  fHighGainFactors[x][z] = 0.005;
	  fLowGainFactors[x][z] = 0.08;
	  fDigitThresholds[x][z][HIGH_GAIN] = 2;
	  fDigitThresholds[x][z][LOW_GAIN] = 2;
	  fBadChannelMask[x][z][HIGH_GAIN] = 1;
	  fBadChannelMask[x][z][LOW_GAIN] = 1; 
	}
    }

    for(int x = 0; x < N_XCOLUMNS_RCU; x++)
      {
	for(int z = 0; z < N_ZROWS_RCU; z++)
	  {
	    for(int gain = 0; gain < N_GAINS; gain++)
	      {
		fEnergyArray[x][z][gain] = 0;
	      }
	  }
      }
	  
  
}
   
AliHLTPHOSDigitMaker::~AliHLTPHOSDigitMaker() 
{
  //See header file for documentation
}

Int_t
AliHLTPHOSDigitMaker::MakeDigits(AliHLTPHOSChannelDataHeaderStruct* channelDataHeader)
{
  //See header file for documentation
  
  Int_t j = 0;
  //  Float_t tmpEnergy = 0;

//   Int_t xMod = -1;
//   Int_t zMod = -1;
  
  AliHLTPHOSMapper mapper;
  UShort_t coord1[4];
  UShort_t coord2[4];
  
  
  AliHLTPHOSChannelDataStruct* currentchannel = 0;
  AliHLTPHOSChannelDataStruct* currentchannelLG = 0;  
  AliHLTPHOSChannelDataStruct* tmpchannel = 0;

  Reset();

  fShmPtr->SetMemory(channelDataHeader);
  currentchannel = fShmPtr->NextChannel();
  
//   while(currentchannel != 0)
//     {
//       fEnergyArray[currentchannel->fX][currentchannel->fZ][currentchannel->fGain] = currentchannel->fEnergy; 
//       currentchannel = fShmPtr->NextChannel();
//     }

 
  while(currentchannel != 0)
    {
      AliHLTPHOSMapper::GetChannelCoord(currentchannel->fChannelID, coord1);
      if(fOrdered)
	{
	  tmpchannel = currentchannel;
	  if(coord1[2] == HIGH_GAIN) 
	    {

	      if(currentchannel->fEnergy < MAX_BIN_VALUE)
		{
		  fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
		  fDigitStructPtr->fX = coord1[0];
		  fDigitStructPtr->fZ = coord1[1];
		  fDigitStructPtr->fAmplitude = currentchannel->fEnergy;
		  fDigitStructPtr->fEnergy = currentchannel->fEnergy * fHighGainFactors[coord1[0]][coord1[1]];
		  //TODO: fix time //CRAP
		  fDigitStructPtr->fTime = currentchannel->fTime * 0.0000001;
		  fDigitStructPtr->fCrazyness = currentchannel->fCrazyness;
		  fDigitStructPtr->fModule = coord1[3];
		  j++;	      
		  currentchannel = fShmPtr->NextChannel();
		  AliHLTPHOSMapper::GetChannelCoord(currentchannel->fChannelID, coord2);
		  if(coord1[0] == coord2[0] && coord1[1] == coord2[1])
		    {
		      fShmPtr->NextChannel();
		    }
		}
	      else
		{
		  currentchannel = fShmPtr->NextChannel();
		  AliHLTPHOSMapper::GetChannelCoord(currentchannel->fChannelID, coord2);
		  if(coord2[2] == LOW_GAIN)
		    {
		      fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
		      fDigitStructPtr->fX = coord2[0];
		      fDigitStructPtr->fZ = coord2[1];
		      fDigitStructPtr->fAmplitude = currentchannel->fEnergy;
		      fDigitStructPtr->fEnergy = currentchannel->fEnergy * fHighGainFactors[coord2[0]][coord2[1]];
		      //TODO: fix time //CRAP
		      fDigitStructPtr->fTime = currentchannel->fTime * 0.0000001;
		      fDigitStructPtr->fCrazyness = currentchannel->fCrazyness;
		      fDigitStructPtr->fModule = coord2[3];
		      j++;
		      currentchannel = fShmPtr->NextChannel();		      
		    }
		  else
		    {
		      fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
		      fDigitStructPtr->fX = coord1[0];
		      fDigitStructPtr->fZ = coord1[1];
		      fDigitStructPtr->fAmplitude = tmpchannel->fEnergy;
		      fDigitStructPtr->fEnergy = tmpchannel->fEnergy * fHighGainFactors[coord1[0]][coord1[1]];
		      //TODO: fix time //CRAP
		      fDigitStructPtr->fTime = tmpchannel->fTime * 0.0000001;
		      fDigitStructPtr->fCrazyness = 100;
		      fDigitStructPtr->fModule = coord1[3];
		      j++;	      
		    }
		    
		}
	    }
	  else
	    {      
	      fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
	      fDigitStructPtr->fX = coord1[0];
	      fDigitStructPtr->fZ = coord1[1];
	      fDigitStructPtr->fAmplitude = tmpchannel->fEnergy;
	      fDigitStructPtr->fEnergy = tmpchannel->fEnergy * fHighGainFactors[coord1[0]][coord1[1]];
	      //TODO: fix time //CRAP
	      fDigitStructPtr->fTime = tmpchannel->fTime * 0.0000001;
	      fDigitStructPtr->fCrazyness = 100;
	      fDigitStructPtr->fModule = coord1[3];
	      j++;	      
	      currentchannel = fShmPtr->NextChannel();		      
	    }

	}
      else 
	{
	  if(coord1[2] == LOW_GAIN)
	    {
	      currentchannelLG = currentchannel;
	      currentchannel = fShmPtr->NextChannel();
	      AliHLTPHOSMapper::GetChannelCoord(currentchannel->fChannelID, coord2);
	      
	      if(coord1[0] == coord2[0] && coord1[1] == coord2[1])
		{
		  if(currentchannel->fEnergy < MAX_BIN_VALUE) 
		    {
		      fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
		      fDigitStructPtr->fX = coord2[0];
		      fDigitStructPtr->fZ = coord2[1];
		      fDigitStructPtr->fAmplitude = currentchannel->fEnergy;
		      fDigitStructPtr->fEnergy = currentchannel->fEnergy * fHighGainFactors[coord2[0]][coord2[1]];
		      //TODO: fix time //CRAP
		      fDigitStructPtr->fTime = currentchannel->fTime * 0.0000001;
		      fDigitStructPtr->fCrazyness = currentchannel->fCrazyness;
		      fDigitStructPtr->fModule = coord2[3];
		      j++;
		      currentchannel = fShmPtr->NextChannel();
		    }
		  else
		    {
		      fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
		      fDigitStructPtr->fX = coord1[0];
		      fDigitStructPtr->fZ = coord1[1];
		      fDigitStructPtr->fAmplitude = currentchannelLG->fEnergy;
		      fDigitStructPtr->fEnergy = currentchannelLG->fEnergy * fHighGainFactors[coord1[0]][coord1[1]];
		      //TODO: fix time //CRAP
		      fDigitStructPtr->fTime = currentchannelLG->fTime * 0.0000001;
		      fDigitStructPtr->fCrazyness = currentchannel->fCrazyness;
		      fDigitStructPtr->fModule = coord1[3];
		      j++;
		      currentchannel = fShmPtr->NextChannel();
		    }
		}
	      else
		{
		  fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
		  fDigitStructPtr->fX = coord1[0];
		  fDigitStructPtr->fZ = coord1[1];
		  fDigitStructPtr->fAmplitude = currentchannelLG->fEnergy;
		  fDigitStructPtr->fEnergy = currentchannelLG->fEnergy * fHighGainFactors[coord1[0]][coord1[1]];
		  //TODO: fix time //CRAP
		  fDigitStructPtr->fTime = currentchannelLG->fTime * 0.0000001;
		  fDigitStructPtr->fCrazyness = currentchannel->fCrazyness;
		  fDigitStructPtr->fModule = coord1[3];
		  j++;
		}
	    }
	  else
	    {
	      fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
	      fDigitStructPtr->fX = coord1[0];
	      fDigitStructPtr->fZ = coord1[1];
	      fDigitStructPtr->fAmplitude = currentchannel->fEnergy;
	      fDigitStructPtr->fEnergy = currentchannel->fEnergy * fHighGainFactors[coord1[0]][coord1[1]];
	      //TODO: fix time //CRAP
	      fDigitStructPtr->fTime = currentchannel->fTime * 0.0000001;
	      fDigitStructPtr->fCrazyness = currentchannel->fCrazyness;
	      fDigitStructPtr->fModule = coord2[3];
	      j++;
	      fShmPtr->NextChannel();
	    }
	}
    }




//   for(int x = 0; x < N_XCOLUMNS_RCU; x++)
//     {
//       for(int z = 0; z < N_ZROWS_RCU; z++)  
// 	{
// 	  xMod = x + rcuData->fRcuX * N_XCOLUMNS_RCU;
// 	  zMod = z + rcuData->fRcuZ * N_ZROWS_RCU;
	  
// 	  if(fEnergyArray[x][z][HIGH_GAIN] > fDigitThresholds[xMod][zMod][HIGH_GAIN] && fEnergyArray[x][z][HIGH_GAIN] < MAX_BIN_VALUE && fBadChannelMask[xMod][zMod][HIGH_GAIN])
// 	    {
// 	      fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
// 	      fDigitStructPtr->fX = xMod;
// 	      fDigitStructPtr->fZ = zMod;
// 	      fDigitStructPtr->fAmplitude = fEnergyArray[x][z][HIGH_GAIN];
// 	      fDigitStructPtr->fEnergy = fEnergyArray[x][z][HIGH_GAIN] * fHighGainFactors[xMod][zMod];
// 	      //TODO: fix time //CRAP
// 	      fDigitStructPtr->fTime = fCellDataPtr->fTime * 0.0000001;
// 	      fDigitStructPtr->fCrazyness = fCellDataPtr->fCrazyness;
// 	      fDigitStructPtr->fModule = rcuData->fModuleID;
// 	      j++;
// 	    }
// 	  else if(fEnergyArray[x][z][LOW_GAIN] >= MAX_BIN_VALUE && fBadChannelMask[xMod][zMod][LOW_GAIN])
// 	    {
// 	      if(fEnergyArray[x][z][LOW_GAIN] > fDigitThresholds[xMod][zMod][LOW_GAIN])
// 		{	      
// 		  fDigitStructPtr = &(fDigitContainerStructPtr->fDigitDataStruct[j+fDigitCount]);
// 		  fDigitStructPtr->fX = xMod;
// 		  fDigitStructPtr->fZ = zMod;
// 		  fDigitStructPtr->fAmplitude = fEnergyArray[x][z][LOW_GAIN];
// 		  fDigitStructPtr->fEnergy = fEnergyArray[x][z][LOW_GAIN] * fLowGainFactors[xMod][zMod];
// 		  //TODO: fix time //CRAP
// 		  fDigitStructPtr->fTime = fCellDataPtr->fTime  * 0.0000001;;
// 		  fDigitStructPtr->fCrazyness = fCellDataPtr->fCrazyness;
// 		  fDigitStructPtr->fModule = rcuData->fModuleID;
// 		  j++;
// 		}
// 	    }
// 	}
//     }
  
  fDigitCount += j;
  fDigitContainerStructPtr->fNDigits = fDigitCount;
  return fDigitCount; 
}


void
AliHLTPHOSDigitMaker::Reset()
{ 
  //See header file for documentation

  for(int x = 0; x < N_XCOLUMNS_RCU; x++)
    {
      for(int z = 0; z < N_ZROWS_RCU; z++)  
	{
	  for(int gain = 0; gain < N_GAINS; gain++)
	    {
	      fEnergyArray[x][z][gain] = 0;
	    }
	}
    }

  fDigitCount = 0;
}

void 
AliHLTPHOSDigitMaker::SetGlobalHighGainFactor(Float_t factor)
{
  //See header file for documentation
  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(int z = 0; z < N_ZROWS_MOD; z++)
	{
	  fHighGainFactors[x][z] = factor;
	}
    }
}

void
AliHLTPHOSDigitMaker::SetGlobalLowGainFactor(Float_t factor)
{
  //See header file for documentation
  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(int z = 0; z < N_ZROWS_MOD; z++)
	{
	  fLowGainFactors[x][z] = factor;
	}
    }
}

void 
AliHLTPHOSDigitMaker::SetDigitThresholds(const char* filepath, Int_t nSigmas)
{
  //See header file for documentation
  TFile *histFile = new TFile(filepath);
  
  TH2F *lgHist = (TH2F*)histFile->Get("RMSLGMapHist");
  TH2F *hgHist = (TH2F*)histFile->Get("RMSHGMapHist");

  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(int z = 0; z < N_ZROWS_MOD; z++)
	{
	  fDigitThresholds[x][z][LOW_GAIN] = lgHist->GetBinContent(x, z) * nSigmas;
	  fDigitThresholds[x][z][HIGH_GAIN] = hgHist->GetBinContent(x, z) * nSigmas;
	}
    }
}

void 
AliHLTPHOSDigitMaker::SetDigitThresholds(const Float_t thresholdHG, const Float_t thresholdLG)
{
  //See header file for documentation
  for(int x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(int z = 0; z < N_ZROWS_MOD; z++)
	{
	  fDigitThresholds[x][z][LOW_GAIN] = thresholdHG;
	  fDigitThresholds[x][z][HIGH_GAIN] = thresholdLG;
	}
    }
}


void
AliHLTPHOSDigitMaker::SetBadChannelMask(TH2F* badChannelHGHist, TH2F* badChannelLGHist, Float_t qCut)
{
 for(int x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(int z = 0; z < N_ZROWS_MOD; z++)
	{
	  if(badChannelHGHist->GetBinContent(x, z) < qCut && badChannelHGHist->GetBinContent(x, z) > 0)
	    {
	      fBadChannelMask[x][z][HIGH_GAIN] = 1;
	    }
	  else
	    {
	      fBadChannelMask[x][z][HIGH_GAIN] = 0;
	    }
	  if(badChannelLGHist->GetBinContent(x, z) < qCut && badChannelLGHist->GetBinContent(x, z) > 0)
	    {
	      fBadChannelMask[x][z][LOW_GAIN] = 0;
	    }
	  else
	    {
	      fBadChannelMask[x][z][LOW_GAIN] = 0;
	    }
	}
    }
}
