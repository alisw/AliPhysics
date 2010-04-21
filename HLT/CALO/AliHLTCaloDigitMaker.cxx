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
 /** 
 * @file   AliHLTCALOClusterizer.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Digit maker for CALO HLT  
 */
  

    

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloDigitMaker.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTCaloMapper.h"
#include "AliHLTCaloChannelDataStruct.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloCoordinate.h"
#include "AliHLTCaloSharedMemoryInterfacev2.h" // added by PTH
//#include "AliPHOSEMCAGeometry.h"
#include "TH2F.h"
#include "AliHLTCaloConstants.h"

ClassImp(AliHLTCaloDigitMaker);

//using namespace CaloHLTConst;

AliHLTCaloDigitMaker::AliHLTCaloDigitMaker(TString det) :
  AliHLTCaloConstantsHandler(det),
  fShmPtr(0),
  fDigitStructPtr(0),
  fDigitCount(0),
  fMapperPtr(0),
  fHighGainFactors(0),
  fLowGainFactors(0),
  fBadChannelMask(0),
  fChannelBook(0),
  fMaxEnergy(900)
{
  // See header file for documentation

  fShmPtr = new AliHLTCaloSharedMemoryInterfacev2(det);

  fHighGainFactors = new Float_t*[fCaloConstants->GetNXCOLUMNSMOD()];
  fLowGainFactors = new Float_t*[fCaloConstants->GetNXCOLUMNSMOD()];
  
  fBadChannelMask = new Bool_t**[fCaloConstants->GetNXCOLUMNSMOD()];
  
  fChannelBook= new AliHLTCaloDigitDataStruct**[fCaloConstants->GetNXCOLUMNSMOD()];
  
  for(int x = 0; x < fCaloConstants->GetNXCOLUMNSMOD(); x++)
    {
      fHighGainFactors[x] = new Float_t[fCaloConstants->GetNZROWSMOD()];
      fLowGainFactors[x] = new Float_t[fCaloConstants->GetNZROWSMOD()];

      fBadChannelMask[x] = new Bool_t*[fCaloConstants->GetNZROWSMOD()];

      fChannelBook[x] = new AliHLTCaloDigitDataStruct*[fCaloConstants->GetNZROWSMOD()];

      for(int z = 0; z < fCaloConstants->GetNZROWSMOD(); z++)
	{

	  fHighGainFactors[x][z] = 0.005;
	  fLowGainFactors[x][z] = 0.08;
	 
	  fBadChannelMask[x][z] = new Bool_t[fCaloConstants->GetNGAINS()];
	  fBadChannelMask[x][z][fCaloConstants->GetHIGHGAIN()] = false;
	  fBadChannelMask[x][z][fCaloConstants->GetLOWGAIN()] = false; 
	  
	  fChannelBook[x][z] = 0;
	  
	}
    }	  
}

AliHLTCaloDigitMaker::~AliHLTCaloDigitMaker() 
{
  //See header file for documentation
}

Int_t
AliHLTCaloDigitMaker::MakeDigits(AliHLTCaloChannelDataHeaderStruct* channelDataHeader, AliHLTUInt32_t availableSize)
{
  //See header file for documentation
  
  Reset();

  UInt_t totSize = sizeof(AliHLTCaloDigitDataStruct);
  
//   Int_t xMod = -1;
//   Int_t zMod = -1;
  

  AliHLTCaloCoordinate coord;
  
  
  AliHLTCaloChannelDataStruct* currentchannel = 0;
  
  fShmPtr->SetMemory(channelDataHeader);
  currentchannel = fShmPtr->NextChannel();

  while(currentchannel != 0)
    {
      if(availableSize < totSize) return -1;

      fMapperPtr->ChannelId2Coordinate(currentchannel->fChannelID, coord);
      
      //      fMapperPtr->GetLocalCoord(currentchannel->fChannelID, locCoord);
      if(UseDigit(coord, currentchannel))
      {
	AddDigit(currentchannel, coord);
	//	j++;	      
	totSize += sizeof(AliHLTCaloDigitDataStruct);
      }
      currentchannel = fShmPtr->NextChannel(); // Get the next channel
    }
//       if(currentchannel)
// 	{
// 	  fMapperPtr->GetLocalCoord(currentchannel->fChannelID, locCoord);
// 	  if(UseDigit(coord1, currentchannel))
// 	    {
// 	      AddDigit(currentchannel, coord1, locCoord);
// 	      j++;	      
// 	      totSize += sizeof(AliHLTCaloDigitDataStruct);
// 	    }
// 	  currentchannel = fShmPtr->NextChannel(); // Get the next channel
// 	}
//     }
   
//   fDigitCount += j;
   return fDigitCount; 
}

void 
AliHLTCaloDigitMaker::SetGlobalHighGainFactor(Float_t factor)
{
  //See header file for documentation
  for(int x = 0; x < fCaloConstants->GetNXCOLUMNSMOD(); x++)
    {
      for(int z = 0; z < fCaloConstants->GetNZROWSMOD(); z++)
	{
	  fHighGainFactors[x][z] = factor;
	}
    }
}

void
AliHLTCaloDigitMaker::SetGlobalLowGainFactor(Float_t factor)
{
  //See header file for documentation
  for(int x = 0; x < fCaloConstants->GetNXCOLUMNSMOD(); x++)
    {
      for(int z = 0; z < fCaloConstants->GetNZROWSMOD(); z++)
	{
	  fLowGainFactors[x][z] = factor;
	}
    }
}

void
AliHLTCaloDigitMaker::SetBadChannelMask(TH2F* badChannelHGHist, TH2F* badChannelLGHist, Float_t qCut)
{
 for(int x = 0; x < fCaloConstants->GetNXCOLUMNSMOD(); x++)
    {
      for(int z = 0; z < fCaloConstants->GetNZROWSMOD(); z++)
	{
	  if(badChannelHGHist->GetBinContent(x, z) < qCut && badChannelHGHist->GetBinContent(x, z) > 0)
	    {
	      fBadChannelMask[x][z][fCaloConstants->GetHIGHGAIN()] = true;
	    }
	  else
	    {
	      fBadChannelMask[x][z][fCaloConstants->GetHIGHGAIN()] = false;
	    }
	  if(badChannelLGHist->GetBinContent(x, z) < qCut && badChannelLGHist->GetBinContent(x, z) > 0)
	    {
	      fBadChannelMask[x][z][fCaloConstants->GetLOWGAIN()] = false;
	    }
	  else
	    {
	      fBadChannelMask[x][z][fCaloConstants->GetLOWGAIN()] = false;
	    }
	}
    }
}

void
AliHLTCaloDigitMaker::Reset()
{
  fDigitCount = 0;
  for(int x = 0; x < fCaloConstants->GetNXCOLUMNSMOD(); x++)
    {
      for(int z = 0; z < fCaloConstants->GetNZROWSMOD(); z++)
	{
	  fChannelBook[x][z] = 0;
	}
    }	  

}


void AliHLTCaloDigitMaker::AddDigit(AliHLTCaloChannelDataStruct* channelData, AliHLTCaloCoordinate &coord)
{

  AliHLTCaloDigitDataStruct *tmpDigit = fDigitStructPtr + 1;

  if(fChannelBook[coord.fX][coord.fZ])
    {
      tmpDigit = fDigitStructPtr;
      fDigitStructPtr = fChannelBook[coord.fX][coord.fZ];
      fDigitCount--;
      //      printf("Going to overwrite digit: x = %d, z = %d, gain = %d, energy = %f\n", fDigitStructPtr->fX, fDigitStructPtr->fZ, fDigitStructPtr->fGain, fDigitStructPtr->fEnergy);
    }
  
  fChannelBook[coord.fX][coord.fZ] = fDigitStructPtr;

 
  fDigitStructPtr->fX = coord.fX;
  fDigitStructPtr->fZ = coord.fZ;
  fDigitStructPtr->fGain = coord.fGain;
  fDigitStructPtr->fOverflow = false;
  
  fDigitStructPtr->fID = fDigitStructPtr->fZ * fCaloConstants->GetNXCOLUMNSMOD() + fDigitStructPtr->fX;

  if(coord.fGain == fCaloConstants->GetHIGHGAIN() )
    {
      fDigitStructPtr->fEnergy = channelData->fEnergy*fHighGainFactors[coord.fX][coord.fZ];
      if(channelData->fEnergy >= fMaxEnergy)
	{
	  fDigitStructPtr->fOverflow = true;
	}
            //printf("HG channel (x = %d, z = %d) with amplitude: %f --> Digit with energy: %f \n", coord.fX, coord.fZ, channelData->fEnergy, fDigitStructPtr->fEnergy);
    }
  else
    {
      fDigitStructPtr->fEnergy = channelData->fEnergy*fLowGainFactors[coord.fX][coord.fZ];
      if(channelData->fEnergy >= 1023)
	{
	  fDigitStructPtr->fOverflow = true;
	}
    //        printf("LG channel (x = %d, z = %d) with amplitude: %f --> Digit with energy: %f\n", coord.fX, coord.fZ, channelData->fEnergy, fDigitStructPtr->fEnergy); 
    }
  fDigitStructPtr->fTime = channelData->fTime * 0.0000001; //TODO
  fDigitStructPtr->fCrazyness = channelData->fCrazyness;
  fDigitStructPtr->fModule = coord.fModuleId;
  fDigitStructPtr = tmpDigit;
  //  fDigitStructPtr++;
  fDigitCount++;
}

bool AliHLTCaloDigitMaker::UseDigit(AliHLTCaloCoordinate &channelCoordinates, AliHLTCaloChannelDataStruct *channel) 
{
   if(fBadChannelMask[channelCoordinates.fX][channelCoordinates.fZ][0] == true) return false;
  AliHLTCaloDigitDataStruct *tmpDigit = fChannelBook[channelCoordinates.fX][channelCoordinates.fZ];
  //printf("UseDigit: Got digit, x: %d, z: %d, gain: %d, amp: %f\n", channelCoordinates.fX, channelCoordinates.fZ, channelCoordinates.fGain, channel->fEnergy);
  if(tmpDigit)
    {
      if(channelCoordinates.fGain == fCaloConstants->GetLOWGAIN())
	{
	  	  //printf("UseDigit: Already have digit with, x: %d, z: %d, with high gain \n", channelCoordinates.fX, channelCoordinates.fZ);
	  if(tmpDigit->fOverflow)
	    {
	      	     // printf("But it was in overflow! Let's use this low gain!\n");
	      return true;
	    }
	  return false;
	}
      else
	{
	  	  //printf("UseDigit: Already have digit with, x: %d, z: %d, with low gain: %d\n", channelCoordinates.fX, channelCoordinates.fZ);
	  if(channel->fEnergy > fMaxEnergy )
	    {
	      return false;
	    }
	  return true;
	}
    }
  return true;
}

void AliHLTCaloDigitMaker::SetBadChannel(Int_t x, Int_t z, Bool_t bad)
{
   // See header file for class documentation
   fBadChannelMask[x][z][0] = bad;
   fBadChannelMask[x][z][1] = bad;
}

void AliHLTCaloDigitMaker::SetGain(Int_t x, Int_t z, Float_t ratio, Float_t gain)
{
   // See header file for class documentation
   
   fHighGainFactors[x][z] = gain;
   fLowGainFactors[x][z] = gain * ratio;
   
}
