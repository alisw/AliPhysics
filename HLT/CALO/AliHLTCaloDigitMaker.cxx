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

#include "AliHLTCaloConstants.h"
#include "AliHLTCaloMapper.h"

#include "AliHLTCaloChannelDataStruct.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloSharedMemoryInterfacev2.h" // added by PTH
//#include "AliPHOSEMCAGeometry.h"
#include "TH2F.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTCaloConstants.h"

ClassImp(AliHLTCaloDigitMaker);

//using namespace CaloHLTConst;

AliHLTCaloDigitMaker::AliHLTCaloDigitMaker(TString det) :
  AliHLTCaloConstantsHandler(det),
  fShmPtr(0),
  fDigitStructPtr(0),
  fDigitCount(0),
  fOrdered(true),
  fMapperPtr(0),
  fHighGainFactors(0),
  fLowGainFactors(0),
  fBadChannelMask(0),
  fChannelBook(0)
{
  //BALLE BALLE
  
  //Must set this in the child instance
  
  // See header file for documentation

  fShmPtr = new AliHLTCaloSharedMemoryInterfacev2();

  fHighGainFactors = new Float_t*[fCaloConstants->GetNXCOLUMNSMOD()];
  fLowGainFactors = new Float_t*[fCaloConstants->GetNXCOLUMNSMOD()];

  fBadChannelMask = new Float_t**[fCaloConstants->GetNXCOLUMNSMOD()];

  fChannelBook= new AliHLTCaloDigitDataStruct**[fCaloConstants->GetNXCOLUMNSMOD()];

  for(int x = 0; x < fCaloConstants->GetNXCOLUMNSMOD(); x++)
    {
      fHighGainFactors[x] = new Float_t[fCaloConstants->GetNZROWSMOD()];
      fLowGainFactors[x] = new Float_t[fCaloConstants->GetNZROWSMOD()];

      fBadChannelMask[x] = new Float_t*[fCaloConstants->GetNZROWSMOD()];

      fChannelBook[x] = new AliHLTCaloDigitDataStruct*[fCaloConstants->GetNZROWSMOD()];

      for(int z = 0; z < fCaloConstants->GetNZROWSMOD(); z++)
	{

	  fHighGainFactors[x][z] = 0.005;
	  fLowGainFactors[x][z] = 0.08;
	 
	  fBadChannelMask[x][z] = new Float_t[fCaloConstants->GetNGAINS()];
	  fBadChannelMask[x][z][fCaloConstants->GetHIGHGAIN()] = 1;
	  fBadChannelMask[x][z][fCaloConstants->GetLOWGAIN()] = 1; 
	  
	  fChannelBook[x][z] = 0;
	  
	}
    }	  
  
  //Must be set in child instance
//fMapperPtr = new AliHLTCaloMapper(det);
}
   
AliHLTCaloDigitMaker::~AliHLTCaloDigitMaker() 
{
  //See header file for documentation
}

Int_t
AliHLTCaloDigitMaker::MakeDigits(AliHLTCaloChannelDataHeaderStruct* channelDataHeader, AliHLTUInt32_t availableSize)
{
  //See header file for documentation
  
  Int_t j = 0;
  UInt_t totSize = sizeof(AliHLTCaloDigitDataStruct);
  
//   Int_t xMod = -1;
//   Int_t zMod = -1;
  
  UShort_t coord1[4];
  //  UShort_t coord2[4];
  Float_t locCoord[3];
  
  
  AliHLTCaloChannelDataStruct* currentchannel = 0;
  //  AliHLTCaloChannelDataStruct* currentchannelLG = 0;  
  AliHLTCaloChannelDataStruct* tmpchannel = 0;
  
  fShmPtr->SetMemory(channelDataHeader);
  currentchannel = fShmPtr->NextChannel();

  while(currentchannel != 0)
    {
      if(availableSize < totSize) return -1;

      fMapperPtr->GetChannelCoord(currentchannel->fChannelID, coord1);
      
      tmpchannel = currentchannel;
	  
      if(coord1[2] == fCaloConstants->GetHIGHGAIN()) // We got a completely new crystal
	{
	  fMapperPtr->GetLocalCoord(currentchannel->fChannelID, locCoord);
	  if(UseDigit(coord1, currentchannel))
 	    {
	      AddDigit(currentchannel, coord1, locCoord);
	      j++;	      
	      totSize += sizeof(AliHLTCaloDigitDataStruct);
	    }
	  currentchannel = fShmPtr->NextChannel(); // Get the next channel
	}
      else if(coord1[2] == fCaloConstants->GetLOWGAIN())
	{
	  fMapperPtr->GetLocalCoord(currentchannel->fChannelID, locCoord);
	  if(UseDigit(coord1, currentchannel))
	    {
	      AddDigit(currentchannel, coord1, locCoord);
	      j++;	      
	      totSize += sizeof(AliHLTCaloDigitDataStruct);
	    }
	  currentchannel = fShmPtr->NextChannel(); // Get the next channel
	}
    }

  fDigitCount += j;
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
	      fBadChannelMask[x][z][fCaloConstants->GetHIGHGAIN()] = 1;
	    }
	  else
	    {
	      fBadChannelMask[x][z][fCaloConstants->GetHIGHGAIN()] = 0;
	    }
	  if(badChannelLGHist->GetBinContent(x, z) < qCut && badChannelLGHist->GetBinContent(x, z) > 0)
	    {
	      fBadChannelMask[x][z][fCaloConstants->GetLOWGAIN()] = 0;
	    }
	  else
	    {
	      fBadChannelMask[x][z][fCaloConstants->GetLOWGAIN()] = 0;
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


void AliHLTCaloDigitMaker::AddDigit(AliHLTCaloChannelDataStruct* channelData, UShort_t* channelCoordinates, Float_t* localCoordinates)
{

  fChannelBook[channelCoordinates[0]][channelCoordinates[0]] = fDigitStructPtr;

  fDigitStructPtr->fX = channelCoordinates[0];
  fDigitStructPtr->fZ = channelCoordinates[1];

  fDigitStructPtr->fLocX = localCoordinates[0];
  fDigitStructPtr->fLocZ = localCoordinates[1];

  if(channelCoordinates[2] == fCaloConstants->GetHIGHGAIN() )
    {
      fDigitStructPtr->fEnergy = channelData->fEnergy*fHighGainFactors[channelCoordinates[0]][channelCoordinates[1]];
      if(channelData->fEnergy >= 1023)
	{
	  fDigitStructPtr->fOverflow = true;
	}
      //	printf("HG channel (x = %d, z = %d) with amplitude: %f --> Digit with energy: %f \n", channelCoordinates[0], channelCoordinates[1], channelData->fEnergy, fDigitStructPtr->fEnergy);
    }
  else
    {
      fDigitStructPtr->fEnergy = channelData->fEnergy*fLowGainFactors[channelCoordinates[0]][channelCoordinates[1]];
      if(channelData->fEnergy >= 1023)
	{
	  fDigitStructPtr->fOverflow = true;
	}
      //	printf("LG channel (x = %d, z = %d) with amplitude: %f --> Digit with energy: %f\n", channelCoordinates[0], channelCoordinates[1], channelData->fEnergy, fDigitStructPtr->fEnergy); 
    }
  fDigitStructPtr->fTime = channelData->fTime * 0.0000001; //TODO
  fDigitStructPtr->fCrazyness = channelData->fCrazyness;
  fDigitStructPtr->fModule = channelCoordinates[3];
  fDigitStructPtr++;
}

bool AliHLTCaloDigitMaker::UseDigit(UShort_t *channelCoordinates, AliHLTCaloChannelDataStruct *channel) 
{
  AliHLTCaloDigitDataStruct *tmpDigit = fChannelBook[channelCoordinates[0]][channelCoordinates[1]];
  if(tmpDigit)
    {
      if(channelCoordinates[2] == fCaloConstants->GetLOWGAIN())
	{
	  if(tmpDigit->fOverflow)
	    {
	      return true;
	    }
	  return false;
	}
      else
	{
	  if(channel->fEnergy >= fCaloConstants->GetMAXBINVALUE() )
	    {
	      return false;
	    }
	  return true;
	}
    }
  return true;
}
