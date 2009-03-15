// $Id$

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

#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSMapper.h"

#include "AliHLTPHOSChannelDataStruct.h"
#include "AliHLTPHOSChannelDataHeaderStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSSharedMemoryInterfacev2.h" // added by PTH

#include "TH2F.h"

ClassImp(AliHLTPHOSDigitMaker);

using namespace PhosHLTConst;

AliHLTPHOSDigitMaker::AliHLTPHOSDigitMaker() :
  AliHLTPHOSBase(),
  fShmPtr(0),
  fDigitStructPtr(0),
  fDigitCount(0),
  fOrdered(true),
  fMapperPtr(0)
{
  // See header file for documentation

  fShmPtr = new AliHLTPHOSSharedMemoryInterfacev2();

  for(int x = 0; x < NXCOLUMNSMOD; x++)
    {
      for(int z = 0; z < NZROWSMOD; z++)
	{
	  fHighGainFactors[x][z] = 0.005;
	  fLowGainFactors[x][z] = 0.08;
	  fBadChannelMask[x][z][HIGHGAIN] = 1;
	  fBadChannelMask[x][z][LOWGAIN] = 1; 
	}
    }	  
  fMapperPtr = new AliHLTPHOSMapper();

}
   
AliHLTPHOSDigitMaker::~AliHLTPHOSDigitMaker() 
{
  //See header file for documentation
}

Int_t
AliHLTPHOSDigitMaker::MakeDigits(AliHLTPHOSChannelDataHeaderStruct* channelDataHeader, AliHLTUInt32_t availableSize)
{
  //See header file for documentation
  
  Int_t j = 0;
  UInt_t totSize = sizeof(AliHLTPHOSDigitDataStruct);
  
//   Int_t xMod = -1;
//   Int_t zMod = -1;
  
  UShort_t coord1[4];
  UShort_t coord2[4];
  
  
  AliHLTPHOSChannelDataStruct* currentchannel = 0;
  AliHLTPHOSChannelDataStruct* currentchannelLG = 0;  
  AliHLTPHOSChannelDataStruct* tmpchannel = 0;
  
  fShmPtr->SetMemory(channelDataHeader);
  currentchannel = fShmPtr->NextChannel();

  while(currentchannel != 0)
    {
      if(availableSize < totSize) return -1;

      AliHLTPHOSMapper::GetChannelCoord(currentchannel->fChannelID, coord1);
      
      if(fOrdered) // High gain comes before low gain
      	{
	  tmpchannel = currentchannel;
	  
	  if(coord1[2] == HIGHGAIN) // We got a completely new crystal
	    {

	      if(currentchannel->fEnergy < MAXBINVALUE) // Make sure we don't have signal overflow
		{

		  AddDigit(currentchannel, coord1);
		  j++;	      
		  totSize += sizeof(AliHLTPHOSDigitDataStruct);

		  currentchannel = fShmPtr->NextChannel(); // Get the next channel

		  if(currentchannel != 0) // There was a next channel!
		    {
		      AliHLTPHOSMapper::GetChannelCoord(currentchannel->fChannelID, coord2);
		      if(coord1[0] == coord2[0] && coord1[1] == coord2[1]) // Did we get the low gain channel for this crystal?
			{
			  currentchannel = fShmPtr->NextChannel(); // In that case, jump to next channel
			}
		    }
		}

	      else // Ooops, overflow, we try the next channel... 
		{
		  currentchannel = fShmPtr->NextChannel();
		  if(currentchannel != 0) // There was a next channel
		    {
		      AliHLTPHOSMapper::GetChannelCoord(currentchannel->fChannelID, coord2);
		      if(coord2[0] == coord1[0] && coord2[1] == coord1[1]) // It is a low gain channel with the same coordinates, we may use it
			{
			  
			  AddDigit(currentchannel, coord2);
			  j++;
			  totSize += sizeof(AliHLTPHOSDigitDataStruct);
			  currentchannel = fShmPtr->NextChannel();		      
			}
		      
		      else // No low gain channel with information about the overflow channel so we just use the overflowed one...
			{
			  AddDigit(tmpchannel, coord1);
			  j++;	      
			  totSize += sizeof(AliHLTPHOSDigitDataStruct);
			  // no need to get the next channel here, we already did...
			}
		      
		    }
		}
	    }
	  else // Well, there seem to be missing a high gain channel for this crystal, let's use the low gain one
	    {      
	      AddDigit(tmpchannel, coord1);
	      j++;	      
	      totSize += sizeof(AliHLTPHOSDigitDataStruct);
	      currentchannel = fShmPtr->NextChannel(); 
	    }

	}
      else  //Reversed ordered (low gain before high gain)
	{
	  if(coord1[2] == LOWGAIN) // We got a new channel!
	    {
	      currentchannelLG = currentchannel; // Ok, let's back up the low gain channel and look for the fancy high gain one
	      currentchannel = fShmPtr->NextChannel();

	      if(currentchannel != 0) //There was another channel in the event
		{
		  AliHLTPHOSMapper::GetChannelCoord(currentchannel->fChannelID, coord2);
		  
		  if(coord1[0] == coord2[0] && coord1[1] == coord2[1]) // Aha! Found the high gain channel
		    {
		      if(currentchannel->fEnergy < MAXBINVALUE)  // To overflow or not to overflow?
			{

			  AddDigit(currentchannel, coord2);
			  j++;
			  totSize += sizeof(AliHLTPHOSDigitDataStruct);
			  currentchannel = fShmPtr->NextChannel();
			}
		      else // Oh well, better use the low gain channel then
			{
			  AddDigit(currentchannelLG, coord1);
			  j++;
			  totSize += sizeof(AliHLTPHOSDigitDataStruct);
			  currentchannel = fShmPtr->NextChannel();
			}
		    }
		  else // No available high gain channel for this crystal, adding the low gain one
		    {
		      AddDigit(currentchannelLG, coord1);
		      j++;
		      totSize += sizeof(AliHLTPHOSDigitDataStruct);
		    }
		}
	      else //Fine, no more channels, better add this one...
		{
		  AddDigit(currentchannelLG, coord1);
		  j++;
		  totSize += sizeof(AliHLTPHOSDigitDataStruct);
		}
	    } 
	  else // Cool, no annoying low gain channel for this channel
	    {
	      AddDigit(currentchannel, coord1);
	      j++;
	      currentchannel = fShmPtr->NextChannel();
	    }
	}
    }

  fDigitCount += j;
  return fDigitCount; 
}

void 
AliHLTPHOSDigitMaker::SetGlobalHighGainFactor(Float_t factor)
{
  //See header file for documentation
  for(int x = 0; x < NXCOLUMNSMOD; x++)
    {
      for(int z = 0; z < NZROWSMOD; z++)
	{
	  fHighGainFactors[x][z] = factor;
	}
    }
}

void
AliHLTPHOSDigitMaker::SetGlobalLowGainFactor(Float_t factor)
{
  //See header file for documentation
  for(int x = 0; x < NXCOLUMNSMOD; x++)
    {
      for(int z = 0; z < NZROWSMOD; z++)
	{
	  fLowGainFactors[x][z] = factor;
	}
    }
}

void
AliHLTPHOSDigitMaker::SetBadChannelMask(TH2F* badChannelHGHist, TH2F* badChannelLGHist, Float_t qCut)
{
 for(int x = 0; x < NXCOLUMNSMOD; x++)
    {
      for(int z = 0; z < NZROWSMOD; z++)
	{
	  if(badChannelHGHist->GetBinContent(x, z) < qCut && badChannelHGHist->GetBinContent(x, z) > 0)
	    {
	      fBadChannelMask[x][z][HIGHGAIN] = 1;
	    }
	  else
	    {
	      fBadChannelMask[x][z][HIGHGAIN] = 0;
	    }
	  if(badChannelLGHist->GetBinContent(x, z) < qCut && badChannelLGHist->GetBinContent(x, z) > 0)
	    {
	      fBadChannelMask[x][z][LOWGAIN] = 0;
	    }
	  else
	    {
	      fBadChannelMask[x][z][LOWGAIN] = 0;
	    }
	}
    }
}

