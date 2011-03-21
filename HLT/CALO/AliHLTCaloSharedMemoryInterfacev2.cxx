// $Id: AliHLTCaloSharedMemoryInterfacev2.cxx 35071 2009-09-29 05:26:09Z phille $

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliHLTCaloSharedMemoryInterfacev2.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
#include "AliHLTCaloChannelDataStruct.h"
#include "AliHLTLogging.h"
#include "AliHLTCaloMapper.h"
#include "AliHLTCaloConstants.h"

ClassImp(AliHLTCaloSharedMemoryInterfacev2);

AliHLTCaloSharedMemoryInterfacev2::AliHLTCaloSharedMemoryInterfacev2(TString det): AliHLTCaloConstantsHandler(det),
									fCurrentChannel(0),
									fChannelDataPtr(0),
									fIsSetMemory(false),
									fHasRawData(false),
									fMaxCnt(0),
									fCurrentCnt(0), 
									fRawDataPtr(0),
									fRawData()
									//	fSpecification(0)
{
   // See header file for class documentation
  for(Int_t i = 0; i < 32; i++)
    {
      fMapperPtr[i] = 0;
    }
}

AliHLTCaloSharedMemoryInterfacev2::~AliHLTCaloSharedMemoryInterfacev2()
{

}

AliHLTCaloChannelDataStruct*
AliHLTCaloSharedMemoryInterfacev2::NextChannel()
{
  // Comment
  AliHLTCaloChannelDataStruct* tmpChannelPtr = 0;
  if(fCurrentCnt < fMaxCnt)
    {
      tmpChannelPtr = reinterpret_cast<AliHLTCaloChannelDataStruct*>(fChannelDataPtr);
      fCurrentCnt++;
      fChannelDataPtr += sizeof(AliHLTCaloChannelDataStruct);
      
     if(fHasRawData == true)
      if( false )	
	{
	  fRawData.fEnergy = tmpChannelPtr->fEnergy;
	  fRawData.fTime = tmpChannelPtr->fTime;
	  fRawData.fChannelID = tmpChannelPtr->fChannelID;
	  fRawData.fCrazyness = tmpChannelPtr->fCrazyness;
	  Reset(fRawData);
	  //AliHLTCaloMapper::ChannelId2Coordinate(const UShort_t channelId, AliHLTCaloCoordinate &channelCoord)
 
 
	  AliHLTCaloMapper::ChannelId2Coordinate( fRawData.fChannelID, fRawData.fCoordinate);
	  
	  if( fRawData.fChannelID == fRawDataPtr[0] )
	    {
	      Reset(fRawData);
	      // cout << __FILE__ << __LINE__ << "fRawData.fChannelID == fRawDataPtr[0] = " << fRawDataPtr[0] << endl;
	      // cout << " copying raw dat not yet implemnted " << endl;
	      UShort_t tmpTotSize = fRawDataPtr[1];
	      UShort_t tmpStartBin = fRawDataPtr[2];
	      UShort_t tmpBunchtSize = fRawDataPtr[3];
	      // fRawDataPtr
	      //     UShort_t tmpSamplesLeft = tmpTotSize -4;
	      int tmpSamplesLeft = tmpTotSize -4;

	      fRawData.nSamplesUsed = tmpTotSize + tmpStartBin;

	      if(tmpSamplesLeft > 0 )
		{
		  while(tmpSamplesLeft > 0)
		    {
		      for(int i=0; i < tmpBunchtSize; i++ )
			{
//			  fRawData.fDataPtr[i + tmpStartBin] = fRawDataPtr[ i+ 4];
			  tmpSamplesLeft --;
			}
		    }
		  fRawDataPtr+= tmpTotSize;
		}
	    }
	  else
	    {
	      // cout << __FILE__ << __LINE__ << "ERROR! fRawData.fChannelID = "<< fRawData.fChannelID << " but fRawDataPtr[0] = " << fRawDataPtr[0] << endl;
	    }
 
 
	  // HLTDebug("Raw data interface not yet implemented, ignoring raw data");
	}
      return tmpChannelPtr;
    }
  else
    {
      Reset();
      return 0;
    }
  return 0;
}

void  
AliHLTCaloSharedMemoryInterfacev2::NextRawChannel( )
{
  if(fHasRawData == false )
    {
      cout << __FILE__ << __LINE__<< "ERROR: no raw data present" << endl;
    }
  else
    {
      for(int i = 0; i <  200 ; i++ )

	{
	  cout << fRawDataPtr[i] << "\t";
	  if(i%16 == 0)
	    {
	      cout << endl;
	    }
	}
    }
}

void
//AliHLTCaloSharedMemoryInterfacev2::SetMemory(AliHLTCaloChannelDataHeaderStruct* channelDataHeaderPtr,  const unsigned long specification)
AliHLTCaloSharedMemoryInterfacev2::SetMemory(AliHLTCaloChannelDataHeaderStruct* channelDataHeaderPtr)
{
  //  fSpecification = specification;

  //Shutting up rule checker
  fHasRawData = channelDataHeaderPtr->fHasRawData; 
  fMaxCnt = channelDataHeaderPtr->fNChannels;
  fChannelDataPtr = reinterpret_cast<AliHLTUInt8_t*>(channelDataHeaderPtr) + sizeof(AliHLTCaloChannelDataHeaderStruct);
  

  if(fHasRawData == true)
    {
      fRawDataPtr = reinterpret_cast<  UShort_t* >(channelDataHeaderPtr); 
      int inc =  sizeof (AliHLTCaloChannelDataHeaderStruct) +  fMaxCnt*sizeof(AliHLTCaloChannelDataStruct);
      fRawDataPtr += inc/sizeof(UShort_t );
    }

    fIsSetMemory = true;
}


void
AliHLTCaloSharedMemoryInterfacev2::Reset()
{
  //Shutting up rule checker
  fCurrentCnt = 0;
  fIsSetMemory = false;
  fHasRawData = false;
}

