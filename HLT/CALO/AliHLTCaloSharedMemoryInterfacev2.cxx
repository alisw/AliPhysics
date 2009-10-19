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


AliHLTCaloSharedMemoryInterfacev2::AliHLTCaloSharedMemoryInterfacev2(): fCurrentChannel(0),
									fChannelDataPtr(0),
									fIsSetMemory(false),
									fHasRawData(false),
									fMaxCnt(0),
									fCurrentCnt(0), 
									fRawDataPtr(0),
									fRawData()
									//	fSpecification(0)
{
  //  GetSpecFromDDLIndex
  //  AliHLTCaloMapper  *fMapperPtr[32];

  /*
  for(int i=0; i < 32; i++)
    {
      fMapperPtr[i] = new AliHLTCaloMapper( AliHLTCaloMapper::GetSpecFromDDLIndex(i) ) ;
    }
  */
}



AliHLTCaloSharedMemoryInterfacev2::~AliHLTCaloSharedMemoryInterfacev2()
{

}


/*
struct AliHLTCaloChannelDataStruct
{
  Float_t fEnergy;
  Float_t fTime;
  UShort_t fChannelID;
  Short_t fCrazyness;
  //  Short_t fRawDataSize; //the size of the raw data
};
*/
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
      //     if( false )	
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
			  fRawData.fDataPtr[i + tmpStartBin] = fRawDataPtr[ i+ 4];
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

/*
  AliHLTCaloChannelDataStruct*   
  AliHLTCaloSharedMemoryInterfacev2::NextChannel()
  {
  cout << __FILE__ << __LINE__ << " TP0" << endl;

  AliHLTCaloChannelDataStruct* tmpChannelPtr = 0;
  if(fCurrentCnt < fMaxCnt)
  {
  cout << __FILE__ << __LINE__ << " TP1" << endl;  
  tmpChannelPtr = reinterpret_cast<AliHLTCaloChannelDataStruct*>(fChannelDataPtr);
  fCurrentCnt++;
  fChannelDataPtr += sizeof(AliHLTCaloChannelDataStruct);
   
  //    if(fHasRawData == true)
  if(false)	
  {
  cout << __FILE__ << __LINE__ << " TP2" << endl; 
  fRawData.fEnergy = tmpChannelPtr->fEnergy;
  fRawData.fTime = tmpChannelPtr->fTime;
  fRawData.fChannelID = tmpChannelPtr->fChannelID; 
  fRawData.fCrazyness  = tmpChannelPtr->fCrazyness; 
  Reset(fRawData);
  AliHLTCaloMapper::ChannelId2Coordinate( fRawData.fChannelID, fRawData.fCoordinate);
  cout << __FILE__ << __LINE__ << " TP3" << endl; 
  if( fRawData.fChannelID ==  fRawDataPtr[0]  )
  {
  cout << __FILE__ << __LINE__ << " TP4" << endl; 
  Reset(fRawData);
  cout << __FILE__ << __LINE__ << " TP5" << endl; 
  UShort_t tmpTotSize =  fRawDataPtr[1];
  UShort_t tmpStartBin   =  fRawDataPtr[2];
  UShort_t tmpBunchtSize =  fRawDataPtr[3];
  UShort_t tmpSamplesLeft = tmpTotSize -4; 
  fRawData.nSamplesUsed =  tmpTotSize +  tmpStartBin;
  cout << __FILE__ << __LINE__ << " TP6" << endl; 

  while(tmpSamplesLeft > 0)
  {
  if ( (tmpBunchtSize <  ( ALTROMAXSAMPLES +4)) && tmpBunchtSize > 0   )
  {
  cout << __FILE__ << __LINE__ << " TP7" << endl; 
  for(int i=0; i < tmpBunchtSize; i++ )
  {
  cout << __FILE__ << __LINE__ << "  TMP tot size =" <<  tmpBunchtSize <<  tmpTotSize  << endl;
  cout << __FILE__ << __LINE__ << "  TMP bunch size =" <<  tmpTotSize  << endl;
  cout << __FILE__ << __LINE__ << " TP8, i=" << i <<", tmpStartBin = "<< tmpStartBin  << " fCurrentCnt = "<< fCurrentCnt << endl; 
  cout << __FILE__ << __LINE__ << " acessing fRawData.fDataPtr["<<  i + tmpStartBin  << "]" << endl; 
  fRawData.fDataPtr[ i + tmpStartBin ] = fRawDataPtr[ i+ 4];
  //			    cout << __FILE__ << __LINE__ << " acessing fRawData.fDataPtr["<<  i + tmpStartBin  << "]" << endl; 
  tmpSamplesLeft --;
  }
  fRawDataPtr+= tmpTotSize;
  }
  else
  {
  tmpSamplesLeft = -1;
  }
  }
  }
  }
  return tmpChannelPtr;
  } 
  else
  {
  Reset();
  return 0;
  }
  return 0;
  //   }
  }
*/


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

 
void 
AliHLTCaloSharedMemoryInterfacev2::Reset(AliHLTCaloChannelRawDataStruct &str)
{
  for(int i=0; i< ALTROMAXSAMPLES; i++ )
    {
      str.fDataPtr[i] = 0;
    }
 
}
