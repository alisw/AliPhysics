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



#include "AliHLTPHOSSharedMemoryInterface.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSValidCellDataStruct.h"
#include <iostream>


AliHLTPHOSSharedMemoryInterface::AliHLTPHOSSharedMemoryInterface(): fCurrentChannel(0), 
								    fCellEnergiesPtr(0), 
								    fIsSetMemory(false), 
								    fMaxCnt(0),
								    fCurrentX(0),
								    fCurrentZ(0),
								    fCurrentGain(0),
								    fCurrentCnt(0), 
								    fCharDataOffset(0), 
								    fCharPtr(0), 
								    fIntPtr(0)
								  
 
{
 fCharDataOffset = sizeof(AliHLTPHOSRcuCellEnergyDataStruct);
}


AliHLTPHOSSharedMemoryInterface::~AliHLTPHOSSharedMemoryInterface()
{

}


AliHLTPHOSValidCellDataStruct*   
AliHLTPHOSSharedMemoryInterface::NextChannel()
{
  // Returns the next channel of the current AliHLTPHOSRcuCellEnergyDataStruct
  // Returns zero when all cannels are read
  /*
  if(fCurrentCnt < fMaxCnt)
    {
     fCurrentChannel = &fCellEnergiesPtr->fValidData[fCurrentCnt];
     fCurrentChannel->fData = fIntPtr; 
     fIntPtr +=  fCurrentChannel->fNSamples;
     fCurrentCnt ++;
     return fCurrentChannel;
    }
  */

  // Changed by OD
  if(fCurrentCnt < fMaxCnt)
    {
      for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
	{
	  for(Int_t z = 0; z < N_ZROWS_MOD; z++)
	    {
	      for(Int_t gain = 0; gain < N_GAINS; gain++)
		{
		  fCurrentChannel =  &(fCellEnergiesPtr->fValidData[x][z][gain]);
		  if(fCurrentChannel->fID == fCurrentCnt)
		    {
		      fCurrentChannel->fData = fIntPtr; 
		      fIntPtr +=  fCurrentChannel->fNSamples;
		      fCurrentCnt ++;
		      return fCurrentChannel;
		    }
     		}
	    }
	}
    }
  else
    {
      Reset();
      return 0;
    }
}


void
AliHLTPHOSSharedMemoryInterface::SetMemory(AliHLTPHOSRcuCellEnergyDataStruct *rcuCellEnergyPtr)
{
  //Shutting up rule checker
  fCellEnergiesPtr =  rcuCellEnergyPtr;
  fMaxCnt =  fCellEnergiesPtr->fCnt;
  PingPongPointer();
  fIsSetMemory = true;
}


void
AliHLTPHOSSharedMemoryInterface::Reset()
{
  //Shutting up rule checker
  fMaxCnt =0;
  fCurrentCnt = 0;
  fCurrentX = 0;
  fCurrentZ = 0;
  fCurrentGain = 0;
  fIsSetMemory = false;
}


void 
AliHLTPHOSSharedMemoryInterface::PingPongPointer()
{
  // ping pong ping ping pong ping pong
  fCharPtr = (char *)fCellEnergiesPtr ;
  fCharPtr += fCharDataOffset; 
  fIntPtr = (Int_t *)fCharPtr; 
}
