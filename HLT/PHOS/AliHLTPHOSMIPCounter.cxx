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

#include "AliHLTPHOSMIPCounter.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "TH2I.h"

AliHLTPHOSMIPCounter::AliHLTPHOSMIPCounter()
    : AliHLTPHOSBase(),
    fMIPCountEvent ( 0 ),
    fMIPCountTotal ( 0 ),
    fMIPRate ( 0 ),
    fLowerBound ( 0 ),
    fUpperBound ( 0 ),
    fUpperStartTime ( 0 ),
    fLowerStartTime ( 0 ),
    fZeroThreshold ( 0 ),
    fChannelHistPtr ( 0 )
{ 

}


AliHLTPHOSMIPCounter::~AliHLTPHOSMIPCounter()
{
}

Int_t
AliHLTPHOSMIPCounter::CountMIPs(AliHLTPHOSDigitContainerDataStruct* digitContainerPtr)
{
  fMIPCountEvent = 0;
  Bool_t IsMIP = true;
  Int_t *dataPtr = 0;
  AliHLTPHOSDigitDataStruct *digitPtr;
  for(UInt_t i = 0; i < digitContainerPtr->fNDigits; i++)
  {
    digitPtr = &(digitContainerPtr->fDigitDataStruct[i]);
    dataPtr = digitPtr->fData;
    if(digitPtr->fCrazyness != 0)
    {
      continue;
    }
    if(digitPtr->fAmplitude < fLowerBound || digitPtr->fAmplitude > fUpperBound)
    {
      continue;
    }
    for(Int_t time = (Int_t)(digitPtr->fTime - 2); time < (digitPtr->fTime - 3); time++)
    {
      if((Float_t)dataPtr[time] < (digitPtr->fAmplitude - (digitPtr->fAmplitude)/10))
      {
	IsMIP = false;
	break;
      }
    }
    if(!IsMIP)
      continue;
    for(Int_t sample = 0; sample < fLowerStartTime; sample++)
      {
	if(dataPtr[sample] > fZeroThreshold || dataPtr[sample] < -fZeroThreshold)
	{
	  IsMIP = false;
	  break;
	}
      }
      if(dataPtr[(Int_t)fUpperStartTime + 3] < fZeroThreshold)
	IsMIP = false;
    if(IsMIP)
    {
      fMIPCountEvent++;
      fChannelHistPtr->Fill(digitPtr->fX, digitPtr->fZ);
    }
  }
  fMIPCountTotal += fMIPCountEvent;
  return fMIPCountEvent; 
}
