/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*



 Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

// --- ROOT system ---AliEMCALTriggerRawDigit
#include <Riostream.h>
#include <TMath.h>

#include "AliEMCALTriggerRawDigit.h"
#include "AliLog.h"

ClassImp(AliEMCALTriggerRawDigit)

//____________________________________________________________________________
AliEMCALTriggerRawDigit::AliEMCALTriggerRawDigit() : AliEMCALRawDigit(),
fL0Trigger(0),
fNL0Times(0),
fL0Times(),
fL1TimeSum(-1)
{
	// default ctor 
	for (Int_t i = 0; i < 10; i++) fL0Times[i] = -1;
}

//____________________________________________________________________________
AliEMCALTriggerRawDigit::AliEMCALTriggerRawDigit(Int_t id, Int_t timeSamples[], Int_t nSamples) : AliEMCALRawDigit(id, timeSamples, nSamples),
fL0Trigger(0),
fNL0Times(0),
fL0Times(),
fL1TimeSum(-1)
{
	//
	for (Int_t i = 0; i < 10; i++) fL0Times[i] = -1;
}

//____________________________________________________________________________
AliEMCALTriggerRawDigit::~AliEMCALTriggerRawDigit() 
{
	//
	//delete [] fL0Times;
}

//____________________________________________________________________________
Bool_t AliEMCALTriggerRawDigit::SetL0Time(const Int_t i)
{
	//
	for (Int_t j = 0; j < fNL0Times; j++)
	{
		if (i == fL0Times[j]) 
		{
			AliWarning("L0 time already there! Won't add it twice");
			return kFALSE;
		}
	}
	
	fNL0Times++;
	
	if (fNL0Times > 9)
	{
		AliError("More than 10 L0 times!");
		return kFALSE;
	}
	
	fL0Times[fNL0Times - 1] = i;
	
	return kTRUE;
}

//____________________________________________________________________________
Bool_t AliEMCALTriggerRawDigit::GetL0Time(const Int_t i, Int_t& time) const
{
	//
	if (i < 0 || i > fNL0Times)
	{
		AliError("Bad index!");
		return kFALSE;
	}
		
	time = fL0Times[i];
	
	return kTRUE;
}

//____________________________________________________________________________
Bool_t AliEMCALTriggerRawDigit::GetL0Times(Int_t times[]) const
{
	//
	if (sizeof(times) < (sizeof(Int_t) * fNL0Times))
	{
		AliError("Array size too small!");
		return kFALSE;
	}
	
	for (Int_t i = 0; i < fNL0Times; i++) times[i] = fL0Times[i];
	
	return kTRUE;
}

//____________________________________________________________________________
Int_t AliEMCALTriggerRawDigit::GetL0TimeSum(const Int_t time) const
{
	//
	
	Int_t value = 0;
	
	for (Int_t i = 0; i < fNSamples; i++)
	{
		Int_t timeBin, amp;
		GetTimeSample(i, timeBin, amp);
		
		if (timeBin >= time && timeBin < time + 4) value += amp;
	}
	
	return value;
}

//____________________________________________________________________________
void AliEMCALTriggerRawDigit::Print(const Option_t* /*opt*/) const
{
	//
	printf("===\n| Digit id: %4d / %d Time Samples: \n",fId,fNSamples);
	for (Int_t i=0; i < fNSamples; i++) 
	{
		Int_t timeBin, amp;
		GetTimeSample(i, timeBin, amp);
		printf("| (%d,%d) ",timeBin,amp);
	}	
	printf("\n");
	printf("| L0: %4d / %d Time(s): \n",fL0Trigger,fNL0Times);
	for (Int_t i = 0; i < fNL0Times; i++) 
	{
		Int_t time;
		if (GetL0Time(i, time)) printf("| %d ",time);
	}
	printf("\n");
	printf("| Time sum: %d\n", fL1TimeSum);
}

