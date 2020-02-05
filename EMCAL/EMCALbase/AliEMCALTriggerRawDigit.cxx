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

// --- ROOT system ---
#include <iomanip>
#include <iostream>
#include <TMath.h>

#include "AliEMCALTriggerRawDigit.h"

ClassImp(AliEMCALTriggerRawDigit) ;

AliEMCALTriggerRawDigit::AliEMCALTriggerRawDigit() : AliEMCALRawDigit(),
fTriggerBits(0),
fNL0Times(0),
fL0Times(),
fL1TimeSum(-1),
fL1SubRegion(-1)
{
  for (Int_t i = 0; i < 10; i++) fL0Times[i] = -1;
}

AliEMCALTriggerRawDigit::AliEMCALTriggerRawDigit(Int_t id, Int_t timeSamples[], Int_t nSamples) 
: AliEMCALRawDigit(id, timeSamples, nSamples),
fTriggerBits(0),
fNL0Times(0),
fL0Times(),
fL1TimeSum(-1),
fL1SubRegion(-1)
{	
  for (Int_t i = 0; i < 10; i++) fL0Times[i] = -1;
}

AliEMCALTriggerRawDigit::~AliEMCALTriggerRawDigit() 
{	
  //delete [] fL0Times;
}

Bool_t AliEMCALTriggerRawDigit::SetL0Time(const Int_t i)
{  
  if (fNL0Times > 9)
  {
    AliError("More than 10 L0 times!");
    return kFALSE;
  }

  for (Int_t j = 0; j < fNL0Times; j++)
  {
    if (i == fL0Times[j]) 
    {
      AliDebug(1,Form("Digit id %d: L0 time %d already there! Won't add it twice",fId,i));
      return kFALSE;
    }
  }
  
  fNL0Times++;

  fL0Times[fNL0Times - 1] = i;
  
  return kTRUE;
}

Bool_t AliEMCALTriggerRawDigit::GetL0Time(const Int_t i, Int_t& time) const
{
	
	if (i < 0 || i > fNL0Times)
	{
		AliError("Bad index!");
		return kFALSE;
	}
		
	time = fL0Times[i];
	
	return kTRUE;
}

Bool_t AliEMCALTriggerRawDigit::GetL0Times(Int_t times[]) const
{	
  for (Int_t i = 0; i < fNL0Times; i++) times[i] = fL0Times[i];
  
  return kTRUE;
}

Int_t AliEMCALTriggerRawDigit::GetL0TimeSum(const Int_t time) const
{	
  Int_t value = 0;
  
  for (Int_t i = 0; i < fNSamples; i++)
  {
    Int_t timeBin, amp;
    GetTimeSample(i, timeBin, amp);
    
    if (timeBin >= time && timeBin < time + 4) value += amp;
  }
  
  return value;
}

Int_t AliEMCALTriggerRawDigit::GetTriggerBit(const TriggerType_t type, const Int_t mode) const
{	
  Int_t shift = kTriggerTypeEnd * mode;
  Int_t mask  = 1 << type;
  
  return ((fTriggerBits >> shift) & mask);
}	

void AliEMCALTriggerRawDigit::Print(const Option_t* /*opt*/) const
{	
  std::cout << *this;
}

void AliEMCALTriggerRawDigit::PrintStream(std::ostream &output) const {
  AliEMCALRawDigit::PrintStream(output);
  
  output << "| L0: (" << GetTriggerBit(kL0,1) << "," << GetTriggerBit(kL0,0) << ") / " << fNL0Times << " Time(s): \n";
  for (Int_t i = 0; i < fNL0Times; i++) 
  {
    Int_t time;
    if (GetL0Time(i, time)) output << "| "<< time << " ";
  }
  output << "\n";
  
  output << "| L1: g high (" << GetTriggerBit(kL1GammaHigh,1) << "," << GetTriggerBit(kL1GammaHigh,0) 
         << ") g low (" << GetTriggerBit(kL1GammaLow,1) << "," << GetTriggerBit(kL1GammaLow,0) 
         << ") j high (" << GetTriggerBit(kL1JetHigh,1) << "," << GetTriggerBit(kL1JetHigh,0) 
         << ") j low (" << GetTriggerBit(kL1JetLow,1) << "," << GetTriggerBit(kL1JetLow,0) 
         << ") / Time sum: " << fL1TimeSum << "\n";
}

std::ostream &operator<<(std::ostream &in, const AliEMCALTriggerRawDigit &dig) {
  in << dig;
  return in;
}
