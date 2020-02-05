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

// --- AliRoot header files ---
#include "AliEMCALRawDigit.h"
#include "AliLog.h"

ClassImp(AliEMCALRawDigit) ;

AliEMCALRawDigit::AliEMCALRawDigit() : TObject(),
fId(-1),
fNSamples(0),
fSamples(0x0),
fAmplitude(0),
fTime(0)
{ }

AliEMCALRawDigit::AliEMCALRawDigit(Int_t id, Int_t timeSamples[], Int_t nSamples) : TObject(),
fId(id),
fNSamples(nSamples),
fSamples(0x0),
fAmplitude(0),
fTime(0)
{
  fSamples = new Int_t[fNSamples];
  for (Int_t i = 0; i < fNSamples; i++) fSamples[i] = timeSamples[i];
}

AliEMCALRawDigit::~AliEMCALRawDigit() 
{
  if(fSamples) delete [] fSamples;
  fSamples = NULL;
}

void AliEMCALRawDigit::Clear(Option_t *) 
{
  if(fSamples) delete [] fSamples;
  fSamples = NULL;
}

Bool_t AliEMCALRawDigit::operator==(const AliEMCALRawDigit &rhs) const {
  return Compare(&rhs) == 0;
}

Bool_t AliEMCALRawDigit::operator<(const AliEMCALRawDigit &rhs) const {
  return Compare(&rhs) < 0;
}

Bool_t AliEMCALRawDigit::GetSamples(Int_t samples[], Int_t ns) const
{
  if (ns <= 0) return kFALSE;
  
  int smax=TMath::Min(ns,fNSamples);
  for (int i=0;i<smax;i++) samples[i]=fSamples[i];
  return kTRUE;
}

Bool_t AliEMCALRawDigit::GetTimeSample(const Int_t iSample, Int_t& sample) const
{
  if (iSample > fNSamples || iSample < 0) return kFALSE;

  sample = fSamples[iSample];

  return kTRUE;
}

Bool_t AliEMCALRawDigit::GetTimeSample(const Int_t iSample, Int_t& timeBin, Int_t& amp) const
{  
  if (iSample > fNSamples || iSample < 0) return kFALSE;
  
  amp     = (Short_t)(fSamples[iSample] & 0xFFFF);
  timeBin = (Short_t)(fSamples[iSample] >> 16 );
  
  return kTRUE;
}

void AliEMCALRawDigit::SetTimeSamples(const Int_t timeSamples[], const Int_t nSamples) 
{  
  if (fSamples) 
  {
    AliDebug(1,"Samples already filled: delete first!");
    fNSamples = 0;
    delete [] fSamples;
  }
  
  fNSamples = nSamples;
  fSamples = new Int_t[fNSamples];
  for (Int_t i = 0; i < fNSamples; i++) fSamples[i] = timeSamples[i];
}

Bool_t AliEMCALRawDigit::GetMaximum(Int_t& amplitude, Int_t& time) const
{  
  if (!fNSamples)
  {
    AliDebug(1,"Digit has no time sample");
    return kFALSE;
  }
		
  amplitude = 0;
  for (Int_t i = 0; i < fNSamples; i++)
  {
    Int_t t, a;
    if (GetTimeSample(i, t, a))
    {
      if (a > amplitude)
      {
        amplitude = a;
        time      = t;
      }
    }
  }
  
  return kTRUE;
}

Int_t AliEMCALRawDigit::Compare(const TObject * obj) const
{	
  Int_t rv=0;
  
  AliEMCALRawDigit* digit = (AliEMCALRawDigit*)obj; 
  
  Int_t iddiff = fId - digit->GetId();
  
  if (iddiff > 0) 
    rv =  1;
  else if (iddiff < 0)
    rv = -1; 
  else
    rv =  0;
  
  return rv; 
}

void AliEMCALRawDigit::Print(const Option_t* /*opt*/) const
{  
  std::cout << *this;
}

void AliEMCALRawDigit::PrintStream(std::ostream &stream) const {
  stream << "===\n| Digit id: " << std::setw(4) << fId << " / " << fNSamples << " Time Samples: \n";
  
  for (Int_t i=0; i < fNSamples; i++) 
  {
    Int_t timeBin=-1, amp=0;
    GetTimeSample(i, timeBin, amp);
    stream << "| (" << timeBin << "," << amp << ") ";
  }
  
  stream << "\n";
} 

std::ostream &operator<<(std::ostream &stream, const AliEMCALRawDigit &dig) {
  stream << dig;
  return stream;
}