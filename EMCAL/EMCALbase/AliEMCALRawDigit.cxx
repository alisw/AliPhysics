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
#include <Riostream.h>
#include <TMath.h>

// --- AliRoot header files ---
#include "AliEMCALRawDigit.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliEMCALRawDigit) ;
/// \endcond

///
/// Default constructor 
//____________________________________________________________________________
AliEMCALRawDigit::AliEMCALRawDigit() : TObject(),
fId(-1),
fNSamples(0),
fSamples(0x0),
fAmplitude(0),
fTime(0)
{ }

///
/// Constructor 
//____________________________________________________________________________
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

///
/// Destructor, delete array of time samples
//____________________________________________________________________________
AliEMCALRawDigit::~AliEMCALRawDigit() 
{
  if(fSamples) delete [] fSamples;
  fSamples = NULL;
}

///
/// Clear, delete array of time samples
//____________________________________________________________________________
void AliEMCALRawDigit::Clear(Option_t *) 
{
  if(fSamples) delete [] fSamples;
  fSamples = NULL;
}

///
/// \return the time and amplitude of a given time sample and if the sample was ok
//____________________________________________________________________________
Bool_t AliEMCALRawDigit::GetTimeSample(const Int_t iSample, Int_t& timeBin, Int_t& amp) const
{  
  if (iSample > fNSamples || iSample < 0) return kFALSE;
  
  amp     = (Short_t)(fSamples[iSample] & 0xFFFF);
  timeBin = (Short_t)(fSamples[iSample] >> 16 );
  
  return kTRUE;
}

///
/// Sets the time samples
//____________________________________________________________________________
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

///
/// Checks the maximum amplitude in the time sample
//____________________________________________________________________________
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

///
/// Compares two digits with respect to its Id
/// to sort according increasing Id
//____________________________________________________________________________
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

///
/// Print digit stored info
//____________________________________________________________________________
void AliEMCALRawDigit::Print(const Option_t* /*opt*/) const
{  
  printf("===\n| Digit id: %4d / %d Time Samples: \n",fId,fNSamples);
  
  for (Int_t i=0; i < fNSamples; i++) 
  {
    Int_t timeBin=-1, amp=0;
    GetTimeSample(i, timeBin, amp);
    printf("| (%d,%d) ",timeBin,amp);
  }
  
  printf("\n");
}
