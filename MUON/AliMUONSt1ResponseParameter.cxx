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

/* $Id$ */

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONSt1ResponseParameter
// ---------------------------------
// Describes a set of filters to be applied to a digital value
// in order to simulate electronics characteristics 
// (pedestal, noise, sticky bits, etc....)
// Threshold levels for the MANU zero supression algorithm are included.
// Included in AliRoot 2003/01/28

#include <fstream>

#include <TRandom.h>
#include <TString.h>

#include "AliMUONSt1ResponseParameter.h"
#include "AliLog.h"

ClassImp(AliMUONSt1ResponseParameter);

//_________________________________________________________________________
AliMUONSt1ResponseParameter::AliMUONSt1ResponseParameter()
  :TNamed()
{
// default constructor
  fPedestalMode = kNone;
  fNoiseMode = kNone;
  fState=1;
  fNofSigma=3;
  fStickyOn=fStickyOff=0;
}

//_________________________________________________________________________
AliMUONSt1ResponseParameter::AliMUONSt1ResponseParameter(const TString& name,const TString& title)
:TNamed(name,title)
{
// normal constructor
  fPedestalMode = kNone;
  fNoiseMode = kNone;
  fState=1;
  fNofSigma=3;
  fStickyOn=fStickyOff=0;
}

//_________________________________________________________________________
AliMUONSt1ResponseParameter::~AliMUONSt1ResponseParameter()
{
// destructor
}

//_________________________________________________________________________
void AliMUONSt1ResponseParameter::SetState(Bool_t state)
{
// If set to off, no information will be available from the electronics
// ---

  fState=state;
}

//_________________________________________________________________________
void AliMUONSt1ResponseParameter::SetPedestal(Double_t val)
{
// Set pedestal values to a constant 
// ---

  fPedestalMode = kValue;
  fPedestalParam.value = val;
}

//_________________________________________________________________________
void AliMUONSt1ResponseParameter::SetPedestal(Double_t mean,Double_t sigma)
{
// Set pedestal values to a parameterized gaussian
// ---

  fPedestalMode = kGauss;
  fPedestalParam.gauss.mean  = mean;
  fPedestalParam.gauss.sigma = sigma;
}
//_________________________________________________________________________
void AliMUONSt1ResponseParameter::SetPedestal(const TString& fileName)
{
// Set pedestal values to those given in a file
// ---

  ifstream file(fileName.Data());
  if (file.good()){
    fPedestalMode = kFile;
    for (Int_t ch=0;ch<fgkNofChannels;ch++) {
      Float_t value;
      file>>value;
      fPedestalParam.values[ch] = value;
      //cout<<"Pedestal["<<i<<"]["<<ch<<"]="<<value<<endl;
    }
    file.close();
  } else {
    AliWarning(Form("Can't read file %s",fileName.Data()));
    SetPedestal(150.,10.);
  }
}

//_________________________________________________________________________
void AliMUONSt1ResponseParameter::UnSetPedestal()
{
// Set pedestal values to 0.
// ---

  fPedestalMode=kNone;
}

//_________________________________________________________________________
void AliMUONSt1ResponseParameter::SetNoise(Double_t val)
{
// Set Noise values to a constant value
// ---

  fNoiseMode = kValue;
  fNoiseParam.value = val;
}

//_________________________________________________________________________
void AliMUONSt1ResponseParameter::SetNoise(Double_t mean,Double_t sigma)
{
// Set Noise values to a parameterized gaussian
// ---

  fNoiseMode = kGauss;
  fNoiseParam.gauss.mean  = mean;
  fNoiseParam.gauss.sigma = sigma;
}

//_________________________________________________________________________
void AliMUONSt1ResponseParameter::SetNoise(const TString& fileName)
{
// Set Noise values to those given in a file
// ---

  ifstream file(fileName.Data());
  if (file.good()){
    fNoiseMode = kFile;
    for (Int_t ch=0;ch<fgkNofChannels;ch++) {
      Float_t value;
      file>>value;
      fNoiseParam.values[ch] = value;
      //cout<<"Noise["<<i<<"]["<<ch<<"]="<<value<<endl;
    }
    file.close();
  } else {
    AliWarning(Form("Can't read file %s",fileName.Data()));
    SetNoise(150.,10.);
  }
}

//_________________________________________________________________________
void AliMUONSt1ResponseParameter::SetNofSigma(Int_t nofSigma)
{
// set Number of sigmas to be applied as threshold (for zero suppression)
// ---

  fNofSigma = nofSigma;
}

//_________________________________________________________________________
void AliMUONSt1ResponseParameter::SetStickyBitOn (Int_t bit,Int_t val)
{
// In the response, this bit will always be set to 1 (unless <State> is off)
// ---

  if (val)
    fStickyOn |= (1<<bit);
  else
    fStickyOn &= ~(1<<bit);
}

//_________________________________________________________________________
void AliMUONSt1ResponseParameter::SetStickyBitOff(Int_t bit,Int_t val)
{
// In the response, this bit will always be set to 0
// ---

  if (val)
    fStickyOff |= (1<<bit);
  else
    fStickyOff &= ~(1<<bit);
}

//_________________________________________________________________________
Int_t AliMUONSt1ResponseParameter::ApplyPedestal(Int_t base,Int_t GC) const
{
// calculate the response to <base>, with respect to the current pedestal
// parameters
// --
  Double_t ped = Choose(fPedestalMode,fPedestalParam,GC);
  Double_t nse = Choose(fNoiseMode,fNoiseParam,GC);
  Double_t noise     = gRandom->Gaus(0, nse);
  base+=(Int_t)(noise + ped);
  if (base-ped-noise*fNofSigma<0) base=0;

  return base;
}
//_________________________________________________________________________
Int_t AliMUONSt1ResponseParameter::ApplyStickyBits(Int_t base) const
{
// set the response to <base>, with respect to the current stickyBits
// parameters
// --
  base |= fStickyOn;
  base &= (~fStickyOff);
  return base;
}
//////////////////// Privates methods
//_________________________________________________________________________

Double_t AliMUONSt1ResponseParameter::Choose(TMode mode,TParam param,Int_t GC) const
{
// Choose a (pedestal/noise) value to be applied following the parameterization rule
// ---

  switch (mode){
    case kNone  : return 0;
    case kValue : return param.value;
    case kGauss : return gRandom->Gaus(param.gauss.mean,param.gauss.sigma);
    case kFile  : return param.values[GC];
  }
  AliFatal("No mode is given");
  return 0;
}
