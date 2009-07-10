/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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

/* $Id: $ */

// This class extracts the signal parameters (energy, time, quality)
// from ALTRO samples. Energy is in ADC counts, time is in time bin units.
// A coarse algorithm takes the energy as the maximum
// sample, time is the first sample index, and the quality is the number
// of bunches in the signal.
// 
//     AliPHOSRawFitterv0 *fitterv0=new AliPHOSRawFitterv0();
//     fitterv0->SetSamples(sig,sigStart,sigLength);
//     fitterv0->SetNBunches(nBunches);
//     fitterv0->SetChannelGeo(module,cellX,cellZ,caloFlag);
//     fitterv0->SetCalibData(fgCalibData) ;
//     fitterv0->Eval();
//     Double_t amplitude = fitterv0.GetEnergy();
//     Double_t time      = fitterv0.GetTime();
//     Bool_t   isLowGain = fitterv0.GetCaloFlag()==0;

// Author: Yuri Kharlov

// --- ROOT system ---
#include "TArrayI.h"
#include "TMath.h"
#include "TObject.h"

// --- AliRoot header files ---
#include "AliPHOSRawFitterv0.h"
#include "AliPHOSCalibData.h"
#include "AliLog.h"

ClassImp(AliPHOSRawFitterv0)

//-----------------------------------------------------------------------------
AliPHOSRawFitterv0::AliPHOSRawFitterv0():
  TObject(),
  fSignal(0),
  fModule(0),
  fCellX(0),
  fCellZ(0),
  fCaloFlag(0),
  fStart(0),
  fLength(0),
  fNBunches(0),
  fPedSubtract(kFALSE),
  fEnergy(-111),
  fTime(-111),
  fQuality(0.),
  fPedestalRMS(0.),
  fAmpOffset(0),
  fAmpThreshold(0),
  fOverflow(kFALSE),
  fCalibData(0)
{
  //Default constructor
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv0::~AliPHOSRawFitterv0()
{
  //Destructor
  delete [] fSignal;
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv0::AliPHOSRawFitterv0(const AliPHOSRawFitterv0 &phosFitter ):
  TObject(),
  fSignal      (phosFitter.fSignal),
  fModule      (phosFitter.fModule),
  fCellX       (phosFitter.fCellX),
  fCellZ       (phosFitter.fCellZ),
  fCaloFlag    (phosFitter.fCaloFlag),
  fStart       (phosFitter.fStart),
  fLength      (phosFitter.fLength),
  fNBunches    (phosFitter.fNBunches),
  fPedSubtract (phosFitter.fPedSubtract),
  fEnergy      (phosFitter.fEnergy),
  fTime        (phosFitter.fTime),
  fQuality     (phosFitter.fQuality),
  fPedestalRMS (phosFitter.fPedestalRMS),
  fAmpOffset   (phosFitter.fAmpOffset),
  fAmpThreshold(phosFitter.fAmpThreshold),
  fOverflow    (phosFitter.fOverflow),
  fCalibData   (phosFitter.fCalibData)
{
  //Copy constructor
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv0& AliPHOSRawFitterv0::operator = (const AliPHOSRawFitterv0 &phosFitter)
{
  //Assignment operator.

  if(this != &phosFitter) {
    fSignal       = phosFitter.fSignal;
    fModule       = phosFitter.fModule;
    fCellX        = phosFitter.fCellX;
    fCellZ        = phosFitter.fCellZ;
    fCaloFlag     = phosFitter.fCaloFlag;
    fStart        = phosFitter.fStart;
    fLength       = phosFitter.fLength;
    fNBunches     = phosFitter.fNBunches;
    fPedSubtract  = phosFitter.fPedSubtract;
    fEnergy       = phosFitter.fEnergy;
    fTime         = phosFitter.fTime;
    fQuality      = phosFitter.fQuality;
    fPedestalRMS  = phosFitter.fPedestalRMS;
    fAmpOffset    = phosFitter.fAmpOffset;
    fAmpThreshold = phosFitter.fAmpThreshold;
    fOverflow     = phosFitter.fOverflow;
    fCalibData    = phosFitter.fCalibData;
  }

  return *this;
}

//-----------------------------------------------------------------------------

void AliPHOSRawFitterv0::SetSamples(const UShort_t *sig, Int_t sigStart, Int_t sigLength)
{
  // Set the sample array, its start and length in time bin units

  fStart   = sigStart;
  fLength  = sigLength;
  fSignal  = new UShort_t[fLength];
  for (Int_t i=0; i<fLength; i++) {
    fSignal[i] = sig[i];
  }
}
//-----------------------------------------------------------------------------

void AliPHOSRawFitterv0::SetChannelGeo(const Int_t module, const Int_t cellX,
				     const Int_t cellZ,  const Int_t caloFlag)
{
  // Set geometry address of the channel
  // (for a case if fitting parameters are different for different channels)
  
  fModule   = module;
  fCellX    = cellX;
  fCellZ    = cellZ;
  fCaloFlag = caloFlag;
}
//-----------------------------------------------------------------------------

Bool_t AliPHOSRawFitterv0::Eval()
{
  // Calculate signal parameters (energy, time, quality) from array of samples
  // Energy is a maximum sample minus pedestal 9
  // Time is the first time bin
  // Signal overflows is there are at least 3 samples of the same amplitude above 900

  fEnergy  = 0;
  if (fNBunches > 1) {
    fQuality = 1000;
    return kTRUE;
  }
  
  const Float_t kBaseLine   = 1.0;
  const Int_t   kPreSamples = 10;

  Float_t  pedMean   = 0;
  Float_t  pedRMS    = 0;
  Int_t    nPed      = 0;
  UShort_t maxSample = 0;
  Int_t    nMax      = 0;

  for (Int_t i=0; i<fLength; i++) {
    if (i<kPreSamples) {
      nPed++;
      pedMean += fSignal[i];
      pedRMS  += fSignal[i]*fSignal[i] ;
    }
    
    if(fSignal[i] > maxSample) maxSample = fSignal[i];
    if(fSignal[i] == maxSample) nMax++;

    if(fPedSubtract) {
      if( (fSignal[i]-(Float_t)(pedMean/nPed)) >kBaseLine ) fTime = (Double_t)i;
    }
    else //ZS
      if( (fSignal[i]-(Float_t)fAmpOffset) >kBaseLine ) fTime = (Double_t)i;
  }
  
  fEnergy = (Double_t)maxSample;
  if (maxSample > 900 && nMax > 2) fOverflow = kTRUE;

  if (fPedSubtract) {
    if (nPed > 0) {
      fPedestalRMS=(pedRMS - pedMean*pedMean/nPed)/nPed ;
      if(fPedestalRMS > 0.) 
	fPedestalRMS = TMath::Sqrt(fPedestalRMS) ;
      fEnergy -= (Double_t)(pedMean/nPed); // pedestal subtraction
    }
    else
      return kFALSE;
  }
  else {
    //take pedestals from DB
    Double_t pedestal = (Double_t) fAmpOffset ;
    if (fCalibData) {
      Float_t truePed       = fCalibData->GetADCpedestalEmc(fModule, fCellZ, fCellX) ;
      Int_t   altroSettings = fCalibData->GetAltroOffsetEmc(fModule, fCellZ, fCellX) ;
      pedestal += truePed - altroSettings ;
    }
    else{
      AliWarning(Form("Can not read data from OCDB")) ;
    }
    fEnergy-=pedestal ;
  }
  if (fEnergy < kBaseLine) fEnergy = 0;
  
  return kTRUE;

}
