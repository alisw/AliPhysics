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



#include "AliITSCalibrationSPD.h"

const Double_t AliITSCalibrationSPD::fgkThreshDefault = 3000.;
const Double_t AliITSCalibrationSPD::fgkSigmaDefault = 250.;
const Double_t AliITSCalibrationSPD::fgkCouplColDefault = 0.;
const Double_t AliITSCalibrationSPD::fgkCouplRowDefault = 0.047;
const Double_t AliITSCalibrationSPD::fgkBiasVoltageDefault = 18.182;

ClassImp(AliITSCalibrationSPD)	
///////////////////////////////////////////////////////////////////////////
//  Calibration class for set:ITS                   
//  Specific subdetector implementation for         
//  Silicon pixels                                  
//
//  Modified by D. Elia, G.E. Bruno, H. Tydesjo
//  March-April 2006
//
///////////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliITSCalibrationSPD::AliITSCalibrationSPD():
AliITSCalibration(),
fBaseline(0.0),
fNoise(0.0),
fThresh(fgkThreshDefault),
fSigma(fgkSigmaDefault),
fCouplCol(fgkCouplColDefault),
fCouplRow(fgkCouplRowDefault),
fBiasVoltage(fgkBiasVoltageDefault),
fNrDead(0),
fDeadChannels(0),
fNrNoisy(0),
fNoisyChannels(0){
  // constructor

   SetThresholds(fgkThreshDefault,fgkSigmaDefault);
   SetCouplingParam(fgkCouplColDefault,fgkCouplRowDefault);
   SetBiasVoltage(fgkBiasVoltageDefault);
   SetNoiseParam(0.,0.);
   SetDataType("simulated");
}
//_________________________________________________________________________

void AliITSCalibrationSPD::AddDead(UInt_t col, UInt_t row) {
  //
  // Add a dead channel to fDeadChannel array
  //
  fDeadChannels.Set(fNrDead*2+2);
  fDeadChannels.AddAt(col,fNrDead*2);
  fDeadChannels.AddAt(row,fNrDead*2+1);
  fNrDead++;
}
//_________________________________________________________________________
Int_t AliITSCalibrationSPD::GetDeadColAt(UInt_t index) {
  // 
  // Returns column of index-th dead channel
  //
  if (index<fNrDead) {
    return fDeadChannels.At(index*2);
  }
  return -1;
}
//_________________________________________________________________________
Int_t AliITSCalibrationSPD::GetDeadRowAt(UInt_t index) {
  // 
  // Returns row of index-th dead channel
  //
  if (index<fNrDead) {
    return fDeadChannels.At(index*2+1);
  }
  return -1;
}
//_________________________________________________________________________
Bool_t AliITSCalibrationSPD::IsPixelDead(Int_t col, Int_t row) const {
  //
  // Check if pixel (col,row) is dead
  //
  for (UInt_t i=0; i<fNrDead; i++) { 
    if (fDeadChannels.At(i*2)==col && fDeadChannels.At(i*2+1)==row) {
      return true;
    }
  }
  return false;
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//___________________________________________________________________________
// THIS METHOD SHOULD BE DELETED AS SOON AS POSSIBLE!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Bool_t AliITSCalibrationSPD::IsPixelDead(Int_t mod,Int_t ix,Int_t iz) const {
 // Returns kTRUE if pixel is dead
  // Inputs:
  //    Int_t mod      module number
  //    Int_t ix       x pixel number
  //    Int_t iz       z pixel number
  // Outputs:
  //    none.
  // Return:
  //    kFALSE if pixel is alive, or kTRUE if pixel is dead.

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Double_t fDeadPixels = 0.01; // fix to keep AliITSsimulationSPDdubna alive!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Bool_t  dead = kFALSE;
  Int_t   seed;
  static TRandom ran; // don't use gRandom. This must not be a true randome
  // sequence. These sequence must be random one and then fully repetable.

  seed = mod*256*256+iz*256+ix;
  ran.SetSeed(seed);
  if(ran.Rndm(0)<fDeadPixels) dead = kTRUE;
  return dead;
}
//____________________________________________________________________________
void AliITSCalibrationSPD::AddNoisy(UInt_t col, UInt_t row) {
  //
  // add noisy pixel 
  //
  fDeadChannels.Set(fNrNoisy*2+2);
  fNoisyChannels.AddAt(col,fNrNoisy*2);
  fNoisyChannels.AddAt(row,fNrNoisy*2+1);
  fNrNoisy++;
}
//____________________________________________________________________________
Int_t AliITSCalibrationSPD::GetNoisyColAt(UInt_t index) {
  //
  // Get column of index-th noisy pixel
  //
  if (index<fNrNoisy) {
    return fNoisyChannels.At(index*2);
  }
  return -1;
}
//____________________________________________________________________________
Int_t AliITSCalibrationSPD::GetNoisyRowAt(UInt_t index) {
  //
  // Get row of index-th noisy pixel
  //
  if (index<fNrNoisy) {
    return fNoisyChannels.At(index*2+1);
  }
  return -1;
}
//____________________________________________________________________________
Bool_t AliITSCalibrationSPD::IsPixelNoisy(Int_t col, Int_t row) const {
  //
  // Check if pixel (col,row) is noisy
  //
  for (UInt_t i=0; i<fNrNoisy; i++) { 
    if (fNoisyChannels.At(i*2)==col && fNoisyChannels.At(i*2+1)==row) {
      return true;
    }
  }
  return false;
}
