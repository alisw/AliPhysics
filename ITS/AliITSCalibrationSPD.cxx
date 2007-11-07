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
///////////////////////////////////////////////////////////////////////////
//  Calibration class for set:ITS                   
//  Specific subdetector implementation for         
//  Silicon pixels                                  
//
//  Modified by D. Elia, G.E. Bruno, H. Tydesjo
//  March-April 2006
//  Last mod:  H. Tydesjo  Oct 2007
//  September   2007: CouplingRowDefault = 0.055 (was 0.047)
//
///////////////////////////////////////////////////////////////////////////
const Double_t AliITSCalibrationSPD::fgkThreshDefault = 3000.;
const Double_t AliITSCalibrationSPD::fgkSigmaDefault = 250.;
const Double_t AliITSCalibrationSPD::fgkCouplColDefault = 0.;
const Double_t AliITSCalibrationSPD::fgkCouplRowDefault = 0.055;
const Double_t AliITSCalibrationSPD::fgkBiasVoltageDefault = 18.182;

ClassImp(AliITSCalibrationSPD)	

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
fNrBad(0),
fBadChannels(0){
  // constructor

   SetThresholds(fgkThreshDefault,fgkSigmaDefault);
   SetCouplingParam(fgkCouplColDefault,fgkCouplRowDefault);
   SetBiasVoltage(fgkBiasVoltageDefault);
   SetNoiseParam(0.,0.);
   SetDataType("simulated");
}
//____________________________________________________________________________
void AliITSCalibrationSPD::AddBad(UInt_t col, UInt_t row) {
  //
  // add bad pixel 
  //
  fBadChannels.Set(fNrBad*2+2);
  fBadChannels.AddAt(col,fNrBad*2);
  fBadChannels.AddAt(row,fNrBad*2+1);
  fNrBad++;
}
//____________________________________________________________________________
Int_t AliITSCalibrationSPD::GetBadColAt(UInt_t index) {
  //
  // Get column of index-th bad pixel
  //
  if (index<fNrBad) {
    return fBadChannels.At(index*2);
  }
  return -1;
}
//____________________________________________________________________________
Int_t AliITSCalibrationSPD::GetBadRowAt(UInt_t index) {
  //
  // Get row of index-th bad pixel
  //
  if (index<fNrBad) {
    return fBadChannels.At(index*2+1);
  }
  return -1;
}
//____________________________________________________________________________
Bool_t AliITSCalibrationSPD::IsPixelBad(Int_t col, Int_t row) const {
  //
  // Check if pixel (col,row) is bad
  //
  for (UInt_t i=0; i<fNrBad; i++) { 
    if (fBadChannels.At(i*2)==col && fBadChannels.At(i*2+1)==row) {
      return true;
    }
  }
  return false;
}
