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

//////////////////////////////////////////////////////
//  Base response class forITS                      //
//  It is used to set static data members           //
//  connected to parameters equal for all           //
//  the modules                                     //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

#include <TMath.h>

#include "AliITSresponseSDD.h"

const Int_t AliITSresponseSDD::fgkMaxAdcDefault = 1024;
const Float_t AliITSresponseSDD::fgkDynamicRangeDefault = 132.;
const Float_t AliITSresponseSDD::fgkfChargeLossDefault = 0;
const Float_t AliITSresponseSDD::fgkDiffCoeffDefault = 3.23;
const Float_t AliITSresponseSDD::fgkDiffCoeff1Default = 30.;
const TString AliITSresponseSDD::fgkParam1Default = "same";
const TString AliITSresponseSDD::fgkParam2Default = "same";
const TString AliITSresponseSDD::fgkOptionDefault = "ZS";
const Float_t AliITSresponseSDD::fgkDriftSpeedDefault = 7.3;
const Float_t AliITSresponseSDD::fgkTimeOffsetDefault = 53.57;
const Float_t AliITSresponseSDD::fgkADC2keVDefault = 5.243;
const Float_t AliITSresponseSDD::fgkNsigmasDefault = 3.;
const Int_t AliITSresponseSDD::fgkNcompsDefault = 121;

ClassImp(AliITSresponseSDD)

//_________________________________________________________________________
AliITSresponseSDD::AliITSresponseSDD():
AliITSresponse(),
fJitterError(0.),
fDynamicRange(0.),
fChargeLoss(0.),
fTimeOffset(fgkTimeOffsetDefault),
fADC2keV(fgkADC2keVDefault),
fElectronics(0),
fMaxAdc(fgkMaxAdcDefault),
fNsigmas(fgkNsigmasDefault),
fGaus(),
fNcomps(fgkNcompsDefault),
fOption(),
fParam1(),
fParam2() {
  // default constructor
  fGaus = 0;
  SetDiffCoeff(fgkDiffCoeffDefault,fgkDiffCoeff1Default);
  //  SetNLookUp(fgkNcompsDefault);

  SetJitterError();
  SetElectronics();
  SetDynamicRange(fgkDynamicRangeDefault);
  SetChargeLoss(fgkfChargeLossDefault);
  SetParamOptions(fgkParam1Default.Data(),fgkParam2Default.Data());
  SetZeroSupp(fgkOptionDefault);
  SetOutputOption();
}


//______________________________________________________________________
AliITSresponseSDD::~AliITSresponseSDD() { 

  if(fGaus) delete fGaus;
}


//________________________________________________________________________
void AliITSresponseSDD::SetNLookUp(Int_t p1){
  // Set number of sigmas over which cluster disintegration is performed
  fNcomps=p1;
  if (fGaus) delete fGaus;
  fGaus = new TArrayF(fNcomps+1);
  for(Int_t i=0; i<=fNcomps; i++) {
    Float_t x = -fNsigmas + (2.*i*fNsigmas)/(fNcomps-1);
    (*fGaus)[i] = exp(-((x*x)/2));
  }
}
