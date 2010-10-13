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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for monitoring data                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalMonitoring.h"

ClassImp(AliTRDCalMonitoring)

//_____________________________________________________________________________
AliTRDCalMonitoring::AliTRDCalMonitoring()
  :TNamed()
  ,fDriftVelocity(0)
  ,fGasComposition(0)
  ,fEnvironmentTemperature(0)
{
  //
  // AliTRDCalMonitoring default constructor
  //
  
  for (Int_t i=0; i<540; ++i) {
    fAnodeCurrentsMin[i] = 0;
    fAnodeCurrentsMax[i] = 0;
    fDriftCurrentsMin[i] = 0;
    fDriftCurrentsMax[i] = 0;
    fAnodeVoltagesMin[i] = 0;
    fAnodeVoltagesMax[i] = 0;
    fDriftVoltagesMin[i] = 0;
    fDriftVoltagesMax[i] = 0;
  }
  for (Int_t i=0; i<360; ++i) {
    fLVVoltage[i] = 0;
    fLVCurrent[i] = 0;
  }
  for (Int_t i=0; i<6700; ++i) {
    fADCTresholds[i] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalMonitoring::AliTRDCalMonitoring(const Text_t *name, const Text_t *title)
  :TNamed(name,title)
  ,fDriftVelocity(0)
  ,fGasComposition(0)
  ,fEnvironmentTemperature(0)
{
  //
  // AliTRDCalMonitoring constructor
  //

  for (Int_t i=0; i<540; ++i) {
    fAnodeCurrentsMin[i] = 0;
    fAnodeCurrentsMax[i] = 0;
    fDriftCurrentsMin[i] = 0;
    fDriftCurrentsMax[i] = 0;
    fAnodeVoltagesMin[i] = 0;
    fAnodeVoltagesMax[i] = 0;
    fDriftVoltagesMin[i] = 0;
    fDriftVoltagesMax[i] = 0;
  }
  for (Int_t i=0; i<360; ++i) {
    fLVVoltage[i] = 0;
    fLVCurrent[i] = 0;
  }
  for (Int_t i=0; i<6700; ++i) {
    fADCTresholds[i] = 0;
  }

}

