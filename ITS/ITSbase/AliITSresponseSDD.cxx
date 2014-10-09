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
#include <AliITSgeomTGeo.h>

const Float_t AliITSresponseSDD::fgkTimeOffsetDefault = 54.30;
const Float_t AliITSresponseSDD::fgkADC2keVDefault = 3.34;
const Float_t AliITSresponseSDD::fgkChargevsTimeDefault = 0.00355;
const Float_t AliITSresponseSDD::fgkADCvsDrTimeDefault = 0.0101;
const Float_t AliITSresponseSDD::fgkCarlosRXClockPeriod = 25.;
ClassImp(AliITSresponseSDD)

//_________________________________________________________________________
AliITSresponseSDD::AliITSresponseSDD():
TObject(),
  fTimeOffset(fgkTimeOffsetDefault),
  fADC2keV(fgkADC2keVDefault),
  fChargevsTime(fgkChargevsTimeDefault)
{
  // default constructor
  for(Int_t i=0; i<kNSDDmods;i++){
    fTimeZero[i]=fgkTimeOffsetDefault;
    fDeltaVDrift[i] = fDeltaVDrift[i+kNSDDmods] = 0.;
    fADCtokeV[i]=fgkADC2keVDefault;
    fADCvsDriftTime[i]=fgkADCvsDrTimeDefault;
  }  
  SetVDCorr2Side(kTRUE); // default for new objects will be separate corrections for 2 sides (bwd compatible)
  //  SetVDCorrMult(kTRUE); // default for new objects will have multiplicative correction v'=(1+corr)*v (bwd compatible)
}
//_________________________________________________________________________
void AliITSresponseSDD::SetHalfLadderATimeZero(Int_t lay, Int_t lad, Float_t tzero){
  // Sets time Zero for all modules of a ladder on side A (Z>0)
  Int_t minMod,maxMod;
  if(lay==3){
    minMod=1; 
    maxMod=3;
    if(lad>kNLaddersLay3){
      AliError(Form("Ladder number %d out of range",lad));
      return;
    }
  }else if(lay==4){
    minMod=1; 
    maxMod=4;
    if(lad>kNLaddersLay4){
      AliError(Form("Ladder number %d out of range",lad));
      return;
    }
  }else{
    AliError(Form("Layer number %d out of range",lay));
    return;
  }
  for(Int_t iMod=minMod; iMod<=maxMod; iMod++){
    Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(lay,lad,iMod);
    SetModuleTimeZero(modIndex,tzero);
  }
}
//_________________________________________________________________________
void AliITSresponseSDD::SetHalfLadderCTimeZero(Int_t lay, Int_t lad, Float_t tzero){
  // Sets time Zero for all modules of a ladder on side C (Z<0)
  Int_t minMod,maxMod;
  if(lay==3){
    minMod=4; 
    maxMod=6;
    if(lad>kNLaddersLay3){
      AliError(Form("Ladder number %d out of range",lad));
      return;
    }
  }else if(lay==4){
    minMod=5; 
    maxMod=8;
    if(lad>kNLaddersLay4){
      AliError(Form("Ladder number %d out of range",lad));
      return;
    }
  }else{
    AliError(Form("Layer number %d out of range",lay));
    return;
  }
  for(Int_t iMod=minMod; iMod<=maxMod; iMod++){
    Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(lay,lad,iMod);
    SetModuleTimeZero(modIndex,tzero);
  }
}
//_________________________________________________________________________
void AliITSresponseSDD::PrintChargeCalibrationParams() const{
  // Dump charge calibration parameters

  printf("ADC vs. drift time corr=%f\n",GetChargevsTime());
  printf("-------------------------------------\n");
  printf("Layer 3\n");
  for(Int_t ilad=1; ilad<=14; ilad++){
    for(Int_t idet=1; idet<=6;idet++){
      Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(3,ilad,idet);
      Float_t tz=GetADCtokeV(modIndex);
      printf("%7.2f   ",tz);
    }
    printf("\n");
  }
  printf("\n");
  printf("Layer 4\n");
  for(Int_t ilad=1; ilad<=22; ilad++){
    for(Int_t idet=1; idet<=8;idet++){
      Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(4,ilad,idet);
      Float_t tz=GetADCtokeV(modIndex);
      printf("%7.2f   ",tz);
    }
    printf("\n");
  }  
}
//_________________________________________________________________________
void AliITSresponseSDD::PrintTimeZeroes() const{
  // Dump time zero values

  printf("Layer 3\n");
  for(Int_t ilad=1; ilad<=14; ilad++){
    for(Int_t idet=1; idet<=6;idet++){
      Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(3,ilad,idet);
      Float_t tz=GetTimeZero(modIndex);
      printf("%7.2f   ",tz);
    }
    printf("\n");
  }
  printf("\n");
  printf("Layer 4\n");
  for(Int_t ilad=1; ilad<=22; ilad++){
    for(Int_t idet=1; idet<=8;idet++){
      Int_t modIndex=AliITSgeomTGeo::GetModuleIndex(4,ilad,idet);
      Float_t tz=GetTimeZero(modIndex);
      printf("%7.2f   ",tz);
    }
    printf("\n");
  }
  
}
//_________________________________________________________________________
void AliITSresponseSDD::PrintVdriftCorerctions() const{
  // Dump corrections to vdrift

  for(Int_t iMod=240; iMod<500; iMod++){
    printf("Module %d   dVleft=%f   dVright=%f\n",iMod,GetDeltaVDrift(iMod,0),GetDeltaVDrift(iMod,1));
  }
}
