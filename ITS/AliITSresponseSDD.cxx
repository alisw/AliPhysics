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
const Float_t AliITSresponseSDD::fgkADC2keVDefault = 5.243;
const Float_t AliITSresponseSDD::fgkCarlosRXClockPeriod = 25.;
ClassImp(AliITSresponseSDD)

//_________________________________________________________________________
AliITSresponseSDD::AliITSresponseSDD():
TObject(),
fTimeOffset(fgkTimeOffsetDefault),
fADC2keV(fgkADC2keVDefault){
  // default constructor
  for(Int_t i=0; i<kNSDDmods;i++){
    fTimeZero[i]=fgkTimeOffsetDefault;
  }  
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
