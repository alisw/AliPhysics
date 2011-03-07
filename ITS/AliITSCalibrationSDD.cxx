/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

#include <Riostream.h>
#include <TRandom.h>
#include "AliITSCalibrationSDD.h"
#include "AliLog.h"

//////////////////////////////////////////////////////
//  Calibration class for set:ITS                   //
//  Specific subdetector implementation             //
//  for silicon drift detectors                     //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

const Float_t AliITSCalibrationSDD::fgkTemperatureDefault = 296.;
const Float_t AliITSCalibrationSDD::fgkNoiseDefault = 2.38;
const Float_t AliITSCalibrationSDD::fgkGainDefault = 1.;
const Float_t AliITSCalibrationSDD::fgkBaselineDefault = 20.;
//______________________________________________________________________
ClassImp(AliITSCalibrationSDD)

AliITSCalibrationSDD::AliITSCalibrationSDD():
AliITSCalibration(),
fZeroSupp(kTRUE),
fAMAt20MHz(kFALSE),
fDeadChips(0),
fDeadChannels(0),
fIsBad(kFALSE),
fBadChannels(),
fMapAW0(0),
fMapAW1(0),
fMapTW0(0),
fMapTW1(0),
fDrSpeed0(0),
fDrSpeed1(0)
{
  // default constructor

  SetDeadChannels();
  for(Int_t ian=0;ian<fgkWings*fgkChannels*fgkChips;ian++){
    fBaseline[ian]=fgkBaselineDefault;
    fNoise[ian]=fgkNoiseDefault;
    fGain[ian]=1.;
    SetNoiseAfterElectronics(ian);
  }
  for(Int_t iw=0;iw<fgkWings;iw++){
    SetZSLowThreshold(iw);
    SetZSHighThreshold(iw);
    for(Int_t icp=0;icp<fgkChips;icp++){
      Int_t chipindex=iw*fgkChips+icp;
      fIsChipBad[chipindex]=kFALSE;
    }
  }
  SetTemperature(fgkTemperatureDefault);
  SetDataType();
 }
//______________________________________________________________________
AliITSCalibrationSDD::AliITSCalibrationSDD(const char *dataType):
AliITSCalibration(),
fZeroSupp(kTRUE),
fAMAt20MHz(kFALSE),
fDeadChips(0),
fDeadChannels(0),
fIsBad(kFALSE),
fBadChannels(),
fMapAW0(0),
fMapAW1(0),
fMapTW0(0),
fMapTW1(0),
fDrSpeed0(0),
fDrSpeed1(0)
{
  // constructor

  SetDeadChannels();
  for(Int_t ian=0;ian<fgkWings*fgkChannels*fgkChips;ian++){
    fBaseline[ian]=fgkBaselineDefault;
    fNoise[ian]=fgkNoiseDefault;
    fGain[ian]=1.;
    SetNoiseAfterElectronics(ian);
  }  
  for(Int_t iw=0;iw<fgkWings;iw++){
    SetZSLowThreshold(iw);
    SetZSHighThreshold(iw);
    for(Int_t icp=0;icp<fgkChips;icp++){
      Int_t chipindex=iw*fgkChips+icp;
      fIsChipBad[chipindex]=kFALSE;
    }
  }

  SetTemperature(fgkTemperatureDefault);
  SetDataType(dataType);
 }
//_____________________________________________________________________
AliITSCalibrationSDD::~AliITSCalibrationSDD(){

  //destructor
  if(fMapAW0) delete fMapAW0;
  if(fMapAW1) delete fMapAW1;
  if(fMapTW0) delete fMapTW0;
  if(fMapTW1) delete fMapTW1;
  if(fDrSpeed0) delete fDrSpeed0;
  if(fDrSpeed1) delete fDrSpeed1;
}

//______________________________________________________________________
void AliITSCalibrationSDD::GiveCompressParam(Int_t  cp[4]) const {
  // give compression param
  cp[0]=fZSTH[0];
  cp[1]=fZSTL[0];
  cp[2]=fZSTH[1];
  cp[3]=fZSTL[1];
}
//_____________________________________________________________________
void AliITSCalibrationSDD::SetBadChannel(Int_t i,Int_t anode){
  //Set bad anode (set gain=0 for these channels);

  if(anode<0 || anode >fgkChannels*fgkChips*fgkWings-1){
    AliError("Wrong anode number");
    return;
  }
  fBadChannels[i]=anode;
  fGain[anode]=0;
}
//______________________________________________________________________
void AliITSCalibrationSDD::GetCorrections(Float_t z, Float_t x, Float_t &devz, Float_t &devx, AliITSsegmentationSDD* seg){
  //correction of coordinates using the maps stored in the DB
  Int_t nSide=seg->GetSideFromLocalX(x);
  devz=0;
//     if(nSide==0) devz=fMapAW0->GetCorrection(z,x,seg);
//     else devz=fMapAW1->GetCorrection(z,x,seg);
  devx=0;
  if(nSide==0) devx=fMapTW0->GetCorrection(z,x,seg);
  else devx=fMapTW1->GetCorrection(z,x,seg);
  return;
}
//______________________________________________________________________
void AliITSCalibrationSDD::GetShiftsForSimulation(Float_t z, Float_t x, Float_t &devz, Float_t &devx, AliITSsegmentationSDD* seg){
  //correction of coordinates using the maps stored in the DB
  Int_t nSide=seg->GetSideFromLocalX(x);
  devz=0;
//     if(nSide==0) devz=fMapAW0->GetCorrection(z,x,seg);
//     else devz=fMapAW1->GetCorrection(z,x,seg);
  devx=0;
  if(nSide==0) devx=fMapTW0->GetShiftForSimulation(z,x,seg);
  else devx=fMapTW1->GetShiftForSimulation(z,x,seg);
  return;
}
//______________________________________________________________________
void AliITSCalibrationSDD::PrintGains() const{
  // Print Gains

  if( GetDeadChips() == 0 && 
      GetDeadChannels() == 0 )
    return;  

  // Print Electronics Gains
  cout << "**************************************************" << endl; 
  cout << "             Print Electronics Gains              " << endl;
  cout << "**************************************************" << endl;

  // Print SDD electronic gains
  for(Int_t ian=0; ian<fgkWings*fgkChips*fgkChannels;ian++){
    printf("Gain for channel %d = %f\n",ian,fGain[ian]);
  }
}

//______________________________________________________________________
void AliITSCalibrationSDD::Print(){
  // Print SDD response Parameters

  cout << "**************************************************" << endl;
  cout << "   Silicon Drift Detector Response Parameters    " << endl;
  cout << "**************************************************" << endl;
  cout << "Hardware compression parameters: " << endl; 
  cout << "Noise before electronics (arbitrary units): " << fNoise[0] << endl;
  cout << "Baseline (ADC units): " << fBaseline[0] << endl;
  cout << "Noise after electronics (ADC units): " << fNoiseAfterEl[0] << endl;
  cout << "Temperature: " << Temperature() << " K " << endl;
  PrintGains();

}
