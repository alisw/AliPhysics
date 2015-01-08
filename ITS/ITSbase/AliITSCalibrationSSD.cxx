/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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
/* $Id$  */

#include "AliITSCalibrationSSD.h"
//////////////////////////////////////////////////////
//  Calibration class for set:ITS                   //
//  Specific subdetector implementation             //
//  for silicon strips detectors                    //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

const Int_t AliITSCalibrationSSD::fgkNParDefault = 6;

ClassImp(AliITSCalibrationSSD)

//______________________________________________________________________
AliITSCalibrationSSD::AliITSCalibrationSSD():
  fModule(0),
fNPar(0),
fDetPar(0),
fNoise(0),
fPedestal(),
fGain(0),
fBadChannels(0),
  fIsBad(kFALSE),
fSSDADCpereV(0.),
  fKeVperADC(0)
{
    // Default Constructor

    for(Int_t i=0;i<fgkChipsPerModule;i++){
      fIsChipBad[i]=kFALSE;
    }
  SetSSDADCpereV();
    SetKeVperADC();
}
//______________________________________________________________________
AliITSCalibrationSSD::AliITSCalibrationSSD(const char *dataType):
  fModule(0),
fNPar(0),
fDetPar(0),
fNoise(0),
fPedestal(0),
fGain(0),
fBadChannels(0),
fIsBad(kFALSE) ,
fSSDADCpereV(0.),
fKeVperADC(0){
    // constructor

    SetDataType(dataType);
    SetNDetParam(fgkNParDefault);   // Sets fNPar=6 by default.
    fDetPar = new Double_t[fNPar];
    if (fNPar==6) {
	fDetPar[0]=10.;
	fDetPar[1]=5.;
	fDetPar[2]=0.02;
	fDetPar[3]=0.02;
	fDetPar[4]=0.02;
	fDetPar[5]=0.03;
    } // end if
    for(Int_t i=0;i<fgkChipsPerModule;i++){
      fIsChipBad[i]=kFALSE;
    }
    SetSSDADCpereV();
    SetKeVperADC();
}
//______________________________________________________________________
AliITSCalibrationSSD::~AliITSCalibrationSSD(){
    // destructor
 
    delete [] fDetPar;
    if(fNoise)delete fNoise;
    if(fPedestal)delete fPedestal;
    if(fGain)delete fGain;
    if(fBadChannels)delete fBadChannels;
}
//______________________________________________________________________
void AliITSCalibrationSSD::SetDetParam(Double_t  *par){
    // set det param
    Int_t i;

    for (i=0; i<fNPar; i++) {
	fDetPar[i]=par[i];
	//printf("\n CompressPar %d %d \n",i,fCPar[i]);    
    } // end for i
}
//______________________________________________________________________
void AliITSCalibrationSSD::GetDetParam(Double_t  *par) const {
    // get det param
    Int_t i;

    for (i=0; i<fNPar; i++) {
	par[i]=fDetPar[i];
    } // end for i
}

//______________________________________________________________________
void AliITSCalibrationSSD::FillBadChipMap() {

  // reset
  fIsBad=kFALSE;
  for(Int_t i=0;i<fgkChipsPerModule;i++){
    fIsChipBad[i]=kFALSE;
  }


  Int_t mc=0;
  Int_t cc[fgkChipsPerModule];

  // P-side
  for(Int_t i=0; i<fgkChipsPerModule/2; i++){
    cc[i]=0;
    for(Int_t j=0; j<ChannelsPerChip(); j++) {
      if(IsPChannelBad(i*ChannelsPerChip()+j)) cc[i]++;
    }
    if(cc[i]==ChannelsPerChip()) { SetChipBad(i); mc++; }
  }
  
  // N-side
  for(Int_t i=fgkChipsPerModule/2; i<fgkChipsPerModule; i++){
    cc[i]=0;
    for(Int_t j=0; j<ChannelsPerChip(); j++) {
      if(IsNChannelBad(1535-i*ChannelsPerChip()-j)) cc[i]++;      
    }
    if(cc[i]==ChannelsPerChip()) { SetChipBad(i); mc++; }
  }
  
  if(mc==fgkChipsPerModule) fIsBad=kTRUE;
}
