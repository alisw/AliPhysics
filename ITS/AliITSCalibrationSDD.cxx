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
const Float_t AliITSCalibrationSDD::fgkNoiseDefault = 10.;
const Float_t AliITSCalibrationSDD::fgkGainDefault = 1.;
const Float_t AliITSCalibrationSDD::fgkBaselineDefault = 20.;
const Float_t AliITSCalibrationSDD::fgkMinValDefault  = 4;
//______________________________________________________________________
ClassImp(AliITSCalibrationSDD)

AliITSCalibrationSDD::AliITSCalibrationSDD():
AliITSCalibration(),
fDeadChips(0),
fDeadChannels(0),
fMinVal(fgkMinValDefault),
fIsDead(kFALSE),
fBadChannels()
{
  // default constructor

  SetDeadChannels();
  for(Int_t ian=0;ian<fgkWings*fgkChannels*fgkChips;ian++){
    fBaseline[ian]=fgkBaselineDefault;
    fNoise[ian]=fgkNoiseDefault;
    SetNoiseAfterElectronics(ian);
  }
  for(Int_t iw=0;iw<fgkWings;iw++){
    for(Int_t icp=0;icp<fgkChips;icp++){
      for(Int_t ich=0;ich<fgkChannels;ich++)
	fGain[iw][icp][ich]=1.;
    }
  }
  SetThresholds(fgkMinValDefault,0.);
  SetTemperature(fgkTemperatureDefault);
  SetDataType();
  for(Int_t i=0;i<fgkChips*fgkChannels;i++){
    for(Int_t j=0;j<fgkMapTimeNBin;j++){
      fMapA[i][j]=0;
      fMapT[i][j]=0;
    }
  }
 }
//______________________________________________________________________
AliITSCalibrationSDD::AliITSCalibrationSDD(const char *dataType):
AliITSCalibration(),
fDeadChips(0),
fDeadChannels(0),
fMinVal(fgkMinValDefault),
fIsDead(kFALSE),
fBadChannels(){
  // constructor

  SetDeadChannels();
  for(Int_t ian=0;ian<fgkWings*fgkChannels*fgkChips;ian++){
    fBaseline[ian]=fgkBaselineDefault;
      fNoise[ian]=fgkNoiseDefault;
      SetNoiseAfterElectronics(ian);
  }  
  for(Int_t iw=0;iw<fgkWings;iw++){
    for(Int_t icp=0;icp<fgkChips;icp++){
      for(Int_t ich=0;ich<fgkChannels;ich++)
	fGain[iw][icp][ich]=1.;
    }
  }

  SetThresholds(fgkMinValDefault,0.);
  SetTemperature(fgkTemperatureDefault);
  SetDataType(dataType);
  for(Int_t i=0;i<fgkChips*fgkChannels;i++){
    for(Int_t j=0;j<fgkMapTimeNBin;j++){
      fMapA[i][j]=0;
      fMapT[i][j]=0;
    }
  }

 }


//______________________________________________________________________
void AliITSCalibrationSDD::GiveCompressParam(Int_t  cp[8],Int_t ian) const {
  // give compression param

  cp[0]=(Int_t) fBaseline[ian];
  cp[1]=(Int_t) fBaseline[ian];
  cp[2]=(Int_t)(2.*fNoiseAfterEl[ian] + 0.5);
  cp[3]=(Int_t)(2.*fNoiseAfterEl[ian] + 0.5);
  cp[4]=0;
  cp[5]=0;
  cp[6]=0;
  cp[7]=0;
}
//_____________________________________________________________________
void AliITSCalibrationSDD::SetBadChannel(Int_t i,Int_t anode){
  //Set bad anode (set gain=0 for these channels);

  if(anode<0 || anode >fgkChannels*fgkChips*fgkWings-1)AliError("Wrong anode number");
  Int_t wing=0;
  Int_t chip,channel;
  chip=anode/fgkChannels;
  channel=anode-(chip*fgkChannels);
  if(anode>=fgkChips*fgkChannels) wing=1;
  if(wing==1)chip-=fgkChips;
  fBadChannels[i]=anode;
  fGain[wing][chip][channel]=0;
}
//_____________________________________________________________________
Bool_t AliITSCalibrationSDD::IsBadChannel(Int_t anode){
  //returns kTRUE if the anode i (0-512) has fGain=0
  if(anode<0 || anode >fgkChannels*fgkChips*fgkWings-1)AliError("Wrong anode number");
  Int_t wing=0;
  Int_t chip,channel;
  chip=anode/fgkChannels;
  channel=anode-(chip*fgkChannels);
  if(anode>=fgkChips*fgkChannels) wing=1;
  if(wing==1)chip-=fgkChips;
  if(fGain[wing][chip][channel]==0) return kTRUE;
  else return kFALSE;
} 
/*
//______________________________________________________________________
void AliITSCalibrationSDD::SetDeadChannels(Int_t nchip, Int_t nchan){
  // Set fGain to zero to simulate a random distribution of 
  // dead modules, dead chips and single dead channels

  for( Int_t m=0; m<fgkWings; m++ ) 
    for( Int_t n=0; n<fgkChips; n++ ) 
      for( Int_t p=0; p<fgkChannels; p++ ) 
	fGain[m][n][p] = 1.;
                 
  //fDeadModules  = nmod;  
  fDeadChips    = nchip;  
  fDeadChannels = nchan; 
  fBadChannels.Set(fDeadChannels);  
  // nothing to do
  //if( nmod == 0 && nchip == 0 && nchan == 0 ) return;

  if( nchip == 0 && nchan == 0 ) return;
  // if( nmod < 0 || nmod > fgkModules ) 
  //  { 
  //    cout << "Wrong number of dead modules: " << nmod << endl; 
  //    return; 
  //  }
  
  Int_t nmax = fgkWings*fgkChips; 
  if( nchip < 0 || nchip > nmax ) 
    { 
      cout << "Wrong number of dead chips: " << nchip << endl; 
      return; 
    }
  nmax = (fgkWings*fgkChips - nchip)*fgkChannels; 
  if( nchan < 0 || nchan > nmax ) 
    { 
      cout << "Wrong number of dead channels: " << nchan << endl; 
      return; 
    }
  
  TRandom *gran = new TRandom();
  //  cout << "chips" << endl;
  Int_t * chip     = new Int_t[nchip];
  Int_t i = 0;
  while( i < nchip ) 
    {
      Int_t wing = (Int_t) (fgkWings*gran->Uniform() + 1.);
      if( wing <=0 || wing > fgkWings ) Error("SetDeadChannels","Wrong wing");
        
      Int_t chi = (Int_t) (fgkChips*gran->Uniform() + 1.);
      if( chi <=0 || chi > fgkChips ) Error("SetDeadChannels","Wrong chip:%d\n",chi);
      i++;
      chip[i-1] = chi; 
      for( Int_t m=0; m<fgkChannels; m++ ) 
	fGain[wing-1][chi-1][m] = 0.;
    }

  Int_t * channel      = new Int_t[nchan];
  Int_t * channelChip = new Int_t[nchan];
  i = 0;
  while( i < nchan ) 
    {
      Int_t k; //loop variable
      Int_t wing = (Int_t) (fgkWings*gran->Uniform() + 1.);
      if( wing <=0 || wing > fgkWings ) Error("SetDeadChannels","Wrong wing:%d\n",wing);
      Int_t chipp = (Int_t) (fgkChips*gran->Uniform() + 1.);
      if( chipp <=0 || chipp > fgkChips ) Error("SetDeadChannels","Wrong chip:%d",chipp);
      Int_t flagChip = 0;
      for( k=0; k<nchip; k++) 
	if( chipp == chip[k] ) { 
	  flagChip = 1; break; }
      if( flagChip == 1 ) continue;
      i++;
      channel[i-1] = (Int_t) (fgkChannels*gran->Uniform() + 1.); 
      if( channel[i-1] <=0 || channel[i-1] > fgkChannels ) 
	Error("SetDeadChannels","Wrong channel:%d\n",channel[i-1]);
      channelChip[i-1] = chipp;
      fGain[wing-1][chipp-1][channel[i-1]-1] = 0.;
    }
    
  delete [] chip;
  delete [] channel;
  delete [] channelChip;
}
*/
//______________________________________________________________________
void AliITSCalibrationSDD::PrintGains() const{
  //

  if( GetDeadChips() == 0 && 
      GetDeadChannels() == 0 )
    return;  

  // Print Electronics Gains
  cout << "**************************************************" << endl; 
  cout << "             Print Electronics Gains              " << endl;
  cout << "**************************************************" << endl;

  // Print SDD electronic gains
  for(Int_t t=0; t<fgkWings;t++)
    for(Int_t u=0; u<fgkChips;u++)
      for(Int_t v=0; v<fgkChannels;v++)
	{
	  if( fGain[t][u][v] != 1.0 )
	    cout << "Gain for wing: " << t+1 << ", Chip " << u+1 << 
	      ", Channel " << v+1 << " = " << fGain[t][u][v] << endl;
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
  cout << "Min. Value: " << fMinVal << endl;
  PrintGains();

}
