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

//////////////////////////////////////////////////////
//  Calibration class for set:ITS                   //
//  Specific subdetector implementation             //
//  for silicon drift detectors                     //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

const Double_t AliITSCalibrationSDD::fgkTemperatureDefault = 296.;
const Double_t AliITSCalibrationSDD::fgkNoiseDefault = 10.;
const Double_t AliITSCalibrationSDD::fgkBaselineDefault = 20.;
const Double_t AliITSCalibrationSDD::fgkMinValDefault  = 4;
//______________________________________________________________________
ClassImp(AliITSCalibrationSDD)

  AliITSCalibrationSDD::AliITSCalibrationSDD(){
  // default constructor

  SetDeadChannels();
  fBaseline=0;
  fNoise=0;

  SetNoiseParam(fgkNoiseDefault,fgkBaselineDefault);
  SetNoiseAfterElectronics();
  SetThresholds(fgkMinValDefault,0.);
  SetTemperature(fgkTemperatureDefault);
  SetDataType();
  fCPar[1]=(Int_t) fBaseline;
  fCPar[2]=(Int_t)(2.*fNoiseAfterEl + 0.5);
  fCPar[3]=(Int_t)(2.*fNoiseAfterEl + 0.5);
  fCPar[4]=0;
  fCPar[5]=0;
  fCPar[6]=0;
  fCPar[7]=0;
 }
//______________________________________________________________________
AliITSCalibrationSDD::AliITSCalibrationSDD(const char *dataType){
  // constructor

  SetDeadChannels();
  fBaseline=0;
  fNoise=0;

  SetNoiseParam(fgkNoiseDefault,fgkBaselineDefault);
  SetNoiseAfterElectronics();
  SetThresholds(fgkMinValDefault,0.);
  SetTemperature(fgkTemperatureDefault);
  SetDataType(dataType);
  fCPar[0]=(Int_t) fBaseline;
  fCPar[1]=(Int_t) fBaseline;
  fCPar[2]=(Int_t)(2.*fNoiseAfterEl + 0.5);
  fCPar[3]=(Int_t)(2.*fNoiseAfterEl + 0.5);
  fCPar[4]=0;
  fCPar[5]=0;
  fCPar[6]=0;
  fCPar[7]=0;
 }
//______________________________________________________________________
AliITSCalibrationSDD::AliITSCalibrationSDD(const AliITSCalibrationSDD &ob) : AliITSCalibration(ob) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSCalibrationSDD","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliITSCalibrationSDD& AliITSCalibrationSDD::operator=(const AliITSCalibrationSDD& /* ob */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
  return *this;
}

//______________________________________________________________________
void AliITSCalibrationSDD::SetCompressParam(Int_t  cp[8]){
  // set compression param

  Int_t i;
  for (i=0; i<8; i++) {
    fCPar[i]=cp[i];
    //printf("\n CompressPar %d %d \n",i,fCPar[i]);    
  } // end for i
}
//______________________________________________________________________
void AliITSCalibrationSDD::GiveCompressParam(Int_t  cp[8]) const {
  // give compression param

  Int_t i;
  for (i=0; i<8; i++) {
    cp[i]=fCPar[i];
  } // end for i
}

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
  /*
  //  cout << "modules" << endl;
  Int_t * mod = new Int_t [nmod];
  Int_t i; //loop variable
  for( i=0; i<nmod; i++ ) 
    {
      mod[i] = (Int_t) (1.+fgkModules*gran->Uniform());
      cout << i+1 << ": Dead module nr: " << mod[i] << endl;
      for(Int_t n=0; n<fResponseSDD->Chips(); n++)
	for(Int_t p=0; p<fResponseSDD->Channels(); p++)
	  fGain[mod[i]-1][n][p] = 0.;
    }
  */
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
  for(Int_t i=0; i<8; i++) cout << "fCPar[" << i << "] = " << fCPar[i] <<endl;
  cout << "Noise before electronics (arbitrary units): " << fNoise << endl;
  cout << "Baseline (ADC units): " << fBaseline << endl;
  cout << "Noise after electronics (ADC units): " << fNoiseAfterEl << endl;
  cout << "Temperature: " << Temperature() << " K " << endl;
  cout << "Min. Value: " << fMinVal << endl;
  PrintGains();

}



