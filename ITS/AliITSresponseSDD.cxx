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

#include <Riostream.h>
#include <TRandom.h>

#include "AliITSresponseSDD.h"

//////////////////////////////////////////////////
//  Response class for set:ITS                      //
//  Specific subdetector implementation             //
//  for silicon drift detectors                     //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

const Int_t AliITSresponseSDD::fgkModules;   
const Int_t AliITSresponseSDD::fgkChips;  
const Int_t AliITSresponseSDD::fgkChannels; 
const Int_t AliITSresponseSDD::fgkMaxAdcDefault = 1024;
const Double_t AliITSresponseSDD::fgkDynamicRangeDefault = 132.;
const Double_t AliITSresponseSDD::fgkfChargeLossDefault = 0;
const Double_t AliITSresponseSDD::fgkDiffCoeffDefault = 3.23;
const Double_t AliITSresponseSDD::fgkDiffCoeff1Default = 30.;
const Double_t AliITSresponseSDD::fgkTemperatureDefault = 296.;
const TString AliITSresponseSDD::fgkParam1Default = "same";
const TString AliITSresponseSDD::fgkParam2Default = "same";
const Double_t AliITSresponseSDD::fgkNoiseDefault = 10.;
const Double_t AliITSresponseSDD::fgkBaselineDefault = 20.;
const TString AliITSresponseSDD::fgkOptionDefault = "1D";
const Double_t AliITSresponseSDD::fgkMinValDefault  = 4;
const Double_t AliITSresponseSDD::fgkDriftSpeedDefault = 7.3;
const Double_t AliITSresponseSDD::fgkNsigmasDefault = 3.;
const Int_t AliITSresponseSDD::fgkNcompsDefault = 121;
//______________________________________________________________________
ClassImp(AliITSresponseSDD)

  AliITSresponseSDD::AliITSresponseSDD(){
  // default constructor
  fGaus = 0;
  SetDeadChannels();
  SetMaxAdc(fgkMaxAdcDefault);
  SetDiffCoeff(fgkDiffCoeffDefault,fgkDiffCoeff1Default);
  SetDriftSpeed(fgkDriftSpeedDefault);
  SetNSigmaIntegration(fgkNsigmasDefault);
  SetNLookUp(fgkNcompsDefault);
  // SetClock();
  SetNoiseParam(fgkNoiseDefault,fgkBaselineDefault);
  SetNoiseAfterElectronics();
  SetJitterError();
  SetElectronics();
  SetDynamicRange(fgkDynamicRangeDefault);
  SetChargeLoss(fgkfChargeLossDefault);
  SetThresholds(fgkMinValDefault,0.);
  SetParamOptions(fgkParam1Default.Data(),fgkParam2Default.Data());
  SetTemperature(fgkTemperatureDefault);
  SetZeroSupp(fgkOptionDefault);
  SetDataType();
  SetFilenames();
  SetOutputOption();
  SetDo10to8();
  // set the default zero suppression parameters
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
AliITSresponseSDD::AliITSresponseSDD(const char *dataType){
  // constructor
  fGaus = 0;
  SetDeadChannels();
  SetMaxAdc(fgkMaxAdcDefault);
  SetDiffCoeff(fgkDiffCoeffDefault,fgkDiffCoeff1Default);
  SetDriftSpeed(fgkDriftSpeedDefault);
  SetNSigmaIntegration(fgkNsigmasDefault);
  SetNLookUp(fgkNcompsDefault);
  // SetClock();
  SetNoiseParam(fgkNoiseDefault,fgkBaselineDefault);
  SetNoiseAfterElectronics();
  SetJitterError();
  SetElectronics();
  SetDynamicRange(fgkDynamicRangeDefault);
  SetChargeLoss(fgkfChargeLossDefault);
  SetThresholds(fgkMinValDefault,0.);
  SetParamOptions(fgkParam1Default.Data(),fgkParam2Default.Data());
  SetTemperature(fgkTemperatureDefault);
  SetZeroSupp(fgkOptionDefault);
  SetDataType(dataType);
  SetFilenames();
  SetOutputOption();
  SetDo10to8();
  // set the default zero suppression parameters
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
AliITSresponseSDD::AliITSresponseSDD(const AliITSresponseSDD &ob) : AliITSresponse(ob) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSresponseSDD","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliITSresponseSDD& AliITSresponseSDD::operator=(const AliITSresponseSDD& /* ob */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
  return *this;
}

//______________________________________________________________________
AliITSresponseSDD::~AliITSresponseSDD() { 

  if(fGaus) delete fGaus;
}

//______________________________________________________________________
Int_t AliITSresponseSDD::Convert8to10(Int_t signal) const {
  // Undo the lossive 10 to 8 bit compression.
  // code from Davide C. and Albert W.
  if(Do10to8()){  // kTRUE if the compression is active
    if (signal < 0 || signal > 255) {
      Warning("Convert8to10","out of range signal=%d",signal);
      return 0;
    } // end if signal <0 || signal >255

    if (signal < 128) return signal;
    if (signal < 192) {
      if (TMath::Odd(signal)) return (128+((signal-128)<<1));
      else  return (128+((signal-128)<<1)+1);
    } // end if signal < 192
    if (signal < 224) {
      if (TMath::Odd(signal)) return (256+((signal-192)<<3)+3);
      else  return (256+((signal-192)<<3)+4);
    } // end if signal < 224
    if (TMath::Odd(signal)) return (512+((signal-224)<<4)+7);
    return (512+((signal-224)<<4)+8);
  }
  else {  
    return signal;
  }
}

//______________________________________________________________________
void AliITSresponseSDD::SetCompressParam(Int_t  cp[8]){
  // set compression param

  Int_t i;
  for (i=0; i<8; i++) {
    fCPar[i]=cp[i];
    //printf("\n CompressPar %d %d \n",i,fCPar[i]);    
  } // end for i
}
//______________________________________________________________________
void AliITSresponseSDD::GiveCompressParam(Int_t  cp[8]) const {
  // give compression param

  Int_t i;
  for (i=0; i<8; i++) {
    cp[i]=fCPar[i];
  } // end for i
}
//______________________________________________________________________
void AliITSresponseSDD::SetNLookUp(Int_t p1){
  // Set number of sigmas over which cluster disintegration is performed
  fNcomps=p1;
  fGaus = new TArrayF(fNcomps+1);
  for(Int_t i=0; i<=fNcomps; i++) {
    Double_t x = -fNsigmas + (2.*i*fNsigmas)/(fNcomps-1);
    (*fGaus)[i] = exp(-((x*x)/2));
    //     cout << "fGaus[" << i << "]: " << fGaus->At(i) << endl;
  }
}
//______________________________________________________________________
void AliITSresponseSDD::SetDeadChannels(Int_t nmod, Int_t nchip, Int_t nchan){
  // Set fGain to zero to simulate a random distribution of 
  // dead modules, dead chips and single dead channels

  for( Int_t m=0; m<fgkModules; m++ ) 
    for( Int_t n=0; n<fgkChips; n++ ) 
      for( Int_t p=0; p<fgkChannels; p++ ) 
	fGain[m][n][p] = 1.;
                 
  fDeadModules  = nmod;  
  fDeadChips    = nchip;  
  fDeadChannels = nchan; 
    
  // nothing to do
  if( nmod == 0 && nchip == 0 && nchan == 0 ) return;

  if( nmod < 0 || nmod > fgkModules ) 
    { 
      cout << "Wrong number of dead modules: " << nmod << endl; 
      return; 
    }
  Int_t nmax = (fgkModules-nmod)*fgkChips; 
  if( nchip < 0 || nchip > nmax ) 
    { 
      cout << "Wrong number of dead chips: " << nchip << endl; 
      return; 
    }
  nmax = ((fgkModules - nmod)*fgkChips - nchip)*fgkChannels; 
  if( nchan < 0 || nchan > nmax ) 
    { 
      cout << "Wrong number of dead channels: " << nchan << endl; 
      return; 
    }
  
  TRandom *gran = new TRandom();
  
  //  cout << "modules" << endl;
  Int_t * mod = new Int_t [nmod];
  Int_t i; //loop variable
  for( i=0; i<nmod; i++ ) 
    {
      mod[i] = (Int_t) (1.+fgkModules*gran->Uniform());
      cout << i+1 << ": Dead module nr: " << mod[i] << endl;
      for(Int_t n=0; n<fgkChips; n++)
	for(Int_t p=0; p<fgkChannels; p++)
	  fGain[mod[i]-1][n][p] = 0.;
    }

  //  cout << "chips" << endl;
  Int_t * chip     = new Int_t[nchip];
  Int_t * chipMod = new Int_t[nchip];
  i = 0;
  while( i < nchip ) 
    {
      Int_t module = (Int_t) (fgkModules*gran->Uniform() + 1.);
      if( module <=0 || module > fgkModules ) 
	cout << "Wrong module: " << module << endl;
      Int_t flagMod = 0;
      for( Int_t k=0; k<nmod; k++ ) 
	if( module == mod[k] ) { flagMod = 1; break; }
      if( flagMod == 1 ) continue;
        
      Int_t chi = (Int_t) (fgkChips*gran->Uniform() + 1.);
      if( chi <=0 || chi > fgkChips ) cout << "Wrong chip: " << chi << endl;
      i++;
      chip[i-1] = chi; 
      chipMod[i-1] = module;
      for( Int_t m=0; m<fgkChannels; m++ ) 
	fGain[module-1][chi-1][m] = 0.;
      cout << i << ": Dead chip nr. " << chip[i-1] << " in module nr: " 
	   << chipMod[i-1] << endl;
    }

  //  cout << "channels" << endl;
  Int_t * channel      = new Int_t[nchan];
  Int_t * channelChip = new Int_t[nchan];
  Int_t * channelMod  = new Int_t[nchan];
  i = 0;
  while( i < nchan ) 
    {
      Int_t k; //loop variable
      Int_t module = (Int_t) (fgkModules*gran->Uniform() + 1.);
      if( module <=0 || module > fgkModules ) 
	cout << "Wrong module: " << module << endl;
      Int_t flagMod = 0;
      for( k=0; k<nmod; k++ ) 
	if( module == mod[k] ) { flagMod = 1; break; }
      if( flagMod == 1 ) continue;
      Int_t chipp = (Int_t) (fgkChips*gran->Uniform() + 1.);
      if( chipp <=0 || chipp > fgkChips ) cout << "Wrong chip: "<< chipp<<endl;
      Int_t flagChip = 0;
      for( k=0; k<nchip; k++) 
	if( chipp == chip[k] && module == chipMod[k] ) { 
	  flagChip = 1; break; }
      if( flagChip == 1 ) continue;
      i++;
      channel[i-1] = (Int_t) (fgkChannels*gran->Uniform() + 1.); 
      if( channel[i-1] <=0 || channel[i-1] > fgkChannels ) 
	cout << "Wrong channel: " << channel[i-1] << endl;
      channelChip[i-1] = chipp;
      channelMod[i-1] = module;
      fGain[module-1][chipp-1][channel[i-1]-1] = 0.;
      cout << i << ": Dead channel nr. " << channel[i-1] << " in chip nr. " 
	   << channelChip[i-1] << " in module nr: " << channelMod[i-1] 
	   << endl;
    }
    
  delete [] mod;
  delete [] chip;
  delete [] chipMod;
  delete [] channel;
  delete [] channelMod;
  delete [] channelChip;
}
//______________________________________________________________________
void AliITSresponseSDD::PrintGains(){
  //

  if( GetDeadModules() == 0 && 
      GetDeadChips() == 0 && 
      GetDeadChannels() == 0 )
    return;  

  // Print Electronics Gains
  cout << "**************************************************" << endl; 
  cout << "             Print Electronics Gains              " << endl;
  cout << "**************************************************" << endl;

  // Print SDD electronic gains
  for(Int_t t=0; t<fgkModules;t++)
    for(Int_t u=0; u<fgkChips;u++)
      for(Int_t v=0; v<fgkChannels;v++)
	{
	  if( fGain[t][u][v] != 1.0 )
	    cout << "Gain for Module: " << t+1 << ", Chip " << u+1 << 
	      ", Channel " << v+1 << " = " << fGain[t][u][v] << endl;
	}
}
//______________________________________________________________________
void AliITSresponseSDD::Print(){
  // Print SDD response Parameters

  cout << "**************************************************" << endl;
  cout << "   Silicon Drift Detector Response Parameters    " << endl;
  cout << "**************************************************" << endl;
  cout << "Diffusion Coefficients: "<< fDiffCoeff<< ", "<<fDiffCoeff1 << endl;

  cout << "Hardware compression parameters: " << endl; 
  for(Int_t i=0; i<8; i++) cout << "fCPar[" << i << "] = " << fCPar[i] <<endl;
  cout << "Noise before electronics (arbitrary units): " << fNoise << endl;
  cout << "Baseline (ADC units): " << fBaseline << endl;
  cout << "Noise after electronics (ADC units): " << fNoiseAfterEl << endl;

  cout << "Dynamic Range: " << fDynamicRange << endl;
  cout << "Charge Loss: " << fChargeLoss << endl;
  cout << "Temperature: " << Temperature() << " K " << endl;
  cout << "Drift Speed: " << fDriftSpeed << endl;
  cout << "Electronics (1=PASCAL, 2=OLA): " << fElectronics << endl;

  cout << "N. of Sigma for signal integration: " << fNsigmas << endl;
  cout << "N. of bins in lookup table: " << fNcomps << endl;

  cout << "Max. ADC Value: " << fMaxAdc << endl;
  cout << "Min. Value: " << fMinVal << endl;

  PrintGains();

}



