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
/*
$Id$
$Log$
Revision 1.18.4.1  2002/06/10 17:51:15  hristov
Merged with v3-08-02

Revision 1.19  2002/04/24 22:02:31  nilsen
New SDigits and Digits routines, and related changes,  (including new
noise values).

*/

#include <TString.h>
#include <TRandom.h>

#include "AliITSresponseSDD.h"


//______________________________________________________________________
ClassImp(AliITSresponseSDD)

AliITSresponseSDD::AliITSresponseSDD(){
  // default constructor
   fGaus = 0;
   SetDeadChannels();
   SetMaxAdc();
   SetDiffCoeff();
   SetDriftSpeed();
   SetNSigmaIntegration();
   SetNLookUp();
   // SetClock();
   SetNoiseParam();
   SetNoiseAfterElectronics();
   SetJitterError();
   SetElectronics();
   SetDynamicRange();
   SetChargeLoss();
   SetMinVal();
   SetParamOptions();
   SetTemperature();
   SetZeroSupp();
   SetDataType();
   SetFilenames();
   SetOutputOption();
   SetDo10to8();
   // set the default zero suppression parameters
   fCPar[0]=0;
   fCPar[1]=0;
   fCPar[2]=(Int_t)(fBaseline + 2.*fNoiseAfterEl + 0.5);
   fCPar[3]=(Int_t)(fBaseline + 2.*fNoiseAfterEl + 0.5);
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
   SetMaxAdc();
   SetDiffCoeff();
   SetDriftSpeed();
   SetNSigmaIntegration();
   SetNLookUp();
   // SetClock();
   SetNoiseParam();
   SetNoiseAfterElectronics();
   SetJitterError();
   SetElectronics();
   SetDynamicRange();
   SetChargeLoss();
   SetMinVal();
   SetParamOptions();
   SetTemperature();
   SetZeroSupp();
   SetDataType(dataType);
   SetFilenames();
   SetOutputOption();
   SetDo10to8();
   // set the default zero suppression parameters
   fCPar[0]=0;
   fCPar[1]=0;
   fCPar[2]=(Int_t)(fBaseline + 2.*fNoiseAfterEl + 0.5);
   fCPar[3]=(Int_t)(fBaseline + 2.*fNoiseAfterEl + 0.5);
   fCPar[4]=0;
   fCPar[5]=0;
   fCPar[6]=0;
   fCPar[7]=0;
}
//______________________________________________________________________
AliITSresponseSDD::~AliITSresponseSDD() { 

  if(fGaus) delete fGaus;
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
void AliITSresponseSDD::GiveCompressParam(Int_t  cp[8]){
  // give compression param

    Int_t i;
    for (i=0; i<8; i++) {
        cp[i]=fCPar[i];
    } // end for i
}
//______________________________________________________________________
void AliITSresponseSDD::SetDeadChannels(Int_t nmod, Int_t nchip, Int_t nchan){
    // Set fGain to zero to simulate a random distribution of 
    // dead modules, dead chips and single dead channels

    for( Int_t m=0; m<fModules; m++ ) 
        for( Int_t n=0; n<fChips; n++ ) 
            for( Int_t p=0; p<fChannels; p++ ) 
                 fGain[m][n][p] = 1.;
                 
    fDeadModules  = nmod;  
    fDeadChips    = nchip;  
    fDeadChannels = nchan; 
    
    // nothing to do
    if( nmod == 0 && nchip == 0 && nchan == 0 ) return;

    if( nmod < 0 || nmod > fModules ) 
    { 
        cout << "Wrong number of dead modules: " << nmod << endl; 
        return; 
    }
    Int_t nmax = (fModules-nmod)*fChips; 
    if( nchip < 0 || nchip > nmax ) 
    { 
        cout << "Wrong number of dead chips: " << nchip << endl; 
        return; 
    }
    nmax = ((fModules - nmod)*fChips - nchip)*fChannels; 
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
        mod[i] = (Int_t) (1.+fModules*gran->Uniform());
        cout << i+1 << ": Dead module nr: " << mod[i] << endl;
        for(Int_t n=0; n<fChips; n++)
            for(Int_t p=0; p<fChannels; p++)
                fGain[mod[i]-1][n][p] = 0.;
    }

    //  cout << "chips" << endl;
    Int_t * chip     = new Int_t[nchip];
    Int_t * chip_mod = new Int_t[nchip];
    i = 0;
    while( i < nchip ) 
    {
        Int_t module = (Int_t) (fModules*gran->Uniform() + 1.);
        if( module <=0 || module > fModules ) 
            cout << "Wrong module: " << module << endl;
        Int_t flag_mod = 0;
        for( Int_t k=0; k<nmod; k++ ) 
            if( module == mod[k] ) { flag_mod = 1; break; }
        if( flag_mod == 1 ) continue;
        
        Int_t chi = (Int_t) (fChips*gran->Uniform() + 1.);
        if( chi <=0 || chi > fChips ) cout << "Wrong chip: " << chi << endl;
        i++;
        chip[i-1] = chi; 
        chip_mod[i-1] = module;
        for( Int_t m=0; m<fChannels; m++ ) 
            fGain[module-1][chi-1][m] = 0.;
        cout << i << ": Dead chip nr. " << chip[i-1] << " in module nr: " 
	     << chip_mod[i-1] << endl;
    }

    //  cout << "channels" << endl;
    Int_t * channel      = new Int_t[nchan];
    Int_t * channel_chip = new Int_t[nchan];
    Int_t * channel_mod  = new Int_t[nchan];
    i = 0;
    while( i < nchan ) 
    {
        Int_t k; //loop variable
        Int_t module = (Int_t) (fModules*gran->Uniform() + 1.);
        if( module <=0 || module > fModules ) 
            cout << "Wrong module: " << module << endl;
        Int_t flag_mod = 0;
        for( k=0; k<nmod; k++ ) 
            if( module == mod[k] ) { flag_mod = 1; break; }
        if( flag_mod == 1 ) continue;
        Int_t chipp = (Int_t) (fChips*gran->Uniform() + 1.);
        if( chipp <=0 || chipp > fChips ) cout << "Wrong chip: "<< chipp<<endl;
        Int_t flag_chip = 0;
        for( k=0; k<nchip; k++) 
            if( chipp == chip[k] && module == chip_mod[k] ) { 
		flag_chip = 1; break; }
        if( flag_chip == 1 ) continue;
        i++;
        channel[i-1] = (Int_t) (fChannels*gran->Uniform() + 1.); 
        if( channel[i-1] <=0 || channel[i-1] > fChannels ) 
            cout << "Wrong channel: " << channel[i-1] << endl;
        channel_chip[i-1] = chipp;
        channel_mod[i-1] = module;
        fGain[module-1][chipp-1][channel[i-1]-1] = 0.;
        cout << i << ": Dead channel nr. " << channel[i-1] << " in chip nr. " 
	     << channel_chip[i-1] << " in module nr: " << channel_mod[i-1] 
	     << endl;
    }
    
    delete [] mod;
    delete [] chip;
    delete [] chip_mod;
    delete [] channel;
    delete [] channel_mod;
    delete [] channel_chip;
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
  for(Int_t t=0; t<fModules;t++)
    for(Int_t u=0; u<fChips;u++)
      for(Int_t v=0; v<fChannels;v++)
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
   cout << "Temperature: " << fTemperature << endl;
   cout << "Drift Speed: " << fDriftSpeed << endl;
   cout << "Electronics (1=PASCAL, 2=OLA): " << fElectronics << endl;

   cout << "N. of Sigma for signal integration: " << fNsigmas << endl;
   cout << "N. of bins in lookup table: " << fNcomps << endl;

   cout << "Max. ADC Value: " << fMaxAdc << endl;
   cout << "Min. Value: " << fMinVal << endl;

   PrintGains();

}



