/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include <TString.h>

#include "AliITSresponseSDD.h"


//___________________________________________
ClassImp(AliITSresponseSDD)	

AliITSresponseSDD::AliITSresponseSDD()
{
  // constructor
   SetMaxAdc();
   SetDiffCoeff();
   SetDriftSpeed();
   SetNSigmaIntegration();
   SetNLookUp();
   // SetClock();
   SetNoiseParam();
   SetNoiseAfterElectronics();
   SetElectronics();
   SetDynamicRange();
   SetChargeLoss();
   SetMinVal();
   SetParamOptions();
   SetZeroSupp();
   SetDataType();
   SetFilenames();
   SetOutputOption();
   SetDo10to8();
}
AliITSresponseSDD::~AliITSresponseSDD() { 

  if(fGaus) delete fGaus;

}

void AliITSresponseSDD::SetCompressParam(Int_t  cp[8])
{
  // set compression param

    Int_t i;
    for (i=0; i<8; i++) {
	fCPar[i]=cp[i];
	//printf("\n CompressPar %d %d \n",i,fCPar[i]);
	
    }
}
void AliITSresponseSDD::GiveCompressParam(Int_t  cp[8])
{
  // give compression param

    Int_t i;
    for (i=0; i<8; i++) {
	cp[i]=fCPar[i];
    }
}

void AliITSresponseSDD::Print()
{
  // Print SDD response Parameters

   cout << "**************************************************" << endl;
   cout << "   Silicon Drift Detector Response Parameters    " << endl;
   cout << "**************************************************" << endl;
   cout << "Diffusion Coefficients: " << fDiffCoeff << ", " << fDiffCoeff1 << endl;

   cout << "Hardware compression parameters: " << endl; 
   for(Int_t i=0; i<8; i++) cout << "fCPar[" << i << "] = " << fCPar[i] << endl;
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

   cout << "**************************************************" << endl;
  


}



