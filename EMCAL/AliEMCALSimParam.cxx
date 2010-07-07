/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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
 //
 // Base class for the EMCAL simulation parameters.
 //
 //
 */

// --- Root header files ---
#include "TMath.h"
// --- AliRoot header files ---
#include "AliEMCALSimParam.h"
#include "AliLog.h"


ClassImp(AliEMCALSimParam)

AliEMCALSimParam  * AliEMCALSimParam::fgSimParam = 0 ;
//-----------------------------------------------------------------------------
AliEMCALSimParam::AliEMCALSimParam() :
TNamed(),    
fDigitThreshold(0),
fMeanPhotonElectron(0),
fPinNoise(0),
fTimeDelay(0),
fTimeResolution(0),
//fTimeThreshold(0),    
//fTimeSignalLength(0),
fNADCEC(0),//Digitizer
fA(0.),
fB(0.),
fECPrimThreshold(0.) //SDigitizer   
{
	//Constructor 
	
	//Parameters in Digitizer
	fMeanPhotonElectron = 4400;  // electrons per GeV 
	fPinNoise           = 0.012; // pin noise in GeV from analysis test beam data 
	fDigitThreshold     = 3; // 3 ADC counts not anymore cut in energy: //fPinNoise * 3; // 3 * sigma
	fTimeResolution     = 0.6e-9 ; // 600 pc
	fTimeDelay          = 600e-9 ; // 600 nc

	//fTimeSignalLength   = 1.0e-9 ;
	fNADCEC             = (Int_t) TMath::Power(2,16) ; // number of channels in Tower ADC - 65536
	//fTimeThreshold      = 0.001*10000000 ; // Means 1 MeV in terms of SDigits amplitude ??
	
	//SDigitizer
	fA                  = 0;
	fB                  = 1.e+6; // Dynamic range now 2 TeV
	fECPrimThreshold    = 0.05;  // GeV	// threshold for deposit energy of hit
	
}


//-----------------------------------------------------------------------------
AliEMCALSimParam::AliEMCALSimParam(const AliEMCALSimParam& ):
TNamed(),
fDigitThreshold(0),
fMeanPhotonElectron(0),
fPinNoise(0),
fTimeDelay(0),
fTimeResolution(0),
//fTimeThreshold(0),    
//fTimeSignalLength(0),//Digitizer
fNADCEC(0),
fA(0.),
fB(0.),
fECPrimThreshold(0.)//SDigitizer
{
  //Copy constructor.
  AliError("Should not use copy constructor for singleton") ;

  fgSimParam = this ;
	
}

//-----------------------------------------------------------------------------                                                            
AliEMCALSimParam * AliEMCALSimParam::GetInstance(){
// Get Instance

	if(!fgSimParam){
		fgSimParam = new AliEMCALSimParam() ;
	}
	
	return fgSimParam ;
	
}

//-----------------------------------------------------------------------------
AliEMCALSimParam& AliEMCALSimParam::operator = (const AliEMCALSimParam& simParam)
{
  //Assignment operator.

  if(this != &simParam) {
    AliError("Should not use operator= for singleton\n") ;
  }

  return *this;
}

//-----------------------------------------------------------------------------
void AliEMCALSimParam::Print(Option_t *) const
{
	// Print simulation parameters to stdout
	
	printf("=== Parameters in Digitizer === \n");
	printf("\t Electronics noise in EMC (fPinNoise)       = %f\n", fPinNoise) ;
	printf("\t Threshold  in EMC  (fDigitThreshold)       = %d\n", fDigitThreshold)  ;
	printf("\t Time Resolution (fTimeResolution)          = %g\n", fTimeResolution) ;
	printf("\t Time Delay (fTimeDelay)                    = %g\n", fTimeDelay) ;
	printf("\t Mean Photon-Electron (fMeanPhotonElectron) = %d\n", fMeanPhotonElectron)  ;
	printf("\t N channels in EC section ADC (fNADCEC)     = %d\n", fNADCEC) ;

	printf("\n");
	
	printf("=== Parameters in SDigitizer === \n");
	printf("\t sdigitization parameters       A = %f\n",     fA);
	printf("\t                                B = %f\n",     fB);
	printf("\t Threshold for EC Primary assignment  = %f\n", fECPrimThreshold);
	
}

