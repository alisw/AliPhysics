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


/////////////////////////////////////////////////////////////////
//							       //
//    Class containing the parameters that are configured      //
//    to trigger events with the ZDC (in A-A collisions)       //
//    Use: store the set of parameters needed to calculate     //
//    the trigger function sent to the CTP                     //
//							       //
//    Author: Chiara.Oppedisano@to.infn.it                     //
//							       //
/////////////////////////////////////////////////////////////////                                                             

#include "AliZDCTriggerParameters.h"

ClassImp(AliZDCTriggerParameters)

//________________________________________________________________
AliZDCTriggerParameters::AliZDCTriggerParameters() :
TObject(),
fADCZEMCentralityThr(0),
fADCMBThreshold(0),
fDiscZEMCentralityThr(0),
fDiscMBThreshold(0)
{
  // Default constructor
  for(Int_t j=0; j<4; j++){
     fADCEMDWindow[j] = fDiscEMDWindow[j] = 0.;
     if(j<2){
       fADCCentralWindow[j] = fADCSemicentralWindow[j] = 0.;
       fDiscCentralWindow[j] = fDiscSemicentralWindow[j] = 0.;
     }
  }
}  

//________________________________________________________________
AliZDCTriggerParameters::AliZDCTriggerParameters(Float_t *adcParam, 
	Float_t *discParam) :
fADCZEMCentralityThr(adcParam[0]),
fADCMBThreshold(adcParam[1]),
fDiscZEMCentralityThr(discParam[0]),
fDiscMBThreshold(discParam[1])
{
  // Standard constructor
  fADCCentralWindow[0] = adcParam[2];
  fADCCentralWindow[1] = adcParam[3];
  fADCSemicentralWindow[0] = adcParam[4];
  fADCSemicentralWindow[1] = adcParam[5];
  fADCEMDWindow[0] = adcParam[6];
  fADCEMDWindow[1] = adcParam[7];
  fADCEMDWindow[2] = adcParam[8];
  fADCEMDWindow[3] = adcParam[9];
  //
  fDiscCentralWindow[0] = discParam[2];
  fDiscCentralWindow[1] = discParam[3];
  fDiscSemicentralWindow[0] = discParam[4];
  fDiscSemicentralWindow[1] = discParam[5];
  fDiscEMDWindow[0] = discParam[6];
  fDiscEMDWindow[1] = discParam[7];
  fDiscEMDWindow[2] = discParam[8];
  fDiscEMDWindow[3] = discParam[9];
}

//____________________________________________________________________________
AliZDCTriggerParameters::AliZDCTriggerParameters(const AliZDCTriggerParameters& oldTrigPar) :
  TObject(),
  fADCZEMCentralityThr(oldTrigPar.fADCZEMCentralityThr),
  fADCMBThreshold(oldTrigPar.fADCMBThreshold),
  fDiscZEMCentralityThr(oldTrigPar.fDiscZEMCentralityThr),
  fDiscMBThreshold(oldTrigPar.fDiscMBThreshold)
{
  // Copy constructor
  fADCCentralWindow[0] = oldTrigPar.fADCCentralWindow[0];
  fADCCentralWindow[1] = oldTrigPar.fADCCentralWindow[1];
  fADCSemicentralWindow[0] = oldTrigPar.fADCSemicentralWindow[0];
  fADCSemicentralWindow[1] = oldTrigPar.fADCSemicentralWindow[1];
  fADCEMDWindow[0] = oldTrigPar.fADCEMDWindow[0];
  fADCEMDWindow[1] = oldTrigPar.fADCEMDWindow[1];
  fADCEMDWindow[2] = oldTrigPar.fADCEMDWindow[2];
  fADCEMDWindow[3] = oldTrigPar.fADCEMDWindow[3];
  //
  fDiscCentralWindow[0] = oldTrigPar.fDiscCentralWindow[0];
  fDiscCentralWindow[1] = oldTrigPar.fDiscCentralWindow[1];
  fDiscSemicentralWindow[0] = oldTrigPar.fDiscSemicentralWindow[0];
  fDiscSemicentralWindow[1] = oldTrigPar.fDiscSemicentralWindow[1];
  fDiscEMDWindow[0] = oldTrigPar.fDiscEMDWindow[0];
  fDiscEMDWindow[1] = oldTrigPar.fDiscEMDWindow[1];
  fDiscEMDWindow[2] = oldTrigPar.fDiscEMDWindow[2];
  fDiscEMDWindow[3] = oldTrigPar.fDiscEMDWindow[3];
  
}

//____________________________________________________________________________
AliZDCTriggerParameters &AliZDCTriggerParameters::operator= (const AliZDCTriggerParameters& param) 
{
  //assignement operator
  
  if(&param == this) return *this;

  fADCZEMCentralityThr = param.fADCZEMCentralityThr;
  fADCMBThreshold = param.fADCMBThreshold;
  fDiscZEMCentralityThr = param.fDiscZEMCentralityThr;
  fDiscMBThreshold = param.fDiscMBThreshold;
  for(int i=0; i<2; i++) {
  	fADCCentralWindow[i] = param.fADCCentralWindow[i];
  	fADCSemicentralWindow[i] = param.fADCSemicentralWindow[i];
  	fDiscCentralWindow[i] = param.fDiscCentralWindow[i];
  	fDiscSemicentralWindow[i] = param.fDiscSemicentralWindow[i];
 }
 for(int j=0; j<4; j++){
 	fADCEMDWindow[j] = param.fADCEMDWindow[j];
  	fDiscEMDWindow[j] = param.fDiscEMDWindow[j];
 }
 
 return *this;
}
