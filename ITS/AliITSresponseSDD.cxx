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
   // SetClock();
   SetNoiseParam();
   SetMagicValue();
   SetMinVal();
   SetParamOptions();
   SetZeroSupp();
   SetDataType();
   SetFilenames();
   SetOutputOption();
   SetDo10to8();

}

//__________________________________________________________________________
AliITSresponseSDD::AliITSresponseSDD(const AliITSresponseSDD &source){
  //     Copy Constructor 
  if(&source == this) return;
  Int_t i;
  for(i=0;i<8;i++){this->fCPar[i] = source.fCPar[i];}
  this->fNoise = source.fNoise;
  this->fBaseline = source.fBaseline;
  this->fTopValue = source.fTopValue;
  this->fTemperature = source.fTemperature;
  this->fDriftSpeed = source.fDriftSpeed;
  this->fNsigmas = source.fNsigmas;
  this->fMaxAdc = source.fMaxAdc;
  this->fDiffCoeff = source.fDiffCoeff;
  this->fDiffCoeff1 = source.fDiffCoeff1;
  this->fZeroSuppFlag = source.fZeroSuppFlag;
  this->fMinVal = source.fMinVal;
  this->fWrite = source.fWrite;
  this->fBitComp = source.fBitComp;
  this->fOption = source.fOption;
  this->fParam1 = source.fParam1;
  return;
}

//_________________________________________________________________________
AliITSresponseSDD& 
  AliITSresponseSDD::operator=(const AliITSresponseSDD &source) {
  //    Assignment operator
  if(&source == this) return *this;
  Int_t i;
  for(i=0;i<8;i++){this->fCPar[i] = source.fCPar[i];}
  this->fNoise = source.fNoise;
  this->fBaseline = source.fBaseline;
  this->fTopValue = source.fTopValue;
  this->fTemperature = source.fTemperature;
  this->fDriftSpeed = source.fDriftSpeed;
  this->fNsigmas = source.fNsigmas;
  this->fMaxAdc = source.fMaxAdc;
  this->fDiffCoeff = source.fDiffCoeff;
  this->fDiffCoeff1 = source.fDiffCoeff1;
  this->fZeroSuppFlag = source.fZeroSuppFlag;
  this->fMinVal = source.fMinVal;
  this->fWrite = source.fWrite;
  this->fBitComp = source.fBitComp;
  this->fOption = source.fOption;
  this->fParam1 = source.fParam1;
  return *this;
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
