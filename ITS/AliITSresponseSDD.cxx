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

#include <TMath.h>

#include "AliITSgeom.h"
#include "AliITSresponseSDD.h"
#include "AliITS.h"
#include "AliRun.h"

class AliITS;

//___________________________________________
ClassImp(AliITSresponseSDD)	

AliITSresponseSDD::AliITSresponseSDD()
{
  // constructor
   SetMaxAdc();
   SetDiffCoeff();
   SetQref();
   SetDriftSpeed();
   // SetClock();
   SetNoiseParam();
   SetMagicValue();
   SetMinVal();
   SetParamOptions();
   SetZeroSupp();
   SetDataType();
   SetFilenames();
   SetOutputOption();

}

//__________________________________________________________________________
AliITSresponseSDD::AliITSresponseSDD(const AliITSresponseSDD &source){
  //     Copy Constructor 
  Int_t i;
  if(&source == this) return;
  for(i=0,i<8,i++){this->fCPar[i] = source.fCPar[i];}
  this->fNoise = source.fNoise;
  this->fBaseline = source.fBaseline;
  this->fTopValue = source.fTopValue;
  this->fTemperature = source.fTemperature;
  this->fDriftSpeed = source.fDriftSpeed;
  this->fMaxAdc = source.fMaxAdc;
  this->fDiffCoeff = source.fDiffCoeff;
  this->fQref = source.fQref;
  this->fZeroSuppFlag = source.fZeroSuppFlag;
  this->fMinVal = source.fMinVal;
  this->fWrite = source.fWrite;
  this->fOption = source.fOption;
  this->fParam1 = source.fParam1;
  return;
}

//_________________________________________________________________________
AliITSresponseSDD& 
  AliITSresponseSDD::operator=(const AliITSresponseSDD &source) {
  //    Assignment operator
  Int_t i;
  if(&source == this) return *this;
  for(i=0,i<8,i++){this->fCPar[i] = source.fCPar[i];}
  this->fNoise = source.fNoise;
  this->fBaseline = source.fBaseline;
  this->fTopValue = source.fTopValue;
  this->fTemperature = source.fTemperature;
  this->fDriftSpeed = source.fDriftSpeed;
  this->fMaxAdc = source.fMaxAdc;
  this->fDiffCoeff = source.fDiffCoeff;
  this->fQref = source.fQref;
  this->fZeroSuppFlag = source.fZeroSuppFlag;
  this->fMinVal = source.fMinVal;
  this->fWrite = source.fWrite;
  this->fOption = source.fOption;
  this->fParam1 = source.fParam1;
  return *this;
}

void AliITSresponseSDD::SetCompressParam(Int_t  cp[8])
{
  // set compression param
    Int_t i;
    for(i=0; i<8; i++) {
	fCPar[i]=cp[i];
	//printf("\n CompressPar %d %d \n",i,fCPar[i]);
	
    }
}
void AliITSresponseSDD::GiveCompressParam(Int_t  cp[8])
{
  // give compression param
    Int_t i;
    for(i=0; i<8; i++) {
	cp[i]=fCPar[i];
    }
}
//______________________________________________________________________________
void AliITSresponseSDD::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliITSresponseSDD.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliITSresponse::Streamer(R__b);
      R__b.ReadStaticArray(fCPar);
      R__b >> fNoise;
      R__b >> fBaseline;
      R__b >> fTopValue;
      R__b >> fTemperature;
      R__b >> fDriftSpeed;
      R__b >> fMaxAdc;
      R__b >> fDiffCoeff;
      R__b >> fQref;
      R__b >> fZeroSuppFlag;
      R__b >> fMinVal;
      R__b >> fWrite;
      //R__b.ReadArray(fOption); // Not to be printed out?
      //R__b.ReadArray(fParam1); // Not to be printed out?
      //R__b.ReadArray(fParam2); // Not to be printed out?
      fDataType.Streamer(R__b);
      fFileName1.Streamer(R__b);
      fFileName2.Streamer(R__b);
      fFileName3.Streamer(R__b);
   } else {
      R__b.WriteVersion(AliITSresponseSDD::IsA());
      AliITSresponse::Streamer(R__b);
      R__b.WriteArray(fCPar, 8);
      R__b << fNoise;
      R__b << fBaseline;
      R__b << fTopValue;
      R__b << fTemperature;
      R__b << fDriftSpeed;
      R__b << fMaxAdc;
      R__b << fDiffCoeff;
      R__b << fQref;
      R__b << fZeroSuppFlag;
      R__b << fMinVal;
      R__b << fWrite;
      //R__b.WriteArray(fOption, __COUNTER__); // Not to be printed out?
      //R__b.WriteArray(fParam1, __COUNTER__); // Not to be printed out?
      //R__b.WriteArray(fParam2, __COUNTER__); // Not to be printed out?
      fDataType.Streamer(R__b);
      fFileName1.Streamer(R__b);
      fFileName2.Streamer(R__b);
      fFileName3.Streamer(R__b);
   }
}
