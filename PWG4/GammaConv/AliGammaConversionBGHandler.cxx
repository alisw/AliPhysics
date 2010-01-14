/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt                        *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class for handling of background calculation
//---------------------------------------------
////////////////////////////////////////////////

#include "AliGammaConversionBGHandler.h"
#include "AliKFParticle.h"
#include "AliAnalysisTaskGammaConversion.h"

using namespace std;

ClassImp(AliGammaConversionBGHandler)

AliGammaConversionBGHandler::AliGammaConversionBGHandler() :
  TObject(),
  fNEvents(10),
  fBGEventCounter(NULL),
  fNBinsZ(0),
  fNBinsMultiplicity(0),
  fBinLimitsArrayZ(NULL),
  fBinLimitsArrayMultiplicity(NULL),
  fBGEvents()
{
  // constructor
}

AliGammaConversionBGHandler::AliGammaConversionBGHandler(UInt_t binsZ,UInt_t binsMultiplicity,UInt_t nEvents) :
  TObject(),
  fNEvents(nEvents),
  fBGEventCounter(NULL),
  fNBinsZ(binsZ),
  fNBinsMultiplicity(binsMultiplicity),
  fBinLimitsArrayZ(NULL),
  fBinLimitsArrayMultiplicity(NULL),
  fBGEvents(binsZ,AliGammaConversionMultipicityVector(binsMultiplicity,AliGammaConversionBGEventVector(nEvents)))
{
  // constructor
}

AliGammaConversionBGHandler::AliGammaConversionBGHandler(const AliGammaConversionBGHandler & original) :
  TObject(original),
  fNEvents(original.fNEvents),
  fBGEventCounter(original.fBGEventCounter),
  fNBinsZ(original.fNBinsZ),
  fNBinsMultiplicity(original.fNBinsMultiplicity),
  fBinLimitsArrayZ(original.fBinLimitsArrayZ),
  fBinLimitsArrayMultiplicity(original.fBinLimitsArrayMultiplicity),
  fBGEvents(original.fBGEvents)
{
  //copy constructor	
}

AliGammaConversionBGHandler & AliGammaConversionBGHandler::operator = (const AliGammaConversionBGHandler & /*source*/)
{
  // assignment operator
  return *this;
}

AliGammaConversionBGHandler::~AliGammaConversionBGHandler(){

  //Kenneth remember to clear memory!!!!!!!!!!!!!!!!!!!!!
  if(fBGEventCounter){
    for(Int_t z=0;z<fNBinsZ;z++){
      delete[] fBGEventCounter[z];
    }
    delete[] fBGEventCounter;
    fBGEventCounter = NULL;
  }
  
  if(fBinLimitsArrayZ){
    delete[] fBinLimitsArrayZ;
  }

  if(fBinLimitsArrayMultiplicity){
    delete[] fBinLimitsArrayMultiplicity;
  }
}

void AliGammaConversionBGHandler::Initialize(Double_t *zBinLimitsArray, Double_t *multiplicityBinLimitsArray){
  // see header file for documantation  

  if(zBinLimitsArray){
    fBinLimitsArrayZ = zBinLimitsArray;
  }
  else{
    //Print warning
  }
  
  if(multiplicityBinLimitsArray){
    fBinLimitsArrayMultiplicity = multiplicityBinLimitsArray ;
  }
  else{
    //Print warning
  }
  fBGEventCounter= new Int_t*[fNBinsZ];
  for(Int_t z=0;z<fNBinsZ;z++){
    fBGEventCounter[z]=new Int_t[fNBinsMultiplicity];
  }
}

Int_t AliGammaConversionBGHandler::GetZBinIndex(Double_t zvalue){
  // see header file for documantation
  if(fNBinsZ<2 || zvalue<=fBinLimitsArrayZ[0]){
    return 0;
  }

  for(Int_t i=0; i<fNBinsZ-1 ;i++){
    if(zvalue >= fBinLimitsArrayZ[i] && zvalue <= fBinLimitsArrayZ[i+1]){
      return i;
    }
  }
  return fNBinsZ-1;
}

Int_t AliGammaConversionBGHandler::GetMultiplicityBinIndex(Int_t multiplicity){
  // see header file for documantation  
  if(fNBinsMultiplicity<2){
    return 0;
  }

  for(Int_t i=0; i<fNBinsMultiplicity-1 ;i++){
    if(multiplicity > fBinLimitsArrayMultiplicity[i] && multiplicity < fBinLimitsArrayMultiplicity[i+1]){
      return i;
    }
  }
  return fNBinsMultiplicity-1;
}

void AliGammaConversionBGHandler::AddEvent(TClonesArray * eventGammas, Double_t zvalue, Int_t multiplicity){
  // see header file for documantation  

  Int_t z = GetZBinIndex(zvalue);
  Int_t m = GetMultiplicityBinIndex(multiplicity);

  if(fBGEventCounter[z][m] >= fNEvents -1){
    fBGEventCounter[z][m]=0;
  }
  Int_t eventCounter=fBGEventCounter[z][m];

  //first clear the vector which is present for any gammas
  fBGEvents[z][m][eventCounter].fReconstructedGammas.clear();

  // add the gammas to the vector
  for(Int_t i=0; i< eventGammas->GetEntriesFast();i++){

    fBGEvents[z][m][eventCounter].fReconstructedGammas.push_back((AliKFParticle*)(eventGammas->At(i)));
  }

  //set the z and multiplicity numbers
  fBGEvents[z][m][eventCounter].fZVertexPosition = zvalue;
  fBGEvents[z][m][eventCounter].fChargedTrackMultiplicity = multiplicity;
}

AliGammaConversionKFVector* AliGammaConversionBGHandler::GetBGGoodV0s(Int_t event, Double_t zvalue, Int_t multiplicity){

  Int_t z = GetZBinIndex(zvalue);
  Int_t m = GetMultiplicityBinIndex(multiplicity);

  return &(fBGEvents[z][m][event].fReconstructedGammas);
}
