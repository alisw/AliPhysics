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
  fBGEventENegCounter(NULL),
  fBGProbability(NULL),
  fBGEventVertex(NULL),
  fNBinsZ(0),
  fNBinsMultiplicity(0),
  fBinLimitsArrayZ(NULL),
  fBinLimitsArrayMultiplicity(NULL),
  fBGEvents(),
  fBGEventsENeg()
{
  // constructor
}

AliGammaConversionBGHandler::AliGammaConversionBGHandler(UInt_t binsZ,UInt_t binsMultiplicity,UInt_t nEvents) :
  TObject(),
  fNEvents(nEvents),
  fBGEventCounter(NULL),
  fBGEventENegCounter(NULL),
  fBGProbability(NULL),
  fBGEventVertex(NULL),
  fNBinsZ(binsZ),
  fNBinsMultiplicity(binsMultiplicity),
  fBinLimitsArrayZ(NULL),
  fBinLimitsArrayMultiplicity(NULL),
  fBGEvents(binsZ,AliGammaConversionMultipicityVector(binsMultiplicity,AliGammaConversionBGEventVector(nEvents))),
  fBGEventsENeg(binsZ,AliGammaConversionMultipicityVector(binsMultiplicity,AliGammaConversionBGEventVector(nEvents)))
{
  // constructor
}

AliGammaConversionBGHandler::AliGammaConversionBGHandler(const AliGammaConversionBGHandler & original) :
  TObject(original),
  fNEvents(original.fNEvents),
  fBGEventCounter(original.fBGEventCounter),
  fBGEventENegCounter(original.fBGEventENegCounter),
  fBGProbability(original.fBGProbability),
  fBGEventVertex(original.fBGEventVertex),
  fNBinsZ(original.fNBinsZ),
  fNBinsMultiplicity(original.fNBinsMultiplicity),
  fBinLimitsArrayZ(original.fBinLimitsArrayZ),
  fBinLimitsArrayMultiplicity(original.fBinLimitsArrayMultiplicity),
  fBGEvents(original.fBGEvents),
  fBGEventsENeg(original.fBGEventsENeg)
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

  if(fBGEventVertex){
    for(Int_t z=0;z<fNBinsZ;z++){
      for(Int_t m=0;m<fNBinsMultiplicity;m++){
	delete [] fBGEventVertex[z][m];
      }
      delete [] fBGEventVertex[z];
    }
    delete [] fBGEventVertex;
  }

   if(fBGEventENegCounter){
    for(Int_t z=0;z<fNBinsZ;z++){
      delete[] fBGEventENegCounter[z];
    }
    delete[] fBGEventENegCounter;
    fBGEventENegCounter = NULL;
  }

  if(fBinLimitsArrayZ){
    delete[] fBinLimitsArrayZ;
  }

  if(fBinLimitsArrayMultiplicity){
    delete[] fBinLimitsArrayMultiplicity;
  }
}

void AliGammaConversionBGHandler::Initialize(Double_t * const zBinLimitsArray, Double_t * const multiplicityBinLimitsArray){
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
  if(fBGEventCounter == NULL){
    fBGEventCounter= new Int_t*[fNBinsZ];
  }
  for(Int_t z=0;z<fNBinsZ;z++){
    fBGEventCounter[z]=new Int_t[fNBinsMultiplicity];
  }

  for(Int_t z=0;z<fNBinsZ;z++){
    for(Int_t m=0;m<fNBinsMultiplicity;m++){
      fBGEventCounter[z][m]=0;
    }
  }

  if(fBGEventVertex == NULL){
    fBGEventVertex = new GammaConversionVertex**[fNBinsZ];
  }
  for(Int_t z=0; z < fNBinsZ; z++){
    fBGEventVertex[z]= new GammaConversionVertex*[fNBinsMultiplicity];
  }
  for(Int_t z=0;z<fNBinsZ;z++){
    for(Int_t m=0;m<fNBinsMultiplicity; m++){
      fBGEventVertex[z][m]= new GammaConversionVertex[fNEvents];
    }
  }
  if( fBGEventENegCounter == NULL){
    fBGEventENegCounter = new Int_t*[fNBinsZ];
  }

  for(Int_t z=0; z < fNBinsZ; z++){
    fBGEventENegCounter[z] = new Int_t[fNBinsMultiplicity];
  }

  for(Int_t z=0;z<fNBinsZ;z++){
    for(Int_t m=0;m<fNBinsMultiplicity; m++){
      fBGEventENegCounter[z][m] = 0;
    }
  }

  if(fBGProbability == NULL){
    fBGProbability = new Double_t*[fNBinsZ];
  }
  for(Int_t z=0; z < fNBinsZ; z++){
    fBGProbability[z] = new Double_t[fNBinsMultiplicity];
  }

  for(Int_t z=0;z<fNBinsZ;z++){
    for(Int_t m=0;m<fNBinsMultiplicity; m++){
      fBGProbability[z][m] = 0;
    }
  }
  //filling the probability
  fBGProbability[0][0] = 0.243594;
  fBGProbability[0][1] = 0.279477;
  fBGProbability[0][2] = 0.305104;
  fBGProbability[0][3] = 0.315927;
  fBGProbability[1][0] = 0.241964;
  fBGProbability[1][1] = 0.272995;
  fBGProbability[1][2] = 0.307165;
  fBGProbability[1][3] = 0.292248;
  fBGProbability[2][0] = 0.241059;
  fBGProbability[2][1] = 0.27509;
  fBGProbability[2][2] = 0.283657;
  fBGProbability[2][3] = 0.310512;
  fBGProbability[3][0] = 0.23888;
  fBGProbability[3][1] = 0.283418;
  fBGProbability[3][2] = 0.297232;
  fBGProbability[3][3] = 0.348188;
  fBGProbability[4][0] = 0.245555;
  fBGProbability[4][1] = 0.281218;
  fBGProbability[4][2] = 0.317236;
  fBGProbability[4][3] = 0.323495;
  fBGProbability[5][0] = 0.244572;
  fBGProbability[5][1] = 0.259498;
  fBGProbability[5][2] = 0.278383;
  fBGProbability[5][3] = 0.284696;
  fBGProbability[6][0] = 0.24703;
  fBGProbability[6][1] = 0.275265;
  fBGProbability[6][2] = 0.284004;
  fBGProbability[6][3] = 0.343584;
 
}

Int_t AliGammaConversionBGHandler::GetZBinIndex(Double_t zvalue) const{
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

Int_t AliGammaConversionBGHandler::GetMultiplicityBinIndex(Int_t multiplicity) const{
  // see header file for documantation  
  if(fNBinsMultiplicity<2){
    return 0;
  }

  for(Int_t i=0; i<fNBinsMultiplicity-1 ;i++){
    if(multiplicity >= fBinLimitsArrayMultiplicity[i] && multiplicity < fBinLimitsArrayMultiplicity[i+1]){
      return i;
    }
  }
  return fNBinsMultiplicity-1;
}

void AliGammaConversionBGHandler::AddEvent(TClonesArray * const eventGammas,Double_t xvalue, Double_t yvalue, Double_t zvalue, Int_t multiplicity){
  // see header file for documantation  

  //  cout<<"Entering the AddEvent function"<<endl;

  Int_t z = GetZBinIndex(zvalue);
  Int_t m = GetMultiplicityBinIndex(multiplicity);

  if(fBGEventCounter[z][m] >= fNEvents){
    fBGEventCounter[z][m]=0;
  }
  Int_t eventCounter=fBGEventCounter[z][m];
  
  /*
  if(fBGEventVertex[z][m][eventCounter]){
    delete fBGEventVertex[z][m][eventCounter];
  }
  */
  fBGEventVertex[z][m][eventCounter].fX = xvalue;
  fBGEventVertex[z][m][eventCounter].fY = yvalue;
  fBGEventVertex[z][m][eventCounter].fZ = zvalue;

  //first clear the vector
  // cout<<"Size of vector: "<<fBGEvents[z][m][eventCounter].size()<<endl;
  //  cout<<"Checking the entries: Z="<<z<<", M="<<m<<", eventCounter="<<eventCounter<<endl;

  //  cout<<"The size of this vector is: "<<fBGEvents[z][m][eventCounter].size()<<endl;
  for(UInt_t d=0;d<fBGEvents[z][m][eventCounter].size();d++){
    delete (AliKFParticle*)(fBGEvents[z][m][eventCounter][d]);
  }
  fBGEvents[z][m][eventCounter].clear();
  
  // add the gammas to the vector
  for(Int_t i=0; i< eventGammas->GetEntriesFast();i++){
    //    AliKFParticle *t = new AliKFParticle(*(AliKFParticle*)(eventGammas->At(i)));
    fBGEvents[z][m][eventCounter].push_back(new AliKFParticle(*(AliKFParticle*)(eventGammas->At(i))));
  }
  fBGEventCounter[z][m]++;
}
void AliGammaConversionBGHandler::AddElectronEvent(TClonesArray* const eventENeg, Double_t zvalue, Int_t multiplicity){

  Int_t z = GetZBinIndex(zvalue);
  Int_t m = GetMultiplicityBinIndex(multiplicity);

  if(fBGEventENegCounter[z][m] >= fNEvents){
     fBGEventENegCounter[z][m]=0;
  }
 

   Int_t eventENegCounter=fBGEventENegCounter[z][m];
  
  //first clear the vector
  // cout<<"Size of vector: "<<fBGEvents[z][m][eventCounter].size()<<endl;
  //  cout<<"Checking the entries: Z="<<z<<", M="<<m<<", eventCounter="<<eventCounter<<endl;

  //  cout<<"The size of this vector is: "<<fBGEvents[z][m][eventCounter].size()<<endl;
  for(UInt_t d=0;d<fBGEventsENeg[z][m][eventENegCounter].size();d++){
    delete (AliKFParticle*)(fBGEventsENeg[z][m][eventENegCounter][d]);
  }

  fBGEventsENeg[z][m][eventENegCounter].clear();

  // add the electron to the vector
  for(Int_t i=0; i< eventENeg->GetEntriesFast();i++){
    //    AliKFParticle *t = new AliKFParticle(*(AliKFParticle*)(eventGammas->At(i)));
    fBGEventsENeg[z][m][eventENegCounter].push_back(new AliKFParticle(*(AliKFParticle*)(eventENeg->At(i))));
  }

  fBGEventENegCounter[z][m]++;


}
AliGammaConversionKFVector* AliGammaConversionBGHandler::GetBGGoodV0s(Int_t zbin, Int_t mbin, Int_t event){
  //see headerfile for documentation
  return &(fBGEvents[zbin][mbin][event]);
}
AliGammaConversionKFVector* AliGammaConversionBGHandler::GetBGGoodENeg(Int_t event, Double_t zvalue, Int_t multiplicity){


  //see headerfile for documentation
  Int_t z = GetZBinIndex(zvalue);
  Int_t m = GetMultiplicityBinIndex(multiplicity);
  return &(fBGEventsENeg[z][m][event]);


}
void AliGammaConversionBGHandler::PrintBGArray(){
  //see headerfile for documentation
  for(Int_t z=0;z<fNBinsZ;z++){
    if(z==2){
      cout<<"Getting the data for z bin: "<<z<<endl;
      for(Int_t multiplicity=0;multiplicity<fNBinsMultiplicity;multiplicity++){
	if(multiplicity==2){
	  cout<<"Getting the data for multiplicity bin: "<<multiplicity<<endl;	
	  for(Int_t event=0;event<fNEvents;event++){
	    if(fBGEvents[z][multiplicity][event].size()>0){
	      cout<<"Event: "<<event<<" has: "<<fBGEvents[z][multiplicity][event].size()<<endl;
	    }
	  }
	}
      }
    }
  }
}
