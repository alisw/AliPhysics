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

#include "AliGammaConversionAODBGHandler.h"
#include "AliKFParticle.h"
#include "AliAODConversionPhoton.h"
#include "AliAODConversionMother.h"

using namespace std;

ClassImp(AliGammaConversionAODBGHandler)

AliGammaConversionAODBGHandler::AliGammaConversionAODBGHandler() :
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

AliGammaConversionAODBGHandler::AliGammaConversionAODBGHandler(UInt_t binsZ,UInt_t binsMultiplicity,UInt_t nEvents) :
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


AliGammaConversionAODBGHandler::AliGammaConversionAODBGHandler(UInt_t collisionSystem,UInt_t centMin,UInt_t centMax,UInt_t nEvents, Bool_t useTrackMult) :
  TObject(),
  fNEvents(nEvents),
  fBGEventCounter(NULL),
  fBGEventENegCounter(NULL),
  fBGProbability(NULL),
  fBGEventVertex(NULL),
  fNBinsZ(8),
  fNBinsMultiplicity(5),
  fBinLimitsArrayZ(NULL),
  fBinLimitsArrayMultiplicity(NULL),
  fBGEvents(fNBinsZ,AliGammaConversionMultipicityVector(fNBinsMultiplicity,AliGammaConversionBGEventVector(nEvents))),
  fBGEventsENeg(fNBinsZ,AliGammaConversionMultipicityVector(fNBinsMultiplicity,AliGammaConversionBGEventVector(nEvents)))
{
  // constructor
   
   fBinLimitsArrayZ= new Double_t[8] ;
   if(collisionSystem > 0 && collisionSystem < 8){ // PbPb
      fBinLimitsArrayZ[0] = -50.00;
      fBinLimitsArrayZ[1] = -5.5;
      fBinLimitsArrayZ[2] = -2.9;
      fBinLimitsArrayZ[3] = -0.65;
      fBinLimitsArrayZ[4] = 1.45;
      fBinLimitsArrayZ[5] = 3.65;
      fBinLimitsArrayZ[6] = 6.15;
      fBinLimitsArrayZ[7] = 50;
   }
   else if(collisionSystem == 0){
      fBinLimitsArrayZ[0] = -50.00;
      fBinLimitsArrayZ[1] = -3.375;
      fBinLimitsArrayZ[2] = -1.605;
      fBinLimitsArrayZ[3] = -0.225;
      fBinLimitsArrayZ[4] = 1.065;
      fBinLimitsArrayZ[5] = 2.445;
      fBinLimitsArrayZ[6] = 4.245;
      fBinLimitsArrayZ[7] = 50.00;
   }
   else{ 
      fBinLimitsArrayZ[0] = -50.00;
      fBinLimitsArrayZ[1] = -5.85;
      fBinLimitsArrayZ[2] = -3.35;
      fBinLimitsArrayZ[3] = -1.15;
      fBinLimitsArrayZ[4] = 0.85;
      fBinLimitsArrayZ[5] = 2.95;
      fBinLimitsArrayZ[6] = 5.55;
      fBinLimitsArrayZ[7] = 50.00;
   }



   fBinLimitsArrayMultiplicity= new Double_t[5];
   if(useTrackMult){ // pp
      fBinLimitsArrayMultiplicity[0] = 0;
      fBinLimitsArrayMultiplicity[1] = 8.5;
      fBinLimitsArrayMultiplicity[2] = 16.5;
      fBinLimitsArrayMultiplicity[3] = 27.5;
      fBinLimitsArrayMultiplicity[4] = 200.;
      if(collisionSystem > 0 && collisionSystem < 8){ // PbPb
         if(centMin == 0 && centMax == 5){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 1540.;
            fBinLimitsArrayMultiplicity[2] = 1665.;
            fBinLimitsArrayMultiplicity[3] = 1780.;
            fBinLimitsArrayMultiplicity[4] = 5000.;
         }
         else if(centMin == 0 && centMax == 10){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 1360.;
            fBinLimitsArrayMultiplicity[2] = 1520.;
            fBinLimitsArrayMultiplicity[3] = 1685.;
            fBinLimitsArrayMultiplicity[4] = 5000.;
         }
         else if(centMin == 0 && centMax == 20){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 1110.;
            fBinLimitsArrayMultiplicity[2] = 1360.;
            fBinLimitsArrayMultiplicity[3] = 1600.;
            fBinLimitsArrayMultiplicity[4] = 5000.;
         }
         else if(centMin == 0 && centMax == 80){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 890.;
            fBinLimitsArrayMultiplicity[2] = 1240.;
            fBinLimitsArrayMultiplicity[3] = 1540.;
            fBinLimitsArrayMultiplicity[4] = 5000.;
         }
         else if(centMin == 5 && centMax == 10){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 1250.;
            fBinLimitsArrayMultiplicity[2] = 1345.;
            fBinLimitsArrayMultiplicity[3] = 1445.;
            fBinLimitsArrayMultiplicity[4] = 5000.;
         }
         else if(centMin == 10 && centMax == 20){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 915.;
            fBinLimitsArrayMultiplicity[2] = 1020.;
            fBinLimitsArrayMultiplicity[3] = 1130.;
            fBinLimitsArrayMultiplicity[4] = 5000.;
         }
         else if(centMin == 20 && centMax == 40){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 510.;
            fBinLimitsArrayMultiplicity[2] = 625.;
            fBinLimitsArrayMultiplicity[3] = 730.;
            fBinLimitsArrayMultiplicity[4] = 5000.;
         }
         else if(centMin == 40 && centMax == 80){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 185.;
            fBinLimitsArrayMultiplicity[2] = 250.;
            fBinLimitsArrayMultiplicity[3] = 300.;
            fBinLimitsArrayMultiplicity[4] = 5000.;
         }
         else if(centMin == 60 && centMax == 80){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 55.;
            fBinLimitsArrayMultiplicity[2] = 80.;
            fBinLimitsArrayMultiplicity[3] = 100.;
            fBinLimitsArrayMultiplicity[4] = 5000.;
         }
         else{ // Std 20-40
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 510.;
            fBinLimitsArrayMultiplicity[2] = 625.;
            fBinLimitsArrayMultiplicity[3] = 730.;
            fBinLimitsArrayMultiplicity[4] = 5000.;
         }
      }

      else if(collisionSystem == 8 || collisionSystem == 9){ //pPb

	fBinLimitsArrayMultiplicity[0] = 0.;
	fBinLimitsArrayMultiplicity[1] = 7.5;
	fBinLimitsArrayMultiplicity[2] = 16.5;
	fBinLimitsArrayMultiplicity[3] = 29.5;
	fBinLimitsArrayMultiplicity[4] = 500.;	
	
	if(centMin == 0 && centMax == 20){
	  fBinLimitsArrayMultiplicity[0] = 0.;
	  fBinLimitsArrayMultiplicity[1] = 31.5;
	  fBinLimitsArrayMultiplicity[2] = 40.5;
	  fBinLimitsArrayMultiplicity[3] = 50.5;
	  fBinLimitsArrayMultiplicity[4] = 500.;
         }	
	else if(centMin == 20 && centMax == 40){
	  fBinLimitsArrayMultiplicity[0] = 0.;
	  fBinLimitsArrayMultiplicity[1] = 19.5;
	  fBinLimitsArrayMultiplicity[2] = 25.5;
	  fBinLimitsArrayMultiplicity[3] = 32.5;
	  fBinLimitsArrayMultiplicity[4] = 500.;
         }	
	else if(centMin == 40 && centMax == 60){
	  fBinLimitsArrayMultiplicity[0] = 0.;
	  fBinLimitsArrayMultiplicity[1] = 12.5;
	  fBinLimitsArrayMultiplicity[2] = 16.5;
	  fBinLimitsArrayMultiplicity[3] = 22.5;
	  fBinLimitsArrayMultiplicity[4] = 500.;
         }	
	else if(centMin == 60 && centMax == 80){
	  fBinLimitsArrayMultiplicity[0] = 0.;
	  fBinLimitsArrayMultiplicity[1] = 5.5;
	  fBinLimitsArrayMultiplicity[2] = 9.5;
	  fBinLimitsArrayMultiplicity[3] = 13.5;
	  fBinLimitsArrayMultiplicity[4] = 500.;
         }	
      }
   }
   else{// pp or pPb V0 Mult
      fBinLimitsArrayMultiplicity[0] = 2;
      fBinLimitsArrayMultiplicity[1] = 3;
      fBinLimitsArrayMultiplicity[2] = 4;
      fBinLimitsArrayMultiplicity[3] = 5;
      fBinLimitsArrayMultiplicity[4] = 9999;
      if(collisionSystem > 0 && collisionSystem < 8){ // PbPb
         if(centMin == 0 && centMax == 5){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 27.;
            fBinLimitsArrayMultiplicity[2] = 31.;
            fBinLimitsArrayMultiplicity[3] = 36.;
            fBinLimitsArrayMultiplicity[4] = 100.;
         }
         else if(centMin == 0 && centMax == 10){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 25.;
            fBinLimitsArrayMultiplicity[2] = 30.;
            fBinLimitsArrayMultiplicity[3] = 36.;
            fBinLimitsArrayMultiplicity[4] = 100.;
         }
         else if(centMin == 0 && centMax == 20){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 22.;
            fBinLimitsArrayMultiplicity[2] = 27.;
            fBinLimitsArrayMultiplicity[3] = 33.;
            fBinLimitsArrayMultiplicity[4] = 100.;
         }
         else if(centMin == 0 && centMax == 80){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 18.;
            fBinLimitsArrayMultiplicity[2] = 25.;
            fBinLimitsArrayMultiplicity[3] = 32.;
            fBinLimitsArrayMultiplicity[4] = 100.;
         }
         else if(centMin == 5 && centMax == 10){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 23.;
            fBinLimitsArrayMultiplicity[2] = 27.;
            fBinLimitsArrayMultiplicity[3] = 32.;
            fBinLimitsArrayMultiplicity[4] = 100.;
         }
         else if(centMin == 10 && centMax == 20){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 18.;
            fBinLimitsArrayMultiplicity[2] = 22.;
            fBinLimitsArrayMultiplicity[3] = 27.;
            fBinLimitsArrayMultiplicity[4] = 100.;
         }
         else if(centMin == 20 && centMax == 40){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 11.;
            fBinLimitsArrayMultiplicity[2] = 14.;
            fBinLimitsArrayMultiplicity[3] = 18.;
            fBinLimitsArrayMultiplicity[4] = 100.;
         }
         else if(centMin == 40 && centMax == 80){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 5.;
            fBinLimitsArrayMultiplicity[2] = 7.;
            fBinLimitsArrayMultiplicity[3] = 11.;
            fBinLimitsArrayMultiplicity[4] = 100.;
         }
         else if(centMin == 60 && centMax == 80){
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 2.;
            fBinLimitsArrayMultiplicity[2] = 3.;
            fBinLimitsArrayMultiplicity[3] = 5.;
            fBinLimitsArrayMultiplicity[4] = 100.;
         }
         else{ // Std 20-40
            fBinLimitsArrayMultiplicity[0] = 0.;
            fBinLimitsArrayMultiplicity[1] = 11.;
            fBinLimitsArrayMultiplicity[2] = 14.;
            fBinLimitsArrayMultiplicity[3] = 18.;
            fBinLimitsArrayMultiplicity[4] = 100.;
         }
      }
   }
   
   Initialize(fBinLimitsArrayZ,fBinLimitsArrayMultiplicity);

}

AliGammaConversionAODBGHandler::AliGammaConversionAODBGHandler(const AliGammaConversionAODBGHandler & original) :
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

AliGammaConversionAODBGHandler & AliGammaConversionAODBGHandler::operator = (const AliGammaConversionAODBGHandler & /*source*/)
{
  // assignment operator
  return *this;
}

AliGammaConversionAODBGHandler::~AliGammaConversionAODBGHandler(){

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

void AliGammaConversionAODBGHandler::Initialize(Double_t * const zBinLimitsArray, Double_t * const multiplicityBinLimitsArray){
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

Int_t AliGammaConversionAODBGHandler::GetZBinIndex(Double_t zvalue) const{
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

Int_t AliGammaConversionAODBGHandler::GetMultiplicityBinIndex(Int_t multiplicity) const{
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

void AliGammaConversionAODBGHandler::AddEvent(TList* const eventGammas,Double_t xvalue, Double_t yvalue, Double_t zvalue, Int_t multiplicity){
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
    delete (AliAODConversionPhoton*)(fBGEvents[z][m][eventCounter][d]);
  }
  fBGEvents[z][m][eventCounter].clear();
  
  // add the gammas to the vector
  for(Int_t i=0; i< eventGammas->GetEntries();i++){
    //    AliKFParticle *t = new AliKFParticle(*(AliKFParticle*)(eventGammas->At(i)));
    fBGEvents[z][m][eventCounter].push_back(new AliAODConversionPhoton(*(AliAODConversionPhoton*)(eventGammas->At(i))));
  }
  fBGEventCounter[z][m]++;
}
void AliGammaConversionAODBGHandler::AddElectronEvent(TClonesArray* const eventENeg, Double_t zvalue, Int_t multiplicity){

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
    delete (AliAODConversionPhoton*)(fBGEventsENeg[z][m][eventENegCounter][d]);
  }

  fBGEventsENeg[z][m][eventENegCounter].clear();

  // add the electron to the vector
  for(Int_t i=0; i< eventENeg->GetEntriesFast();i++){
    //    AliKFParticle *t = new AliKFParticle(*(AliKFParticle*)(eventGammas->At(i)));
    fBGEventsENeg[z][m][eventENegCounter].push_back(new AliAODConversionPhoton(*(AliAODConversionPhoton*)(eventENeg->At(i))));
  }

  fBGEventENegCounter[z][m]++;


}
AliGammaConversionAODVector* AliGammaConversionAODBGHandler::GetBGGoodV0s(Int_t zbin, Int_t mbin, Int_t event){
  //see headerfile for documentation
  return &(fBGEvents[zbin][mbin][event]);
}
AliGammaConversionAODVector* AliGammaConversionAODBGHandler::GetBGGoodENeg(Int_t event, Double_t zvalue, Int_t multiplicity){


  //see headerfile for documentation
  Int_t z = GetZBinIndex(zvalue);
  Int_t m = GetMultiplicityBinIndex(multiplicity);
  return &(fBGEventsENeg[z][m][event]);


}
void AliGammaConversionAODBGHandler::PrintBGArray(){
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
