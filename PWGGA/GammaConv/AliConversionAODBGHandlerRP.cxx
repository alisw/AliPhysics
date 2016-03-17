/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                           *
* Authors: Daniel Lohner (Daniel.Lohner@cern.ch), *
*          Lucia Leardini (lucia.leardini@cern.ch)*
* Version 1.0                  *
*                    *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims    *
* about the suitability of this software for any purpose. It is    *
* provided "as is" without express or implied warranty.      *
**************************************************************************/

#if !defined( __CINT__) || defined(__MAKECINT__)

#include <exception>
#include <iostream>
#include "AliLog.h"
#include "AliEventplane.h"
#include "AliConversionAODBGHandlerRP.h"
using namespace std;
#endif

ClassImp(AliConversionAODBGHandlerRP);

//________________________________________________________________________
AliConversionAODBGHandlerRP::AliConversionAODBGHandlerRP(Bool_t IsHeavyIon,Bool_t UseChargedTrackMult,Int_t NEvents) : TObject(),
  fIsHeavyIon(IsHeavyIon),
  fUseChargedTrackMult(UseChargedTrackMult),
  fNEvents(NEvents),
  fBGEventCounter(NULL),
  fNBGEvents(NULL),
  fNBinsRP(8),
  fNBinsZ(7),
  fNBinsMultiplicity(5+Int_t(fUseChargedTrackMult)),
  fBinLimitsArrayRP(NULL),
  fBinLimitsArrayZ(NULL),
  fBinLimitsArrayMultiplicity(NULL),
  fBGEvents(fNBinsRP,AliGammaConversionVertexPositionVector(fNBinsZ,AliGammaConversionBGEventVector(fNEvents)))
//   fBGPool(fNBinsZ,AliGammaConversionMultiplicityVector(fNBinsMultiplicity,AliGammaConversionBGEventVector(fNEvents)))
{
  
  // RP angle Binning  
  fBinLimitsArrayRP = new Double_t[fNBinsRP+1];
  for(Int_t i=0; i < fNBinsRP+1; i++){
    fBinLimitsArrayRP[i] = i*TMath::Pi()/Double_t(fNBinsRP);
  }

  // Vertex Z Binning

  fBinLimitsArrayZ = new Double_t[fNBinsZ+1];
  fBinLimitsArrayZ[0] = -50.00;
  fBinLimitsArrayZ[1] = -3.375;
  fBinLimitsArrayZ[2] = -1.605;
  fBinLimitsArrayZ[3] = -0.225;
  fBinLimitsArrayZ[4] = 1.065;
  fBinLimitsArrayZ[5] = 2.445;
  fBinLimitsArrayZ[6] = 4.245;
  fBinLimitsArrayZ[7] = 50.00;

  // MultiplicityBins
  fBinLimitsArrayMultiplicity= new Double_t[fNBinsMultiplicity+1];

  if(fUseChargedTrackMult){
    // Use Charged Particle Multiplicity
    fBinLimitsArrayMultiplicity[0] = 0;
    fBinLimitsArrayMultiplicity[1] = 8.5;
    fBinLimitsArrayMultiplicity[2] = 16.5;
    fBinLimitsArrayMultiplicity[3] = 27.5;
    fBinLimitsArrayMultiplicity[4] = 41.5;
    fBinLimitsArrayMultiplicity[5] = 200.;

    if(fIsHeavyIon){
        fBinLimitsArrayMultiplicity[0] = 0;
        fBinLimitsArrayMultiplicity[1] = 200.;
        fBinLimitsArrayMultiplicity[2] = 500.;
        fBinLimitsArrayMultiplicity[3] = 1000.;
        fBinLimitsArrayMultiplicity[4] = 1500.;
        fBinLimitsArrayMultiplicity[5] = 5000.;
    }
  } else {
  // Use V0 Multiplicity
    fBinLimitsArrayMultiplicity[0] = 2;
    fBinLimitsArrayMultiplicity[1] = 3;
    fBinLimitsArrayMultiplicity[2] = 4;
    fBinLimitsArrayMultiplicity[3] = 5;
    fBinLimitsArrayMultiplicity[4] = 9999;

    if(fIsHeavyIon){
        fBinLimitsArrayMultiplicity[0] = 2;
        fBinLimitsArrayMultiplicity[1] = 10;
        fBinLimitsArrayMultiplicity[2] = 30;
        fBinLimitsArrayMultiplicity[3] = 50;
        fBinLimitsArrayMultiplicity[4] = 9999;
    }
  }
  Initialize();
}

//________________________________________________________________________
AliConversionAODBGHandlerRP::~AliConversionAODBGHandlerRP()
{
  if(fBinLimitsArrayRP){
    delete[] fBinLimitsArrayRP;
    fBinLimitsArrayRP=0x0;
  }

  if(fBinLimitsArrayZ){
    delete[] fBinLimitsArrayZ;
    fBinLimitsArrayZ=0x0;
  }

  if(fBinLimitsArrayMultiplicity){
    delete[] fBinLimitsArrayMultiplicity;
    fBinLimitsArrayMultiplicity=0x0;
  }

  if(fBGEventCounter){
    for(Int_t psi = 0; psi < fNBinsRP; psi++){
      delete[] fBGEventCounter[psi];
    }
    delete[] fBGEventCounter;
    fBGEventCounter = NULL;
  }

  // Delete pool

  for(Int_t psi = 0; psi < fNBinsRP; psi++){
    for(Int_t z = 0; z < fNBinsZ; z++){
      for(Int_t eventCounter=0; eventCounter < fNBGEvents[psi][z] && eventCounter<fNEvents; eventCounter++){

          for(UInt_t d=0; d < fBGEvents[psi][z][eventCounter].size(); d++){
            delete (AliAODConversionPhoton*)(fBGEvents[psi][z][eventCounter][d]);

          }
      }   
    }
  }

  if(fNBGEvents){
    for(Int_t psi = 0; psi < fNBinsRP; psi++){
      delete[] fNBGEvents[psi];
    }
    delete[] fNBGEvents;
    fNBGEvents = NULL;
  }

}

//________________________________________________________________________
void AliConversionAODBGHandlerRP::Initialize(){

  // Counter

  if(fBGEventCounter == NULL){
    fBGEventCounter = new Int_t*[fNBinsRP];  
  }
  for(Int_t psi = 0; psi < fNBinsRP; psi++){
    fBGEventCounter[psi] = new Int_t[fNBinsZ];
  }
  for(Int_t psi = 0; psi < fNBinsRP; psi++){
    for(Int_t z = 0; z < fNBinsZ; z++){
      fBGEventCounter[psi][z] = 0;
    }
  }

  
  if(fNBGEvents == NULL){
    fNBGEvents = new Int_t*[fNBinsRP];
  }
  for(Int_t psi = 0; psi < fNBinsRP; psi++){
    fNBGEvents[psi] = new Int_t[fNBinsZ];
  }
  for(Int_t psi = 0; psi < fNBinsRP; psi++){
    for(Int_t z = 0; z < fNBinsZ; z++){
      fNBGEvents[psi][z] = 0;
    }
  }
}

//-------------------------------------------------------------
Int_t AliConversionAODBGHandlerRP::GetRPBinIndex(Double_t psivalue) const{

  if(fNBinsRP < 2){
    return 0;
  }

  if(psivalue<=fBinLimitsArrayRP[0]){
    return -1;
  }
  
  for(Int_t i = 0; i < fNBinsRP; i++){
    if(psivalue >= fBinLimitsArrayRP[i] && psivalue <= fBinLimitsArrayRP[i+1]){
      return i;
    }
  }
  return -1;
}

//-------------------------------------------------------------
Int_t AliConversionAODBGHandlerRP::GetZBinIndex(Double_t zvalue) const{

  if(fNBinsZ < 2){
    return 0;
  }

  if(zvalue<=fBinLimitsArrayZ[0]){
    return -1;
  }

  for(Int_t i=0; i < fNBinsZ; i++){
    if(zvalue >= fBinLimitsArrayZ[i] && zvalue <= fBinLimitsArrayZ[i+1]){
      return i;
    }
  }
  return -1;
}

//-------------------------------------------------------------
Int_t AliConversionAODBGHandlerRP::GetMultiplicityBinIndex(Int_t multiplicity) const{
  
  if(fNBinsMultiplicity < 2){
    return 0;
  }

  for(Int_t i=0; i < fNBinsMultiplicity; i++){
    if(multiplicity >= fBinLimitsArrayMultiplicity[i] && multiplicity < fBinLimitsArrayMultiplicity[i+1]){
      return i;
    }
  }
  return -1;
}

//-------------------------------------------------------------
Bool_t AliConversionAODBGHandlerRP::FindBins(TObjArray * const eventGammas,AliVEvent *fInputEvent,Int_t &psibin,Int_t &zbin){
  
  Double_t eventplaneangle;
  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon ==1)eventplaneangle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else eventplaneangle = 0.0;
  psibin = GetRPBinIndex(eventplaneangle);
  
  Double_t vertexz = fInputEvent->GetPrimaryVertex()->GetZ();
  zbin = GetZBinIndex(vertexz);

  if(psibin < fNBinsRP && zbin < fNBinsZ){
    if(psibin >= 0 && zbin >= 0 ){
      return kTRUE;
    }
  }
  //cout<<Form("Requested BG pool does not exist:  z %i m %i",zbin)<<endl;
  return kFALSE;
}

//-------------------------------------------------------------
Bool_t AliConversionAODBGHandlerRP::FindBins(TList * const eventGammas,AliVEvent *fInputEvent,Int_t &psibin,Int_t &zbin){
  
  Double_t eventplaneangle;
  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon ==1)eventplaneangle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else eventplaneangle = 0.0;
  psibin = GetRPBinIndex(eventplaneangle);
  
  Double_t vertexz=fInputEvent->GetPrimaryVertex()->GetZ();
  zbin = GetZBinIndex(vertexz);

  if(psibin < fNBinsRP && zbin < fNBinsZ){
    if(psibin >= 0 && zbin >= 0 ){
      return kTRUE;
    }
  }
  //cout<<Form("Requested BG pool does not exist:  z %i m %i",zbin)<<endl;
  return kFALSE;
}


//-------------------------------------------------------------
void AliConversionAODBGHandlerRP::AddEvent(TObjArray * const eventGammas,AliVEvent *fInputEvent){

  if(eventGammas->GetEntriesFast()==0)return;

  Int_t psi;
  Int_t z;

  if(FindBins(eventGammas,fInputEvent,psi,z)){
    // If Event Stack is full, replace the first entry (First in first out)
    if(fBGEventCounter[psi][z] >= fNEvents){
      fBGEventCounter[psi][z] = 0;
    }

    // Update number of Events stored
    if(fNBGEvents[psi][z] < fNEvents){
      fNBGEvents[psi][z]++;
    }

    Int_t eventCounter = fBGEventCounter[psi][z];

    //clear the vector for old gammas
    for(UInt_t d = 0; d < fBGEvents[psi][z][eventCounter].size(); d++){
      delete (AliAODConversionPhoton*)(fBGEvents[psi][z][eventCounter][d]);
    }

    fBGEvents[psi][z][eventCounter].clear();

    // add the gammas to the vector
    for(Int_t i = 0; i < eventGammas->GetEntriesFast(); i++){
      fBGEvents[psi][z][eventCounter].push_back(new AliAODConversionPhoton(*(AliAODConversionPhoton*)(eventGammas->At(i))));
    }

    fBGEventCounter[psi][z]++;
  }
}
//-------------------------------------------------------------
void AliConversionAODBGHandlerRP::AddEvent(TList * const eventGammas,AliVEvent *fInputEvent){
  if(eventGammas->GetEntries()==0)return;

  Int_t psi;
  Int_t z;

  if(FindBins(eventGammas,fInputEvent,psi,z)){
    // If Event Stack is full, replace the first entry (First in first out)
    if(fBGEventCounter[psi][z] >= fNEvents){
        fBGEventCounter[psi][z]=0;
    }

    // Update number of Events stored
    if(fNBGEvents[psi][z] < fNEvents){
      fNBGEvents[psi][z]++;
    }

    Int_t eventCounter = fBGEventCounter[psi][z];

    //clear the vector for old gammas
    for(UInt_t d = 0; d < fBGEvents[psi][z][eventCounter].size(); d++){
      delete (AliAODConversionPhoton*)(fBGEvents[psi][z][eventCounter][d]);
    }

    fBGEvents[psi][z][eventCounter].clear();

    // add the gammas to the vector
    for(Int_t i = 0; i < eventGammas->GetEntries(); i++){
      fBGEvents[psi][z][eventCounter].push_back(new AliAODConversionPhoton(*(AliAODConversionPhoton*)(eventGammas->At(i))));
    }

    fBGEventCounter[psi][z]++;
  }
}

//-------------------------------------------------------------
AliGammaConversionPhotonVector* AliConversionAODBGHandlerRP::GetBGGoodGammas(TObjArray * const eventGammas,AliVEvent *fInputEvent,Int_t event){
  Int_t psibin;
  Int_t zbin;

  if(FindBins(eventGammas,fInputEvent,psibin,zbin)){
    return &(fBGEvents[psibin][zbin][event]);
  }
  return NULL;
}
//-------------------------------------------------------------
AliGammaConversionPhotonVector* AliConversionAODBGHandlerRP::GetBGGoodGammas(TList * const eventGammas,AliVEvent *fInputEvent,Int_t event){
  Int_t psibin;
  Int_t zbin;

  if(FindBins(eventGammas,fInputEvent,psibin,zbin)){
    return &(fBGEvents[psibin][zbin][event]);
  }
  return NULL;
}

//-------------------------------------------------------------
Int_t AliConversionAODBGHandlerRP::GetNBGEvents(TObjArray * const eventGammas,AliVEvent *fInputEvent){
  Int_t psibin;
  Int_t zbin;

  if(FindBins(eventGammas,fInputEvent,psibin,zbin)){
    return fNBGEvents[psibin][zbin];
  }
  return 0;
}
//-------------------------------------------------------------
Int_t AliConversionAODBGHandlerRP::GetNBGEvents(TList * const eventGammas,AliVEvent *fInputEvent){
  Int_t psibin;
  Int_t zbin;

  if(FindBins(eventGammas,fInputEvent,psibin,zbin)){
    return fNBGEvents[psibin][zbin];
  }
  return 0;
}

