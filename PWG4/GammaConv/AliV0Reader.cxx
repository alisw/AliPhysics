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
// Class used to do analysis on conversion pairs
//---------------------------------------------
////////////////////////////////////////////////

// --- ROOT system ---
#include <TMath.h>

//---- ANALYSIS system ----
#include "AliV0Reader.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliKFVertex.h"

#include "AliStack.h"
#include "AliMCEventHandler.h"


class iostream;
class AliESDv0;
class TFormula;

using namespace std;

ClassImp(AliV0Reader)



  AliV0Reader::AliV0Reader() :
    TObject(),
    fMCStack(NULL),
    fMCTruth(NULL),
    fChain(NULL),
    fESDHandler(NULL),
    fESDEvent(NULL),
    fHistograms(NULL),
    fCurrentV0IndexNumber(0),
    fCurrentV0(NULL),
    fCurrentNegativeKFParticle(NULL),
    fCurrentPositiveKFParticle(NULL),
    fCurrentMotherKFCandidate(NULL),
    fCurrentNegativeESDTrack(NULL),
    fCurrentPositiveESDTrack(NULL),
    fNegativeTrackLorentzVector(NULL),
    fPositiveTrackLorentzVector(NULL),
    fMotherCandidateLorentzVector(NULL),
    fCurrentXValue(0),
    fCurrentYValue(0),
    fCurrentZValue(0),
    fPositiveTrackPID(0),
    fNegativeTrackPID(0),
    fNegativeMCParticle(NULL),
    fPositiveMCParticle(NULL),
    fMotherMCParticle(NULL),
    fMotherCandidateKFMass(0),
    fMotherCandidateKFWidth(0),
    fUseKFParticle(kTRUE),
    fUseESDTrack(kFALSE),
    fDoMC(kFALSE),
    fMaxR(10000),// 100 meter(outside of ALICE)
    fEtaCut(0.),
    fPtCut(0.),
    fChi2CutConversion(0.),
    fChi2CutMeson(0.),
    fPIDProbabilityCutNegativeParticle(0),
    fPIDProbabilityCutPositiveParticle(0),
    fXVertexCut(0.),
    fYVertexCut(0.),
    fZVertexCut(0.),
    fNSigmaMass(0.),
    fUseImprovedVertex(kFALSE),
    fCurrentEventGoodV0s(),
    fPreviousEventGoodV0s()
{

}


AliV0Reader::AliV0Reader(const AliV0Reader & original) :
  TObject(original),
  fMCStack(original.fMCStack),
  fMCTruth(original.fMCTruth),
  fChain(original.fChain),
  fESDHandler(original.fESDHandler),
  fESDEvent(original.fESDEvent),
  fHistograms(original.fHistograms),
  fCurrentV0IndexNumber(original.fCurrentV0IndexNumber),
  fCurrentV0(original.fCurrentV0),
  fCurrentNegativeKFParticle(original.fCurrentNegativeKFParticle),
  fCurrentPositiveKFParticle(original.fCurrentPositiveKFParticle),
  fCurrentMotherKFCandidate(original.fCurrentMotherKFCandidate),
  fCurrentNegativeESDTrack(original.fCurrentNegativeESDTrack),
  fCurrentPositiveESDTrack(original.fCurrentPositiveESDTrack),
  fNegativeTrackLorentzVector(original.fNegativeTrackLorentzVector),
  fPositiveTrackLorentzVector(original.fPositiveTrackLorentzVector),
  fMotherCandidateLorentzVector(original.fMotherCandidateLorentzVector),
  fCurrentXValue(original.fCurrentXValue),
  fCurrentYValue(original.fCurrentYValue),
  fCurrentZValue(original.fCurrentZValue),
  fPositiveTrackPID(original.fPositiveTrackPID),
  fNegativeTrackPID(original.fNegativeTrackPID),
  fNegativeMCParticle(original.fNegativeMCParticle),
  fPositiveMCParticle(original.fPositiveMCParticle),
  fMotherMCParticle(original.fMotherMCParticle),
  fMotherCandidateKFMass(original.fMotherCandidateKFMass),
  fMotherCandidateKFWidth(original.fMotherCandidateKFWidth),
  fUseKFParticle(kTRUE),
  fUseESDTrack(kFALSE),
  fDoMC(kFALSE),
  fMaxR(original.fMaxR),
  fEtaCut(original.fEtaCut),
  fPtCut(original.fPtCut),
  fChi2CutConversion(original.fChi2CutConversion),
  fChi2CutMeson(original.fChi2CutMeson),
  fPIDProbabilityCutNegativeParticle(original.fPIDProbabilityCutNegativeParticle),
  fPIDProbabilityCutPositiveParticle(original.fPIDProbabilityCutPositiveParticle),
  fXVertexCut(original.fXVertexCut),
  fYVertexCut(original.fYVertexCut),
  fZVertexCut(original.fZVertexCut),
  fNSigmaMass(original.fNSigmaMass),
  fUseImprovedVertex(original.fUseImprovedVertex),
  fCurrentEventGoodV0s(original.fCurrentEventGoodV0s),
  fPreviousEventGoodV0s(original.fPreviousEventGoodV0s)
{

}


AliV0Reader & AliV0Reader::operator = (const AliV0Reader & /*source*/)
{
  // assignment operator
  return *this;
}

void AliV0Reader::Initialize(){
  //see header file for documentation

  // Get the input handler from the manager
  fESDHandler = (AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(fESDHandler == NULL){
    //print warning here
  }
  
  // Get pointer to esd event from input handler
  fESDEvent = fESDHandler->GetEvent();
  if(fESDEvent == NULL){
    //print warning here
  }

  //Get pointer to MCTruth
  fMCTruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
  if(fMCTruth == NULL){
    //print warning here
  }

  //Get pointer to the mc stack
  fMCStack = fMCTruth->MCEvent()->Stack();
  if(fMCStack == NULL){
    //print warning here
  }

  AliKFParticle::SetField(fESDEvent->GetMagneticField());

}

AliESDv0* AliV0Reader::GetV0(Int_t index){
  //see header file for documentation

  fCurrentV0 = fESDEvent->GetV0(index);
  UpdateV0Information();
  return fCurrentV0;
}
Bool_t AliV0Reader::CheckForPrimaryVertex(){
  return fESDEvent->GetPrimaryVertex()->GetNContributors()>0;
}

Bool_t AliV0Reader::NextV0(){
  //see header file for documentation

  Bool_t iResult=kFALSE;
  while(fCurrentV0IndexNumber<fESDEvent->GetNumberOfV0s()){
    fCurrentV0 = fESDEvent->GetV0(fCurrentV0IndexNumber);
    
    //checks if on the fly mode is set
    if ( !fCurrentV0->GetOnFlyStatus() ){
      fCurrentV0IndexNumber++;
      if(fHistograms != NULL){
	fHistograms->FillHistogram("V0MassDebugCut1",GetMotherCandidateMass());
      }
      continue;
    }

    if(fESDEvent->GetPrimaryVertex()->GetNContributors()<=0) {//checks if we have a vertex
      fCurrentV0IndexNumber++;
      if(fHistograms != NULL){
	fHistograms->FillHistogram("V0MassDebugCut2",GetMotherCandidateMass());
      }
      continue;
    }

    if(CheckPIDProbability(fPIDProbabilityCutNegativeParticle,fPIDProbabilityCutPositiveParticle)==kFALSE){
      fCurrentV0IndexNumber++;
      if(fHistograms != NULL){
	fHistograms->FillHistogram("V0MassDebugCut3",GetMotherCandidateMass());
      }
      continue;
    }

    fCurrentV0->GetXYZ(fCurrentXValue,fCurrentYValue,fCurrentZValue);
 
    if(GetXYRadius()>fMaxR){ // cuts on distance from collision point
      fCurrentV0IndexNumber++;
      if(fHistograms != NULL){
	fHistograms->FillHistogram("V0MassDebugCut4",GetMotherCandidateMass());
      }
      continue;
    }

    UpdateV0Information();
        
    if(fUseKFParticle){
      if(fCurrentMotherKFCandidate->GetNDF()<=0){
	fCurrentV0IndexNumber++;
	if(fHistograms != NULL){
	  fHistograms->FillHistogram("V0MassDebugCut5",GetMotherCandidateMass());
	}
	continue;
      }
      Double_t chi2V0 = fCurrentMotherKFCandidate->GetChi2()/fCurrentMotherKFCandidate->GetNDF();
      if(chi2V0 > fChi2CutConversion || chi2V0 <=0){
	fCurrentV0IndexNumber++;
	if(fHistograms != NULL){
	  fHistograms->FillHistogram("V0MassDebugCut6",GetMotherCandidateMass());
	}
  	continue;
      }
      
      if(TMath::Abs(fMotherCandidateLorentzVector->Eta())> fEtaCut){
	fCurrentV0IndexNumber++;
	if(fHistograms != NULL){
	  fHistograms->FillHistogram("V0MassDebugCut7",GetMotherCandidateMass());
	}
	continue;
      }
      
      if(fMotherCandidateLorentzVector->Pt()<fPtCut){
	fCurrentV0IndexNumber++;
	if(fHistograms != NULL){
	  fHistograms->FillHistogram("V0MassDebugCut8",GetMotherCandidateMass());
	}
	continue;
      }
    }
    else if(fUseESDTrack){
      //TODO
    }

    iResult=kTRUE;//means we have a v0 who survived all the cuts applied

    fCurrentV0IndexNumber++;
    
    break;
  }
  return iResult; 
}

void AliV0Reader::UpdateV0Information(){
  //see header file for documentation
  
  if(fCurrentNegativeKFParticle != NULL){
    delete fCurrentNegativeKFParticle;
  }
  fCurrentNegativeKFParticle = new AliKFParticle(*(fCurrentV0->GetParamN()),fNegativeTrackPID);
  
  if(fCurrentPositiveKFParticle != NULL){
    delete fCurrentPositiveKFParticle;
  }
  fCurrentPositiveKFParticle = new AliKFParticle(*(fCurrentV0->GetParamP()),fPositiveTrackPID);
    
  if(fCurrentMotherKFCandidate != NULL){
    delete fCurrentMotherKFCandidate;
  }
  fCurrentMotherKFCandidate = new AliKFParticle(*fCurrentNegativeKFParticle,*fCurrentPositiveKFParticle);

  fCurrentNegativeESDTrack = fESDEvent->GetTrack(fCurrentV0->GetNindex());

  fCurrentPositiveESDTrack = fESDEvent->GetTrack(fCurrentV0->GetPindex());

  if(fPositiveTrackPID==-11 && fNegativeTrackPID==11){
    fCurrentMotherKFCandidate->SetMassConstraint(0,fNSigmaMass);
  }

  if(fUseImprovedVertex == kTRUE){
    AliKFVertex primaryVertexImproved(*GetPrimaryVertex());
    primaryVertexImproved+=*fCurrentMotherKFCandidate;
    fCurrentMotherKFCandidate->SetProductionVertex(primaryVertexImproved);
  }

  fCurrentMotherKFCandidate->GetMass(fMotherCandidateKFMass,fMotherCandidateKFWidth);


  if(fNegativeTrackLorentzVector != NULL){
    delete fNegativeTrackLorentzVector;
  }
  if(fUseKFParticle){
    fNegativeTrackLorentzVector = new TLorentzVector(fCurrentNegativeKFParticle->Px(),fCurrentNegativeKFParticle->Py(),fCurrentNegativeKFParticle->Pz());
  }
  else if(fUseESDTrack){
    fNegativeTrackLorentzVector = new TLorentzVector(fCurrentNegativeESDTrack->Px(),fCurrentNegativeESDTrack->Py(),fCurrentNegativeESDTrack->Pz());
  }

  if(fPositiveTrackLorentzVector != NULL){
    delete fPositiveTrackLorentzVector;
  }
  if(fUseKFParticle){
    fPositiveTrackLorentzVector = new TLorentzVector(fCurrentPositiveKFParticle->Px(),fCurrentPositiveKFParticle->Py(),fCurrentPositiveKFParticle->Pz());
  }
  else if(fUseESDTrack){
    fPositiveTrackLorentzVector = new TLorentzVector(fCurrentPositiveESDTrack->Px(),fCurrentPositiveESDTrack->Py(),fCurrentPositiveESDTrack->Pz());
  }

  if(fMotherCandidateLorentzVector != NULL){
    delete fMotherCandidateLorentzVector;
  }
  if(fUseKFParticle){
    fMotherCandidateLorentzVector = new TLorentzVector(*fNegativeTrackLorentzVector + *fPositiveTrackLorentzVector);
  }
  else if(fUseESDTrack){
    fMotherCandidateLorentzVector = new TLorentzVector(*fNegativeTrackLorentzVector + *fPositiveTrackLorentzVector);
  }

  if(fPositiveTrackPID==-11 && fNegativeTrackPID==11){
    fMotherCandidateLorentzVector->SetXYZM(fMotherCandidateLorentzVector->Px() ,fMotherCandidateLorentzVector->Py(),fMotherCandidateLorentzVector->Pz(),0.); 
  }
    
  if(fDoMC == kTRUE){
    fMotherMCParticle= NULL;
    fNegativeMCParticle = fMCStack->Particle(TMath::Abs(fESDEvent->GetTrack(fCurrentV0->GetNindex())->GetLabel()));
    fPositiveMCParticle = fMCStack->Particle(TMath::Abs(fESDEvent->GetTrack(fCurrentV0->GetPindex())->GetLabel()));
    if(fPositiveMCParticle->GetMother(0)>-1){
      fMotherMCParticle = fMCStack->Particle(fPositiveMCParticle->GetMother(0));
    }
  }
  fCurrentEventGoodV0s.push_back(*fCurrentMotherKFCandidate);
}

Bool_t AliV0Reader::HasSameMCMother(){
  //see header file for documentation

  Bool_t iResult = kFALSE;
  if(fDoMC == kTRUE){
    if(fNegativeMCParticle != NULL && fPositiveMCParticle != NULL){
      if(fNegativeMCParticle->GetMother(0) == fPositiveMCParticle->GetMother(0))
	if(fMotherMCParticle){
	  iResult = kTRUE;
	}
    }
  }
  return iResult;
}

Bool_t AliV0Reader::CheckPIDProbability(Double_t negProbCut, Double_t posProbCut){
  //see header file for documentation

  Bool_t iResult=kFALSE;

  Double_t *posProbArray = new Double_t[10];
  Double_t *negProbArray = new Double_t[10];
  AliESDtrack* negTrack  = fESDEvent->GetTrack(fCurrentV0->GetNindex());
  AliESDtrack* posTrack  = fESDEvent->GetTrack(fCurrentV0->GetPindex());
  
  negTrack->GetTPCpid(negProbArray);
  posTrack->GetTPCpid(posProbArray);

  if(negProbArray!=NULL && posProbArray!=NULL){
    if(negProbArray[GetSpeciesIndex(-1)]>=negProbCut && posProbArray[GetSpeciesIndex(1)]>=posProbCut){
      iResult=kTRUE;
    }
  }
  delete [] posProbArray;
  delete [] negProbArray;
  return iResult;
}

void AliV0Reader::GetPIDProbability(Double_t &negPIDProb,Double_t & posPIDProb){

  Double_t *posProbArray = new Double_t[10];
  Double_t *negProbArray = new Double_t[10];
  AliESDtrack* negTrack  = fESDEvent->GetTrack(fCurrentV0->GetNindex());
  AliESDtrack* posTrack  = fESDEvent->GetTrack(fCurrentV0->GetPindex());
  
  negTrack->GetTPCpid(negProbArray);
  posTrack->GetTPCpid(posProbArray);

  if(negProbArray!=NULL && posProbArray!=NULL){
    negPIDProb = negProbArray[GetSpeciesIndex(-1)];
    posPIDProb = posProbArray[GetSpeciesIndex(1)];
  }
  delete [] posProbArray;
  delete [] negProbArray;
}

void AliV0Reader::UpdateEventByEventData(){
  //see header file for documentation

  if(fCurrentEventGoodV0s.size() >0 ){
    fPreviousEventGoodV0s.clear();
    fPreviousEventGoodV0s = fCurrentEventGoodV0s;
  }
  fCurrentEventGoodV0s.clear();
  
  fCurrentV0IndexNumber=0;
}

Double_t AliV0Reader::GetNegativeTrackPhi() const{
  //see header file for documentation

  Double_t offset=0;
  if(fNegativeTrackLorentzVector->Phi()> TMath::Pi()){
    offset = -2*TMath::Pi();
  }
  return fNegativeTrackLorentzVector->Phi()+offset;
}

Double_t AliV0Reader::GetPositiveTrackPhi() const{
  //see header file for documentation

  Double_t offset=0;
  if(fPositiveTrackLorentzVector->Phi()> TMath::Pi()){
    offset = -2*TMath::Pi();
  }
  return fPositiveTrackLorentzVector->Phi()+offset;
}

Double_t AliV0Reader::GetMotherCandidatePhi() const{
  //see header file for documentation

  Double_t offset=0;
  if(fMotherCandidateLorentzVector->Phi()> TMath::Pi()){
    offset = -2*TMath::Pi();
  }
  return fMotherCandidateLorentzVector->Phi()+offset;
}

Int_t AliV0Reader::GetSpeciesIndex(Int_t chargeOfTrack){
  //see header file for documentation

  Int_t iResult = 10; // Unknown particle

  if(chargeOfTrack==-1){ //negative track
    switch(abs(fNegativeTrackPID)){
    case 11:       //electron
      iResult = 0;
      break;
    case 13:       //muon
      iResult = 1;
      break;
    case 211:      //pion
      iResult = 2;
      break;
    case 321:      //kaon
      iResult = 3;
      break;
    case 2212:     //proton
      iResult = 4;
      break;
    case 22:       //photon
      iResult = 5;
      break;
    case 111:      //pi0
      iResult = 6;
      break;
    case 2112:     //neutron
      iResult = 7;
      break;
    case 311:      //K0
      iResult = 8;
      break;
      
      //Put in here for kSPECIES::kEleCon  ????
    }
  }
  else if(chargeOfTrack==1){ //positive track
    switch(abs(fPositiveTrackPID)){
    case 11:       //electron
      iResult = 0;
      break;
    case 13:       //muon
      iResult = 1;
      break;
    case 211:      //pion
      iResult = 2;
      break;
    case 321:      //kaon
      iResult = 3;
      break;
    case 2212:     //proton
      iResult = 4;
      break;
    case 22:       //photon
      iResult = 5;
      break;
    case 111:      //pi0
      iResult = 6;
      break;
    case 2112:     //neutron
      iResult = 7;
      break;
    case 311:      //K0
      iResult = 8;
      break;

      //Put in here for kSPECIES::kEleCon  ????
    }
  }
  else{
    //Wrong parameter.. Print warning
  }
  return iResult;
}
