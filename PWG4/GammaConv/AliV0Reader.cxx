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
#include "AliESDtrack.h"
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
  fMCEvent(NULL),    // for CF
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
  fMaxZ(0.),
  fLineCutZRSlope(0.),
  fLineCutZValue(0.),
  fChi2CutConversion(0.),
  fChi2CutMeson(0.),
  fPIDProbabilityCutNegativeParticle(0),
  fPIDProbabilityCutPositiveParticle(0),
  fXVertexCut(0.),
  fYVertexCut(0.),
  fZVertexCut(0.),
  fNSigmaMass(0.),
  fUseImprovedVertex(kFALSE),
  fUseOwnXYZCalculation(kFALSE),
  fCurrentEventGoodV0s(),
  fPreviousEventGoodV0s()
{
	
}


AliV0Reader::AliV0Reader(const AliV0Reader & original) :
  TObject(original),
  fMCStack(original.fMCStack),
  fMCTruth(original.fMCTruth),
  fMCEvent(original.fMCEvent),  // for CF
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
  fMaxZ(original.fMaxZ),
  fLineCutZRSlope(original.fLineCutZRSlope),
  fLineCutZValue(original.fLineCutZValue),
  fChi2CutConversion(original.fChi2CutConversion),
  fChi2CutMeson(original.fChi2CutMeson),
  fPIDProbabilityCutNegativeParticle(original.fPIDProbabilityCutNegativeParticle),
  fPIDProbabilityCutPositiveParticle(original.fPIDProbabilityCutPositiveParticle),
  fXVertexCut(original.fXVertexCut),
  fYVertexCut(original.fYVertexCut),
  fZVertexCut(original.fZVertexCut),
  fNSigmaMass(original.fNSigmaMass),
  fUseImprovedVertex(original.fUseImprovedVertex),
  fUseOwnXYZCalculation(original.fUseOwnXYZCalculation),
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
	
	
  // for CF
  //Get pointer to the mc event
  fMCEvent = fMCTruth->MCEvent();
  if(fMCEvent == NULL){
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
		
    // moved it up here so that the correction framework can access pt and eta information
    if(UpdateV0Information() == kFALSE){
      fCurrentV0IndexNumber++;
      continue;
    }

    Double_t containerInput[3];
    containerInput[0] = GetMotherCandidatePt();
    containerInput[1] = GetMotherCandidateEta();
    containerInput[2] = GetMotherCandidateMass();
		

    //checks if on the fly mode is set
    if ( !fCurrentV0->GetOnFlyStatus() ){
      if(fHistograms != NULL){
	fHistograms->FillHistogram("ESD_CutGetOnFly_InvMass",GetMotherCandidateMass());
      }
      fCurrentV0IndexNumber++;
      continue;
    }
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepGetOnFly);		// for CF	
		
    //checks if we have a prim vertex
    if(fESDEvent->GetPrimaryVertex()->GetNContributors()<=0) { 
      if(fHistograms != NULL){
	fHistograms->FillHistogram("ESD_CutNContributors_InvMass",GetMotherCandidateMass());
      }
      fCurrentV0IndexNumber++;
      continue;
    }
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepNContributors);		// for CF	
		
		
    //Check the pid probability
    if(CheckPIDProbability(fPIDProbabilityCutNegativeParticle,fPIDProbabilityCutPositiveParticle)==kFALSE){
      if(fHistograms != NULL){
	fHistograms->FillHistogram("ESD_CutPIDProb_InvMass",GetMotherCandidateMass());
      }
      fCurrentV0IndexNumber++;
      continue;
    }
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepTPCPID);			// for CF
		
		
		
		
    fCurrentV0->GetXYZ(fCurrentXValue,fCurrentYValue,fCurrentZValue);
		
		
    if(GetXYRadius()>fMaxR){ // cuts on distance from collision point
      if(fHistograms != NULL){
	fHistograms->FillHistogram("ESD_CutR_InvMass",GetMotherCandidateMass());
      }
      fCurrentV0IndexNumber++;
      continue;
    }		
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepR);			// for CF
		
		
		
    if((TMath::Abs(fCurrentZValue)*fLineCutZRSlope)-fLineCutZValue > GetXYRadius() ){ // cuts out regions where we do not reconstruct
      if(fHistograms != NULL){
	fHistograms->FillHistogram("ESD_CutLine_InvMass",GetMotherCandidateMass());
      }
      fCurrentV0IndexNumber++;
      continue;
    }		
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepLine);			// for CF
		
		
    if(TMath::Abs(fCurrentZValue) > fMaxZ ){ // cuts out regions where we do not reconstruct
      if(fHistograms != NULL){
	fHistograms->FillHistogram("ESD_CutZ_InvMass",GetMotherCandidateMass());
      }
      fCurrentV0IndexNumber++;
      continue;
    }		
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepZ);		// for CF	
		
		
    /* Moved further up so corr framework can work
       if(UpdateV0Information() == kFALSE){
       fCurrentV0IndexNumber++;
       continue;
       }
    */
		
    if(fUseKFParticle){
      if(fCurrentMotherKFCandidate->GetNDF()<=0){
	if(fHistograms != NULL){
	  fHistograms->FillHistogram("ESD_CutNDF_InvMass",GetMotherCandidateMass());
	}
	fCurrentV0IndexNumber++;
	continue;
      }
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepNDF);		// for CF	
			
			
      Double_t chi2V0 = fCurrentMotherKFCandidate->GetChi2()/fCurrentMotherKFCandidate->GetNDF();
      if(chi2V0 > fChi2CutConversion || chi2V0 <=0){
	if(fHistograms != NULL){
	  fHistograms->FillHistogram("ESD_CutChi2_InvMass",GetMotherCandidateMass());
	}
	fCurrentV0IndexNumber++;
	continue;
      }
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepChi2);			// for CF
			
			
      if(TMath::Abs(fMotherCandidateLorentzVector->Eta())> fEtaCut){
	if(fHistograms != NULL){
	  fHistograms->FillHistogram("ESD_CutEta_InvMass",GetMotherCandidateMass());
	}
	fCurrentV0IndexNumber++;
	continue;
      }
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepEta);			// for CF
			
			
      if(fMotherCandidateLorentzVector->Pt()<fPtCut){
	if(fHistograms != NULL){
	  fHistograms->FillHistogram("ESD_CutPt_InvMass",GetMotherCandidateMass());
	}
	fCurrentV0IndexNumber++;
	continue;
      }
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepPt);			// for CF
			
			
    }
    else if(fUseESDTrack){
      //TODO
    }
		
    fCurrentEventGoodV0s.push_back(*fCurrentMotherKFCandidate);
		
    iResult=kTRUE;//means we have a v0 who survived all the cuts applied
		
    fCurrentV0IndexNumber++;
		
    break;
  }
  return iResult; 
}

Bool_t AliV0Reader::UpdateV0Information(){
  //see header file for documentation
	
  Bool_t iResult=kTRUE;		 				// for taking out not refitted, kinks and like sign tracks 
	
  Bool_t switchTracks = kFALSE;
	
  fCurrentNegativeESDTrack = fESDEvent->GetTrack(fCurrentV0->GetNindex());
  fCurrentPositiveESDTrack = fESDEvent->GetTrack(fCurrentV0->GetPindex());
	
  if(fCurrentNegativeESDTrack->GetSign() == fCurrentPositiveESDTrack->GetSign()){             // avoid like sign
    iResult=kFALSE;
    if(fHistograms != NULL){
      fHistograms->FillHistogram("ESD_CutLikeSign_InvMass",GetMotherCandidateMass());
    }
  }
	
  if(fCurrentPositiveESDTrack->GetSign() == -1 && fCurrentNegativeESDTrack->GetSign() == 1){  // switch wrong signed tracks
    fCurrentNegativeESDTrack = fESDEvent->GetTrack(fCurrentV0->GetPindex());
    fCurrentPositiveESDTrack = fESDEvent->GetTrack(fCurrentV0->GetNindex());
    switchTracks = kTRUE;
  }
	
  if( !(fCurrentNegativeESDTrack->GetStatus() & AliESDtrack::kTPCrefit) || 
      !(fCurrentPositiveESDTrack->GetStatus() & AliESDtrack::kTPCrefit) ){
    //  if( !(fCurrentNegativeESDTrack->GetStatus() & AliESDtrack::kITSrefit) || 
    //      !(fCurrentPositiveESDTrack->GetStatus() & AliESDtrack::kITSrefit) ){
		
    iResult=kFALSE;
    if(fHistograms != NULL){
      fHistograms->FillHistogram("ESD_CutRefit_InvMass",GetMotherCandidateMass());
    }
  }
	
  if( fCurrentNegativeESDTrack->GetKinkIndex(0) > 0 || 
      fCurrentPositiveESDTrack->GetKinkIndex(0) > 0) {			
		
    iResult=kFALSE;
    if(fHistograms != NULL){
      fHistograms->FillHistogram("ESD_CutKink_InvMass",GetMotherCandidateMass());
    }
  }
	
  if(fCurrentNegativeKFParticle != NULL){
    delete fCurrentNegativeKFParticle;
  }
  if(switchTracks == kFALSE){
    fCurrentNegativeKFParticle = new AliKFParticle(*(fCurrentV0->GetParamN()),fNegativeTrackPID);
  }
  else{
    fCurrentNegativeKFParticle = new AliKFParticle(*(fCurrentV0->GetParamP()),fNegativeTrackPID);
  }
	
  if(fCurrentPositiveKFParticle != NULL){
    delete fCurrentPositiveKFParticle;
  }
  if(switchTracks == kFALSE){
    fCurrentPositiveKFParticle = new AliKFParticle(*(fCurrentV0->GetParamP()),fPositiveTrackPID);
  }
  else{
    fCurrentPositiveKFParticle = new AliKFParticle(*(fCurrentV0->GetParamN()),fPositiveTrackPID);
  }
    
  if(fCurrentMotherKFCandidate != NULL){
    delete fCurrentMotherKFCandidate;
  }
  fCurrentMotherKFCandidate = new AliKFParticle(*fCurrentNegativeKFParticle,*fCurrentPositiveKFParticle);
	
	
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
	
  //  if(iResult==kTRUE){
  //	fCurrentEventGoodV0s.push_back(*fCurrentMotherKFCandidate); // moved it to NextV0() after all the cuts are applied
  //  }


  // for CF
  Double_t containerInput[3];
  containerInput[0] = GetMotherCandidatePt();
  containerInput[1] = GetMotherCandidateEta();
  containerInput[2] = GetMotherCandidateMass();

  fCFManager->GetParticleContainer()->Fill(containerInput,kStepLikeSign);		// for CF	
  fCFManager->GetParticleContainer()->Fill(containerInput,kStepTPCRefit);		// for CF	
  fCFManager->GetParticleContainer()->Fill(containerInput,kStepKinks);		// for CF	

  return iResult;
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
    //    fPreviousEventGoodV0s.clear();
    //    fPreviousEventGoodV0s = fCurrentEventGoodV0s;
    if(fPreviousEventGoodV0s.size()>19){
      for(UInt_t nCurrent=0;nCurrent<fCurrentEventGoodV0s.size();nCurrent++){
	fPreviousEventGoodV0s.erase(fPreviousEventGoodV0s.begin());
	fPreviousEventGoodV0s.push_back(fCurrentEventGoodV0s.at(nCurrent));
      }
    }
    else{
      for(UInt_t nCurrent=0;nCurrent<fCurrentEventGoodV0s.size();nCurrent++){
	if(fPreviousEventGoodV0s.size()<20){
	  fPreviousEventGoodV0s.push_back(fCurrentEventGoodV0s.at(nCurrent));
	}
	else{
	  fPreviousEventGoodV0s.erase(fPreviousEventGoodV0s.begin());
	  fPreviousEventGoodV0s.push_back(fCurrentEventGoodV0s.at(nCurrent));
	}
      }
    }
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


Double_t AliV0Reader::GetMotherCandidateRapidity() const{
  //see header file for documentation
	
  Double_t rapidity=0;
  if(fMotherCandidateLorentzVector->Energy() - fMotherCandidateLorentzVector->Pz() == 0 || fMotherCandidateLorentzVector->Energy() + fMotherCandidateLorentzVector->Pz() == 0) rapidity=0;
  else rapidity = 0.5*(TMath::Log((fMotherCandidateLorentzVector->Energy() + fMotherCandidateLorentzVector->Pz()) / (fMotherCandidateLorentzVector->Energy()-fMotherCandidateLorentzVector->Pz())));
  return rapidity;
	
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

Bool_t GetHelixCenter(AliESDtrack* track, Double_t b,Int_t charge, Double_t center[2]){
  // see header file for documentation
  
  Double_t pi = 3.14159265358979323846;
  
  Double_t  helix[6];
  track->GetHelixParameters(helix,b);
  
  Double_t xpos =  helix[5];
  Double_t ypos =  helix[0];
  Double_t radius = TMath::Abs(1./helix[4]);
  Double_t phi = helix[2];

  if(phi < 0){
    phi = phi + 2*pi;
  }

  phi -= pi/2.;
  Double_t xpoint =  radius * TMath::Cos(phi);
  Double_t ypoint =  radius * TMath::Sin(phi);

  if(charge > 0){
    xpoint = - xpoint;
    ypoint = - ypoint;
  }

  if(charge < 0){
    xpoint =  xpoint;
    ypoint =  ypoint;
  }
  center[0] =  xpos + xpoint;
  center[1] =  ypos + ypoint;

  return 1;
}

Bool_t GetConvPosXY(AliESDtrack* ptrack,AliESDtrack* ntrack, Double_t b, Double_t convpos[2]){
  //see header file for documentation

  Double_t helixcenterpos[2];
  GetHelixCenter(ptrack,b,ptrack->Charge(),helixcenterpos);

  Double_t helixcenterneg[2];
  GetHelixCenter(ntrack,b,ntrack->Charge(),helixcenterneg);

  Double_t  poshelix[6];
  ptrack->GetHelixParameters(poshelix,b);
  Double_t posradius = TMath::Abs(1./poshelix[4]);

  Double_t  neghelix[6];
  ntrack->GetHelixParameters(neghelix,b);
  Double_t negradius = TMath::Abs(1./neghelix[4]);

  Double_t xpos = helixcenterpos[0];
  Double_t ypos = helixcenterpos[1];
  Double_t xneg = helixcenterneg[0];
  Double_t yneg = helixcenterneg[1];

  convpos[0] = (xpos*negradius + xneg*posradius)/(negradius+posradius);
  convpos[1] = (ypos*negradius+  yneg*posradius)/(negradius+posradius);

  return 1;
}



Double_t GetConvPosZ(AliESDtrack* ptrack,AliESDtrack* ntrack, Double_t b){
  //see header file for documentation

  Double_t  helixpos[6];
  ptrack->GetHelixParameters(helixpos,b);

  Double_t  helixneg[6];
  ntrack->GetHelixParameters(helixneg,b);

  Double_t negtrackradius =  TMath::Abs(1./helixneg[4]);
  Double_t postrackradius =  TMath::Abs(1./helixpos[4]);

  Double_t pi = 3.14159265358979323846;

  Double_t convpos[2];
  GetConvPosXY(ptrack,ntrack,b,convpos);

   Double_t convposx = convpos[0];
   Double_t convposy = convpos[1];

   Double_t helixcenterpos[2];
   GetHelixCenter(ptrack,b,ptrack->Charge(),helixcenterpos);

   Double_t helixcenterneg[2];
   GetHelixCenter(ntrack,b,ntrack->Charge(),helixcenterneg);

   Double_t xpos = helixcenterpos[0];
   Double_t ypos = helixcenterpos[1];
   Double_t xneg = helixcenterneg[0];
   Double_t yneg = helixcenterneg[1];

   Double_t delta_x_pos = convposx -  xpos;
   Double_t delta_y_pos = convposy -  ypos;

   Double_t delta_x_neg = convposx -  xneg;
   Double_t delta_y_neg = convposy -  yneg;

   Double_t alpha_pos =  pi + TMath::ATan2(-delta_y_pos,-delta_x_pos);
   Double_t alpha_neg =  pi + TMath::ATan2(-delta_y_neg,-delta_x_neg);

   Double_t vertex_x_neg =  xneg +  TMath::Abs(negtrackradius)*
   TMath::Cos(alpha_neg);
   Double_t vertex_y_neg =  yneg +  TMath::Abs(negtrackradius)*
   TMath::Sin(alpha_neg);

   Double_t vertex_x_pos =  xpos +  TMath::Abs(postrackradius)*
   TMath::Cos(alpha_pos);
   Double_t vertex_y_pos =  ypos +  TMath::Abs(postrackradius)*
   TMath::Sin(alpha_pos);

   Double_t x0neg =   helixneg[5];
   Double_t y0neg =   helixneg[0];

   Double_t x0pos =   helixpos[5];
   Double_t y0pos =   helixpos[0];

   Double_t d_neg = TMath::Sqrt((vertex_x_neg -  x0neg)*(vertex_x_neg - x0neg)
                               +(vertex_y_neg -  y0neg)*(vertex_y_neg - y0neg));

   Double_t d_pos = TMath::Sqrt((vertex_x_pos -  x0pos)*(vertex_x_pos - x0pos)
                               +(vertex_y_pos -  y0pos)*(vertex_y_pos - y0pos));

   Double_t r_neg =  TMath::Sqrt(negtrackradius*negtrackradius -
   d_neg*d_neg/4.);

   Double_t r_pos = TMath::Sqrt(postrackradius*postrackradius -
   d_pos*d_pos/4.);

   Double_t deltabeta_neg =  2*(pi +   TMath::ATan2(-d_neg/2.,-r_neg));
   Double_t deltabeta_pos = 2*(pi + TMath::ATan2(-d_pos/2.,-r_pos));

   Double_t delta_U_neg = negtrackradius*deltabeta_neg;
   Double_t delta_U_pos = postrackradius*deltabeta_pos;

   Double_t zphase_neg = ntrack->GetZ() +  delta_U_neg * ntrack->GetTgl();
   Double_t zphase_pos = ptrack->GetZ() +  delta_U_pos * ptrack->GetTgl();

   Double_t convposz =
   (zphase_pos*negtrackradius+zphase_neg*postrackradius)/(negtrackradius+postrackradius);

   return convposz;
}
