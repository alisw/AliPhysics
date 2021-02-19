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
 ************************************v**************************************/

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TMath.h"
#include "Riostream.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskCorrelationhhK0s.h"
#include "AliPIDResponse.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisCorrelationEventCollection.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliVVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"


class AliAnalysisTaskCorrelationhhK0s;   
using namespace std;          
ClassImp(AliAnalysisTaskCorrelationhhK0s) 

AliAnalysisTaskCorrelationhhK0s::AliAnalysisTaskCorrelationhhK0s() :AliAnalysisTaskSE(), 
  fAnalysisType("AOD"), 
  fCollidingSystem("pp"), 
  fAOD(0), 
  fPIDResponse(0),
  fEventCuts(0), 			  			
  fOutputList(0), 
  fSignalTree(0), 
  fBkgTree(0), 
  fOutputList2(0),
  fOutputList3(0),  
  fOutputList4(0),  
  fMCEvent(0), 
  fReadMCTruth(0),
  fIshhCorr(0),
  isEfficiency(0),
  isHybridMCTruth(0),
  fEventColl(0x0), 
  fEvt(0x0), 
  fzVertexBins(10), 
  fnMultBins(20),	 
  fMaxFirstMult(50),
  fMaxSecondMult(150),
  fnEventsToMix(50),
  fEtaTrigger(0.8),
  fEtahAssoc(0.8),
  fEtaV0Assoc(0.8),
  fFilterBitValue(128),
  fYear(2016),
  fHistPt(0), 
  fHistDCAxym1(0),
  fHistDCAzm1(0),
  fHistDCAxym2(0),
  fHistDCAzm2(0),
  fHistPtV0(0), 
  fHistPthAssoc(0), 
  fHistPtTMaxBefAllCfrDataMC(0),
  fHistPtTMinBefAll(0),
  fHistPtTMinBefAllMC(0),
  fHistPtvsMult(0), 
  fHistPtvsMultBefAll(0), 
  fHistPtMaxvsMult(0), 
  fHistPtMaxvsMultBefAll(0), 
  fHistPtMaxvsMultBefAllReco(0), 
  fHistPtMaxvsMultBefAllGen(0), 
  fHistZvertex(0),  
  fHistFractionSharedTPCClusters(0),
  fHistNumberChargedAllEvents(0),
  fHistNumberChargedNoTrigger(0),
  fHistNumberChargedTrigger(0),
  fHist_eta_phi(0),  
  fHist_eta_phi_PtMax(0),  
  fHist_multiplicity(0),
  fHist_multiplicity_EvwTrigger(0),
  fHistEventMult(0), 
  fHistEventV0(0), 
  fHistTrack(0),
  fHistTOFBunchCrossing(0), 
  fHistLengthvsCrossedRowsAfterSel(0),
  fHistLengthvsCrossedRows(0),
  fHistIsCommonParton(0),
  fHistCommonParton(0),
  fHistCommonPartonTrueCasc(0),
  fHistTriggerComposition(0), 
  fHistTriggerCompositionMCTruth(0), 
  fHistAssocComposition(0), 
  fHistAssocCompositionMCTruth(0), 
  fHistTrackAssoc(0), 
  fHistPDG(0),
  fHistTrackBufferOverflow(0), 
  fHistSecondParticleAll(0),
  fHistSecondParticleTruthAll(0),
  fHistSecondParticle(0),
  fHistSecondParticleTruth(0),
  fMassV0(0), 
  fDCAxyDaughters(0),
  fDCAzDaughters(0),
  fDCAxyDaughtersBis(0),
  fDCAzDaughtersBis(0),
  fDCAxyPos(0),
  fDCAzPos(0),
  fDCAxyNeg(0),
  fDCAzNeg(0),    
  fHistMultvsV0All(0),
  fHistMultvsV0AllTruth(0),
  fHistMultvsV0MCAll(0),
  fHistTriggerNotLeading(0),
  fHistTriggerNotLeadingMC(0),
  fHistMassvsPt_tagli(0x0),
  fHistMultvsTriggerBefAll(0),
  fHistMultvsTriggerMCTruthBefAll(0),
  fHistMultvsTrigger(0),
  fHistMultvsTriggerMCTruth(0),
  fHistMassPhoton(0),
  fHistMass2Photon(0),
  fHistPtArmvsAlpha(0),
  fHistPtArmvsAlphaAfterSelection(0),
  fHistPtArmvsAlphaAfterPhotonSelection(0),
  fHistPtArmvsAlphaAfterLambdaRejectionSelection(0),
  fHistTrigger(0),
  fHistTriggerMCTruth(0),
  fHistTriggerwV0(0),
  fHistTriggerwV0MCTruth(0),
  fHistMultiplicityVsVertexZ(0),
  fHistTriggervsMult(0),
  fHistTriggervsMultMC(0),
  fHistMultiplicityOfMixedEvent(0),
  fHistGeneratedTriggerPtPhi(0),
  fHistSelectedTriggerPtPhi(0),
  fHistSelectedGenTriggerPtPhi(0),
  fHistGeneratedV0PtTMaxPhi(0),
  fHistCPGeneratedV0PtTMaxPhi(0),
  fHistSelectedV0PtTMaxPhi(0),
  fHistGeneratedTriggerPtEta(0),
  fHistSelectedTriggerPtEta(0),
  fHistSelectedGenTriggerPtEta(0),
  fHistGeneratedV0PtTMaxEta(0),
  fHistCPGeneratedV0PtTMaxEta(0),
  fHistSelectedV0PtTMaxEta(0),
  fHistGeneratedV0PtPtTMax(0),
  fHistCPGeneratedV0PtPtTMax(0),
  fHistSelectedV0PtPtTMax(0),
  fHistSelectedGenV0PtPtTMax(0),
  fHistGeneratedV0PtEta(0),
  fHistSelectedV0PtEta(0),
  fHistReconstructedV0PtMass(0),
  fHistSelectedV0PtMass(0),
  fHistTriggerPtRecovsPtGen(0),
  fHistAssocPtRecovsPtGen(0),
  fHistTriggerPtRecovsPtGenNotPrim(0),
  fHistAssocPtRecovsPtGenNotPrim(0),
  fHistTriggerPtRecovsPtGenPion(0),
  fHistAssocPtRecovsPtGenPion(0),
  fHistTriggerPtRecovsPtGenProton(0),
  fHistAssocPtRecovsPtGenProton(0),
  fHistTriggerPtRecovsPtGenKaon(0),
  fHistAssocPtRecovsPtGenKaon(0),
  fHistResolutionTriggerPt(0),
  fHistResolutionTriggerPhi(0),
  fHistResolutionTriggerEta(0),
  fHistResolutionTriggerPhiPhi(0),
  fHistResolutionTriggerPhiEta(0),
  fHistResolutionV0Pt(0),
  fHistResolutionV0Phi(0),
  fHistResolutionV0PhivsPt(0),
  fHistResolutionV0Eta(0),
  fHistResolutionV0PtvsPt(0),
  fHistPrimaryTrigger(0),
  fHistPrimaryV0(0),
  fminPtj(2), 
  fmaxPtj(10), 
  fV0("K0s"),  
  fminPtV0(0), 
  fmaxPtV0(30),  
  fminPthAssoc(0), 
  fmaxPthAssoc(30),    
  Evcounter(0), 
  Evcounterczero(0),
  fmolt(5),
  farrGT(0), 
  fTrackBufferSize(20200),
  fTreeVariablePtTrigger(0),		      
  fTreeVariableChargeTrigger(0),		      
  fTreeVariableEtaTrigger(0), 		      
  fTreeVariablePhiTrigger(0),		      
  fTreeVariableDCAz(0),			      
  fTreeVariableDCAxy(0),
  fTreeVariableChargeAssoc(0),     
  fTreeVariableSkipAssoc(0),     
//  fTreeVariableEtaDaughterPos(0),
//  fTreeVariableEtaDaughterNeg(0),
  fTreeVariableIsCommonParton(0),
  fTreeVariablePdgCommonParton(0),
  fTreeVariableAssocDCAz(0),
  fTreeVariableAssocDCAxy(0),  
  fTreeVariableisPrimaryTrigger(0),  
  fTreeVariableisPrimaryV0(0),  
  fTreeVariableRapK0Short(0),		      	      
  fTreeVariableDcaV0ToPrimVertex (0),	      	      
  fTreeVariableDcaPosToPrimVertex(0),	      	      
  fTreeVariableDcaNegToPrimVertex(0),	      	      
  fTreeVariableV0CosineOfPointingAngle(0),      	      
  fTreeVariablePtV0(0),			      
  fTreeVariablectau(0),			      
  fTreeVariableInvMassK0s(0),		      
  fTreeVariableInvMassLambda(0),		      
  fTreeVariableInvMassAntiLambda(0),		      
  fTreeVariableEtaV0(0),			      
  fTreeVariablePhiV0(0),			      
  fTreeVariablePtArmenteros(0),                   
  fTreeVariableAlpha(0),
  fTreeVariableDeltaEta(0),			       
  fTreeVariableDeltaPhi(0),
  fTreeVariableDeltaTheta(0),			       
  fTreeVariableMultiplicity(0),                   
  fTreeVariableZvertex(0),
  fTreeVariablePDGCodeTrigger(0),
  fTreeVariablePDGCodeAssoc(0),
  FifoShiftok(kFALSE)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purpos, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorrelationhhK0s::AliAnalysisTaskCorrelationhhK0s(const char* name) : AliAnalysisTaskSE(name),
								 fAnalysisType("AOD"), 
								 fCollidingSystem("pp"), 
								 fAOD(0), 
								 fPIDResponse(0),
								 fEventCuts(0),								 
								 fOutputList(0), 
								 fSignalTree(0), 
								 fBkgTree(0), 
								 fOutputList2(0), 
								 fOutputList3(0), 
								 fOutputList4(0), 
								 fMCEvent(0), 
								 fReadMCTruth(0),
								 fIshhCorr(0),
								 isEfficiency(0),
                                                                 isHybridMCTruth(0),
								 fEventColl(0x0), 
								 fEvt(0x0), 
								 fzVertexBins(10), 
  fnMultBins(20),	 
  fMaxFirstMult(50),
  fMaxSecondMult(150),
  fnEventsToMix(50),
  fEtaTrigger(0.8),
  fEtahAssoc(0.8),
  fEtaV0Assoc(0.8),
  fFilterBitValue(128),
  fYear(2016),
  fHistPt(0), 
  fHistDCAxym1(0),
  fHistDCAzm1(0),
  fHistDCAxym2(0),
  fHistDCAzm2(0),
  fHistPtV0(0), 
  fHistPthAssoc(0), 
  fHistPtTMaxBefAllCfrDataMC(0), 
  fHistPtTMinBefAll(0),
  fHistPtTMinBefAllMC(0),
  fHistPtvsMult(0), 
  fHistPtvsMultBefAll(0), 
  fHistPtMaxvsMult(0), 
  fHistPtMaxvsMultBefAll(0), 
  fHistPtMaxvsMultBefAllReco(0), 
  fHistPtMaxvsMultBefAllGen(0), 
  fHistZvertex(0),  
  fHistFractionSharedTPCClusters(0),
  fHistNumberChargedAllEvents(0),
  fHistNumberChargedNoTrigger(0),
  fHistNumberChargedTrigger(0),
  fHist_eta_phi(0),  
  fHist_eta_phi_PtMax(0),  
  fHist_multiplicity(0),
  fHist_multiplicity_EvwTrigger(0),
  fHistEventMult(0), 
  fHistEventV0(0), 
  fHistTrack(0),
  fHistTOFBunchCrossing(0), 
  fHistLengthvsCrossedRowsAfterSel(0),
  fHistLengthvsCrossedRows(0),
  fHistIsCommonParton(0),
  fHistCommonParton(0),
  fHistCommonPartonTrueCasc(0),
  fHistTriggerComposition(0), 
  fHistTriggerCompositionMCTruth(0), 
  fHistAssocComposition(0), 
  fHistAssocCompositionMCTruth(0), 
  fHistTrackAssoc(0), 
  fHistPDG(0), 
  fHistTrackBufferOverflow(0), 
  fHistSecondParticleAll(0),
  fHistSecondParticleTruthAll(0),
  fHistSecondParticle(0),
  fHistSecondParticleTruth(0),
  fMassV0(0),
  fDCAxyDaughters(0),
  fDCAzDaughters(0),
  fDCAxyDaughtersBis(0),
  fDCAzDaughtersBis(0),
  fDCAxyPos(0),
  fDCAzPos(0),
  fDCAxyNeg(0),
  fDCAzNeg(0),    
  fHistMultvsV0All(0),
  fHistMultvsV0AllTruth(0),
  fHistMultvsV0MCAll(0), 
  fHistTriggerNotLeading(0),
  fHistTriggerNotLeadingMC(0),
  fHistMassvsPt_tagli(0x0),
  fHistMultvsTriggerBefAll(0),
  fHistMultvsTriggerMCTruthBefAll(0),
  fHistMultvsTrigger(0),
  fHistMultvsTriggerMCTruth(0),
  fHistMassPhoton(0),
  fHistMass2Photon(0),
  fHistPtArmvsAlpha(0),
  fHistPtArmvsAlphaAfterSelection(0),
  fHistPtArmvsAlphaAfterPhotonSelection(0),
  fHistPtArmvsAlphaAfterLambdaRejectionSelection(0),
  fHistTrigger(0),
  fHistTriggerMCTruth(0),
  fHistTriggerwV0(0),
  fHistTriggerwV0MCTruth(0),
  fHistMultiplicityVsVertexZ(0),
  fHistTriggervsMult(0),
  fHistTriggervsMultMC(0),
  fHistMultiplicityOfMixedEvent(0),
  fHistGeneratedTriggerPtPhi(0),
  fHistSelectedTriggerPtPhi(0),
  fHistSelectedGenTriggerPtPhi(0),
  fHistGeneratedV0PtTMaxPhi(0),
  fHistCPGeneratedV0PtTMaxPhi(0),
  fHistSelectedV0PtTMaxPhi(0),
  fHistGeneratedTriggerPtEta(0),
  fHistSelectedTriggerPtEta(0),
  fHistSelectedGenTriggerPtEta(0),
  fHistGeneratedV0PtTMaxEta(0),
  fHistCPGeneratedV0PtTMaxEta(0),
  fHistSelectedV0PtTMaxEta(0),
  fHistGeneratedV0PtPtTMax(0),
  fHistCPGeneratedV0PtPtTMax(0),
  fHistSelectedV0PtPtTMax(0),
  fHistSelectedGenV0PtPtTMax(0),
  fHistReconstructedV0PtMass(0),
  fHistGeneratedV0PtEta(0),
  fHistSelectedV0PtEta(0),
  fHistSelectedV0PtMass(0),
  fHistTriggerPtRecovsPtGen(0),
  fHistAssocPtRecovsPtGen(0),
  fHistTriggerPtRecovsPtGenNotPrim(0),
  fHistAssocPtRecovsPtGenNotPrim(0),
  fHistTriggerPtRecovsPtGenPion(0),
  fHistAssocPtRecovsPtGenPion(0),
  fHistTriggerPtRecovsPtGenProton(0),
  fHistAssocPtRecovsPtGenProton(0),
  fHistTriggerPtRecovsPtGenKaon(0),
  fHistAssocPtRecovsPtGenKaon(0),
  fHistResolutionTriggerPt(0),
  fHistResolutionTriggerPhi(0),
  fHistResolutionTriggerEta(0),
  fHistResolutionTriggerPhiPhi(0),
  fHistResolutionTriggerPhiEta(0),
  fHistResolutionV0Pt(0),
  fHistResolutionV0Phi(0),
  fHistResolutionV0PhivsPt(0),
  fHistResolutionV0Eta(0),
  fHistResolutionV0PtvsPt(0),
  fHistPrimaryTrigger(0),
  fHistPrimaryV0(0),
  fminPtj(2), 
  fmaxPtj(10), 
  fV0("K0s"),  
  fminPtV0(0), 
  fmaxPtV0(30),  
  fminPthAssoc(0), 
  fmaxPthAssoc(30),    
  Evcounter(0), 
  Evcounterczero(0),
  fmolt(5),
  farrGT(0), 
  fTrackBufferSize(20200),
  fTreeVariablePtTrigger(0),		      
  fTreeVariableChargeTrigger(0),		      
  fTreeVariableEtaTrigger(0), 		      
  fTreeVariablePhiTrigger(0),		      
  fTreeVariableDCAz(0),			      
  fTreeVariableDCAxy(0),		
  fTreeVariableChargeAssoc(0),     
  fTreeVariableSkipAssoc(0),     
//  fTreeVariableEtaDaughterPos(0),
//  fTreeVariableEtaDaughterNeg(0),
  fTreeVariableIsCommonParton(0),
  fTreeVariablePdgCommonParton(0),
  fTreeVariableAssocDCAz(0),
  fTreeVariableAssocDCAxy(0),  	      
  fTreeVariableisPrimaryTrigger(0),
  fTreeVariableisPrimaryV0(0),
  fTreeVariableRapK0Short(0),		      	      
  fTreeVariableDcaV0ToPrimVertex (0),	      	      
  fTreeVariableDcaPosToPrimVertex(0),	      	      
  fTreeVariableDcaNegToPrimVertex(0),	      	      
  fTreeVariableV0CosineOfPointingAngle(0),      	      
  fTreeVariablePtV0(0),			      
  fTreeVariablectau(0),			      
  fTreeVariableInvMassK0s(0),		      
  fTreeVariableInvMassLambda(0),		      
  fTreeVariableInvMassAntiLambda(0),		      
  fTreeVariableEtaV0(0),			      
  fTreeVariablePhiV0(0),			      
  fTreeVariablePtArmenteros(0),                   
  fTreeVariableAlpha(0),
  fTreeVariableDeltaEta(0),			       
  fTreeVariableDeltaPhi(0),
  fTreeVariableDeltaTheta(0),			       
  fTreeVariableMultiplicity(0),                   
  fTreeVariableZvertex(0),
  fTreeVariablePDGCodeTrigger(0),
  fTreeVariablePDGCodeAssoc(0),
  FifoShiftok(kFALSE)
{
                      
  DefineInput(0, TChain::Class()); 
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());  
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TList::Class());  
  DefineOutput(5, TList::Class());  
  DefineOutput(6, TList::Class());  
}
//_____________________________________________________________________________
AliAnalysisTaskCorrelationhhK0s::~AliAnalysisTaskCorrelationhhK0s()
{

  if(fOutputList) {
    delete fOutputList;  
  }
  if(fSignalTree) {
    delete fSignalTree;
  }
  if(fBkgTree) {
    delete fBkgTree;   
  }
  if(fOutputList2) {
    delete fOutputList2;
  }
  if(fOutputList3) {
    delete fOutputList3;
  }
  if(fOutputList4) {
    delete fOutputList4;
  }
  if (farrGT)
    delete[] farrGT;
  farrGT=0;

  for(unsigned short i=0; i < fzVertexBins; i++){
    for(unsigned short j=0; j < fnMultBins; j++){
      delete fEventColl[i][j];
    }
    delete[] fEventColl[i];
  }
  delete[] fEventColl;
  
}
  
void AliAnalysisTaskCorrelationhhK0s::ProcessMCParticles(Bool_t Generated, AliAODTrack *track, Int_t& labelPrimOrSec, Float_t lPercentiles, Bool_t isV0, Double_t ZAtDCA, Float_t PtTriggMax, Bool_t ishhCorr, Int_t VPdgTrig[], Int_t VParticleTrigLabel[])
{

  Float_t moltep[6]={0,5,10,30,50,100};  //V0M multiplicity intervals
  Int_t PDGCodeAssoc[2]={310, 3122};
  Int_t ParticleType =-999;
  Int_t PAntiP =-999;
  if (fV0=="K0s") {
    ParticleType=0;
    PAntiP=1;
  }
  else   if (fV0=="Lambda") {
    ParticleType=1;
  }

  AliAODMCParticle *VParticle[50]={0};
  Int_t VPdg[50]={0};
  Int_t VParticleLabel[50]={0};
  AliAODMCParticle *VParticleTrig[50]={0};

  TClonesArray* AODMCTrackArraybis =0x0;  
  AODMCTrackArraybis = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArraybis == NULL){
    return;
    Printf("ERROR: stack not available");
  }
  if(Generated){
    // Loop over all generated primary MC particle
    for(Long_t i = 0; i < AODMCTrackArraybis->GetEntriesFast(); i++) {
      
      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArraybis->At(i));
      if (!particle) continue;
      
      //common parton information
      VParticle[0]=particle;
      VPdg[0]=       VParticle[0]->GetPdgCode();
      VParticleLabel[0]=       VParticle[0]->GetLabel();

      for (Int_t i=0; i<50; i++){
        VParticleLabel[i+1]=VParticle[i]->GetMother();
        VParticle[i+1] = static_cast<AliAODMCParticle*>(AODMCTrackArraybis-> At(TMath::Abs(VParticleLabel[i+1])));
        VPdg[i+1] = VParticle[i+1]->GetPdgCode();
	if ((VParticleLabel[i] ==1 || VParticleLabel[i] ==-1) && (VPdg[i]==2212) && (VParticleLabel[i] == VParticleLabel[i+1])) {
	  break;
	}
      }
      //

      if(isV0==kFALSE){ //for trigger particles
	fHistPDG->Fill(particle->GetPdgCode());      
	if((particle->Charge())==0) continue;	
	if(TMath::Abs(particle->Eta())>fEtaTrigger)continue; 
	if (!(particle->IsPhysicalPrimary()))continue; 
	fHistGeneratedTriggerPtPhi->Fill(particle->Pt(), particle->Phi(), lPercentiles);
	fHistGeneratedTriggerPtEta->Fill(particle->Pt(), particle->Eta(), lPercentiles);
      }
      else if(isV0==kTRUE){ //for associated particles

	//let's check if the generated V0 comes from the same parton as the trigger particle                               
        Bool_t IsCommonParton=kFALSE;
        for (Int_t i=1; i<50; i++){ //I start from one since last element cannot be a parton but is a hadron                
          if (IsCommonParton==1) break;
          for (Int_t j=1; j<50; j++){
            if ((VParticleLabel[i] == VParticleTrigLabel[j] ) &&  VParticleTrigLabel[j]!=0 && ( TMath::Abs(VPdg[i]) <=8 || TMath::Abs(VPdg[i]) ==21)) { //both Xi and Trigger particle have a common ancestor which has to be a quark or a gluon-> ther\efore te cascade comes form the jet defined by the trigger particle                                 	      
	      IsCommonParton =1;
	      break;
            }
          }
	}

	if(!ishhCorr){ //if associated particles are K0s
	  if (TMath::Abs(particle->GetPdgCode())!=PDGCodeAssoc[ParticleType]) continue;
	  //	  if(TMath::Abs(particle->Eta())>fEtaV0Assoc) continue;
	  if (!(particle->IsPhysicalPrimary()))continue;
	  if (fV0=="Lambda" && TMath::Abs(particle->GetPdgCode())==PDGCodeAssoc[ParticleType]) PAntiP=1;
	  if (fV0=="Lambda" && TMath::Abs(particle->GetPdgCode())==-PDGCodeAssoc[ParticleType]) PAntiP=-1;

	}
	else if(ishhCorr){ //if associated particles are hadrons
	  if((particle->Charge())==0) continue;	
	  //	  if(TMath::Abs(particle->Eta())>fEtahAssoc)continue; 
	  if (!(particle->IsPhysicalPrimary()))continue; 
	  if ((particle->GetLabel()) == (static_cast<AliAODMCParticle*>(AODMCTrackArraybis->At(TMath::Abs(track->GetLabel()))))->GetLabel() ) continue;
	  if (particle->Pt() == (static_cast<AliAODMCParticle*>(AODMCTrackArraybis->At(TMath::Abs(track->GetLabel()))))->Pt() ) continue; 
	}


	Float_t EtaAssoc = fEtaV0Assoc;
	if (ishhCorr) EtaAssoc = fEtahAssoc;
	if(TMath::Abs(particle->Eta())<EtaAssoc){
	fHistGeneratedV0PtTMaxPhi[0]->Fill(PAntiP*PtTriggMax,particle->Phi(), lPercentiles );
	fHistGeneratedV0PtTMaxEta[0]->Fill(PAntiP*PtTriggMax,particle->Eta(), lPercentiles );
	fHistGeneratedV0PtPtTMax[0]->Fill(particle->Pt(),PAntiP*PtTriggMax, lPercentiles );
	fHistGeneratedV0PtEta[0]->Fill(particle->Pt(), particle->Eta(), lPercentiles );
	if (IsCommonParton){
	  fHistCPGeneratedV0PtTMaxPhi[0]->Fill(PAntiP*PtTriggMax,particle->Phi(), lPercentiles );
	  fHistCPGeneratedV0PtTMaxEta[0]->Fill(PAntiP*PtTriggMax,particle->Eta(), lPercentiles );
	  fHistCPGeneratedV0PtPtTMax[0]->Fill(particle->Pt(),PAntiP*PtTriggMax, lPercentiles );
	}
	}

	if (TMath::Abs(particle->Y())<0.5 ){
	fHistGeneratedV0PtTMaxPhi[1]->Fill(PAntiP*PtTriggMax,particle->Phi(), lPercentiles );
	fHistGeneratedV0PtTMaxEta[1]->Fill(PAntiP*PtTriggMax,particle->Eta(), lPercentiles );
	fHistGeneratedV0PtPtTMax[1]->Fill(particle->Pt(),PAntiP*PtTriggMax, lPercentiles );
	fHistGeneratedV0PtEta[1]->Fill(particle->Pt(), particle->Eta(), lPercentiles );
	if (IsCommonParton){
	  fHistCPGeneratedV0PtTMaxPhi[1]->Fill(PAntiP*PtTriggMax,particle->Phi(), lPercentiles );
	  fHistCPGeneratedV0PtTMaxEta[1]->Fill(PAntiP*PtTriggMax,particle->Eta(), lPercentiles );
	  fHistCPGeneratedV0PtPtTMax[1]->Fill(particle->Pt(),PAntiP*PtTriggMax, lPercentiles );
	}

	}
      }      
    }
  }
  else {
    // Loop over all reconstructed primary MC particle (here implemented only for trigger particles)
    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArraybis->At(TMath::Abs(track->GetLabel())));

    VParticleTrig[0]=particle;
    VPdgTrig[0]=       VParticleTrig[0]->GetPdgCode();
    VParticleTrigLabel[0]=       VParticleTrig[0]->GetLabel();

    for (Int_t i=0; i<50; i++){
      VParticleTrigLabel[i+1]=VParticleTrig[i]->GetMother();
      VParticleTrig[i+1] = static_cast<AliAODMCParticle*>(AODMCTrackArraybis-> At(TMath::Abs(VParticleTrigLabel[i+1])));
      VPdgTrig[i+1] = VParticleTrig[i+1]->GetPdgCode();
      //      cout << VPdgTrig[i] << " (" << VParticleTrigLabel[i] << ") " << "<-" ;                                       
      if ((VParticleTrigLabel[i] ==1 || VParticleTrigLabel[i] ==-1) && (VPdgTrig[i]==2212) && (VParticleTrigLabel[i] == VParticleTrigLabel[i+1] )) {
	//cout << endl;                                                                                                    
	break;
      }
    }

    Float_t TrackLength=0;
    TrackLength = GetLengthInActiveZone(track, /*1,*/ 2.0, 220.0, fAOD->GetMagneticField());

    if(!(particle->IsPhysicalPrimary()))     {

      fHistTriggerComposition->Fill(particle->GetPdgCode(),0.1, track->Pt());
      fHistTriggerPtRecovsPtGenNotPrim->Fill(particle->Pt(), track->Pt());

    }

    if(particle->IsPhysicalPrimary()){

      fHistTriggerComposition->Fill(particle->GetPdgCode(),1, track->Pt());
      fHistTriggerPtRecovsPtGen->Fill(particle->Pt(), track->Pt());
      if (particle->GetPdgCode()==211 || particle->GetPdgCode()==-211)      fHistTriggerPtRecovsPtGenPion->Fill(particle->Pt(), track->Pt());
      if (particle->GetPdgCode()==2212 || particle->GetPdgCode()==-2212)      fHistTriggerPtRecovsPtGenProton->Fill(particle->Pt(), track->Pt());
      if (particle->GetPdgCode()==321 || particle->GetPdgCode()==-321)      fHistTriggerPtRecovsPtGenKaon->Fill(particle->Pt(), track->Pt());

      fHistResolutionTriggerPt->Fill(track->Pt()- particle->Pt(),    track->Pt()); //PtTriggerMax
      fHistResolutionTriggerPhi->Fill(track->Phi()- particle->Phi(), track->Pt());
      fHistResolutionTriggerEta->Fill(track->Eta()- particle->Eta(), track->Pt());
      fHistResolutionTriggerPhiPhi->Fill(track->Phi()- particle->Phi(), track->Phi());
      fHistResolutionTriggerPhiEta->Fill(track->Phi()- particle->Phi(), track->Eta());
    
      if(  (TMath::Abs(ZAtDCA) < 1.)) {
	fHistSelectedTriggerPtPhi[0]->Fill(track->Pt(), track->Phi(), lPercentiles);
	fHistSelectedTriggerPtEta[0]->Fill(track->Pt(), track->Eta(), lPercentiles);    
	fHistSelectedGenTriggerPtEta[0]->Fill(particle->Pt(), particle->Eta(), lPercentiles);    
	fHistSelectedGenTriggerPtPhi[0]->Fill(particle->Pt(), particle->Phi(), lPercentiles);    
      }
      labelPrimOrSec=1;
    }
    else if(particle->IsSecondaryFromWeakDecay())      labelPrimOrSec=2;
    else if(particle->IsSecondaryFromMaterial())      labelPrimOrSec=3;
    else labelPrimOrSec=4;

    for (Int_t m =0; m<5;m++){
      if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	for(Int_t p=1; p<=4; p++){
	  if (labelPrimOrSec==p){
	    if( (TMath::Abs(ZAtDCA) < 1))   fHistPrimaryTrigger[m][0]->Fill(p,track->Pt() );
	    //	    if( (TMath::Abs(ZAtDCA) < 2))   fHistPrimaryTrigger[m][1]->Fill(p,track->Pt() );
	    //	    if( (TMath::Abs(ZAtDCA) < 0.5))   fHistPrimaryTrigger[m][2]->Fill(p,track->Pt() );
	  }
	}
      }
    }
    for(Int_t p=1; p<=4; p++){
      if (labelPrimOrSec==p){
	if( (TMath::Abs(ZAtDCA) < 1))   fHistPrimaryTrigger[5][0]->Fill(p,track->Pt() );
	//	if( (TMath::Abs(ZAtDCA) < 2))   fHistPrimaryTrigger[5][1]->Fill(p,track->Pt() );
	//	if( (TMath::Abs(ZAtDCA) < 0.5))   fHistPrimaryTrigger[5][2]->Fill(p,track->Pt() );
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskCorrelationhhK0s::UserCreateOutputObjects()
{

  fEventColl = new AliAnalysisCorrelationEventCollection **[fzVertexBins]; 
  
  for (unsigned short i=0; i<fzVertexBins; i++) {
    fEventColl[i] = new AliAnalysisCorrelationEventCollection *[fnMultBins];
    for (unsigned short j=0; j<fnMultBins; j++) {
      fEventColl[i][j] = new AliAnalysisCorrelationEventCollection(fnEventsToMix+1, fMaxFirstMult, fMaxSecondMult);
    }
  }

  // Store pointer to global tracks
  farrGT = new Int_t[fTrackBufferSize];

  fOutputList = new TList();         
  fOutputList->SetOwner(kTRUE);     
  fOutputList2 = new TList();         
  fOutputList2->SetOwner(kTRUE);     
  fOutputList3 = new TList();         
  fOutputList3->SetOwner(kTRUE);     
  fOutputList4 = new TList();         
  fOutputList4->SetOwner(kTRUE);     
  
  fSignalTree= new TTree("fSignalTree","fSignalTree");
  fSignalTree->Branch("fTreeVariablePtTrigger",          &fTreeVariablePtTrigger   , "fTreeVariablePtTrigger/D");
  fSignalTree->Branch("fTreeVariableChargeTrigger",      &fTreeVariableChargeTrigger, "fTreeVariableChargeTrigger/I");
  fSignalTree->Branch("fTreeVariableEtaTrigger",         &fTreeVariableEtaTrigger  , "fTreeVariableEtaTrigger/D");
  fSignalTree->Branch("fTreeVariablePhiTrigger",         &fTreeVariablePhiTrigger, "fTreeVariablePhiTrigger/D");
  fSignalTree->Branch("fTreeVariableDCAz",               &fTreeVariableDCAz  , "fTreeVariableDCAz/D");
  fSignalTree->Branch("fTreeVariableDCAxy",              &fTreeVariableDCAxy  , "fTreeVariableDCAxy/D");
  fSignalTree->Branch("fTreeVariableChargeAssoc",        &fTreeVariableChargeAssoc  , "fTreeVariableChargeAssoc/I");
  fSignalTree->Branch("fTreeVariableSkipAssoc",          &fTreeVariableSkipAssoc  , "fTreeVariableSkipAssoc/O");
  //  fSignalTree->Branch("fTreeVariableEtaDaughterPos",          &fTreeVariableEtaDaughterPos  , "fTreeVariableEtaDaughterPos/D");
  //  fSignalTree->Branch("fTreeVariableEtaDaughterNeg",          &fTreeVariableEtaDaughterNeg  , "fTreeVariableEtaDaughterNeg/D");
  fSignalTree->Branch("fTreeVariableIsCommonParton",     &fTreeVariableIsCommonParton, "fTreeVariableIsCommonParton/O");
  fSignalTree->Branch("fTreeVariablePdgCommonParton",    &fTreeVariablePdgCommonParton, "fTreeVariablePdgCommonParton/I");
  fSignalTree->Branch("fTreeVariableAssocDCAz",          &fTreeVariableAssocDCAz  , "fTreeVariableAssocDCAz/D");
  fSignalTree->Branch("fTreeVariableAssocDCAxy",         &fTreeVariableAssocDCAxy  , "fTreeVariableAssocDCAxy/D");
  fSignalTree->Branch("fTreeVariableisPrimaryTrigger",   &fTreeVariableisPrimaryTrigger  , "fTreeVariableisPrimaryTrigger/I");
  fSignalTree->Branch("fTreeVariableisPrimaryV0",        &fTreeVariableisPrimaryV0  , "fTreeVariableisPrimaryV0/I");
  fSignalTree->Branch("fTreeVariableRapK0Short",         &fTreeVariableRapK0Short               ,"fTreeVariableRapK0Short/D");
  fSignalTree->Branch("fTreeVariableDcaV0ToPrimVertex",  &fTreeVariableDcaV0ToPrimVertex 	, "fTreeVariableDcaV0ToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableDcaPosToPrimVertex", &fTreeVariableDcaPosToPrimVertex	, "fTreeVariableDcaPosToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableDcaNegToPrimVertex", &fTreeVariableDcaNegToPrimVertex	, "fTreeVariableDcaNegToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableV0CosineOfPointingAngle", &fTreeVariableV0CosineOfPointingAngle	, "fTreeVariableV0CosineOfPointingAngle/D");
  fSignalTree->Branch("fTreeVariablePtV0",               &fTreeVariablePtV0   , "fTreeVariablePtV0/D");
  fSignalTree->Branch("fTreeVariablectau",               &fTreeVariablectau   , "fTreeVariablectau/D");
  fSignalTree->Branch("fTreeVariableInvMassK0s",         &fTreeVariableInvMassK0s, "fTreeVariableInvMassK0s/D");
  fSignalTree->Branch("fTreeVariableInvMassLambda",      &fTreeVariableInvMassLambda, "fTreeVariableInvMassLambda/D");
  fSignalTree->Branch("fTreeVariableInvMassAntiLambda",  &fTreeVariableInvMassAntiLambda, "fTreeVariableInvMassAntiLambda/D");
  fSignalTree->Branch("fTreeVariableEtaV0",              &fTreeVariableEtaV0  , "fTreeVariableEtaV0/D");
  fSignalTree->Branch("fTreeVariablePhiV0",              &fTreeVariablePhiV0, "fTreeVariablePhiV0/D");
  fSignalTree->Branch("fTreeVariablePtArmenteros",       &fTreeVariablePtArmenteros  , "fTreeVariablePtArmenteros/D");
  fSignalTree->Branch("fTreeVariableAlpha",              &fTreeVariableAlpha  , "fTreeVariableAlpha/D");
  fSignalTree->Branch("fTreeVariableDeltaEta",           &fTreeVariableDeltaEta  , "fTreeVariableDeltaEta/D");
  fSignalTree->Branch("fTreeVariableDeltaPhi",           &fTreeVariableDeltaPhi, "fTreeVariableDeltaPhi/D");
  fSignalTree->Branch("fTreeVariableDeltaTheta",         &fTreeVariableDeltaTheta, "fTreeVariableDeltaTheta/D");
  fSignalTree->Branch("fTreeVariableMultiplicity",       &fTreeVariableMultiplicity , "fTreeVariableMultiplicity/D");
  fSignalTree->Branch("fTreeVariableZvertex",            &fTreeVariableZvertex  , "fTreeVariableZvertex/D");
  fSignalTree->Branch("fTreeVariablePDGCodeTrigger",     &fTreeVariablePDGCodeTrigger  , "fTreeVariablePDGCodeTrigger/I");
  fSignalTree->Branch("fTreeVariablePDGCodeAssoc",       &fTreeVariablePDGCodeAssoc  , "fTreeVariablePDGCodeAssoc/I");


  fBkgTree= new TTree("fBkgTree","fBkgTree");
  fBkgTree->Branch("fTreeVariablePtTrigger",          &fTreeVariablePtTrigger   , "fTreeVariablePtTrigger/D");
  fBkgTree->Branch("fTreeVariableChargeTrigger",      &fTreeVariableChargeTrigger, "fTreeVariableChargeTrigger/I");
  fBkgTree->Branch("fTreeVariableEtaTrigger",         &fTreeVariableEtaTrigger  , "fTreeVariableEtaTrigger/D");
  fBkgTree->Branch("fTreeVariablePhiTrigger",         &fTreeVariablePhiTrigger, "fTreeVariablePhiTrigger/D");
  fBkgTree->Branch("fTreeVariableDCAz",               &fTreeVariableDCAz  , "fTreeVariableDCAz/D");
  fBkgTree->Branch("fTreeVariableDCAxy",              &fTreeVariableDCAxy  , "fTreeVariableDCAxy/D");
  fBkgTree->Branch("fTreeVariableChargeAssoc",        &fTreeVariableChargeAssoc  , "fTreeVariableChargeAssoc/I");
  fBkgTree->Branch("fTreeVariableSkipAssoc",          &fTreeVariableSkipAssoc  , "fTreeVariableSkipAssoc/O");
  fBkgTree->Branch("fTreeVariableIsCommonParton",     &fTreeVariableIsCommonParton, "fTreeVariableIsCommonParton/O");
  fBkgTree->Branch("fTreeVariablePdgCommonParton",    &fTreeVariablePdgCommonParton, "fTreeVariablePdgCommonParton/I");
  fBkgTree->Branch("fTreeVariableAssocDCAz",          &fTreeVariableAssocDCAz  , "fTreeVariableAssocDCAz/D");
  fBkgTree->Branch("fTreeVariableAssocDCAxy",         &fTreeVariableAssocDCAxy  , "fTreeVariableAssocDCAxy/D");
  fBkgTree->Branch("fTreeVariableisPrimaryTrigger",   &fTreeVariableisPrimaryTrigger  , "fTreeVariableisPrimaryTrigger/I");  
  fBkgTree->Branch("fTreeVariableisPrimaryV0",        &fTreeVariableisPrimaryV0  , "fTreeVariableisPrimaryV0/I");  
  fBkgTree->Branch("fTreeVariableRapK0Short",         &fTreeVariableRapK0Short               ,"fTreeVariableRapK0Short/D");
  fBkgTree->Branch("fTreeVariableDcaV0ToPrimVertex",  &fTreeVariableDcaV0ToPrimVertex 	, "fTreeVariableDcaV0ToPrimVertex/D");
  fBkgTree->Branch("fTreeVariableDcaPosToPrimVertex", &fTreeVariableDcaPosToPrimVertex	, "fTreeVariableDcaPosToPrimVertex/D");
  fBkgTree->Branch("fTreeVariableDcaNegToPrimVertex", &fTreeVariableDcaNegToPrimVertex	, "fTreeVariableDcaNegToPrimVertex/D");
  fBkgTree->Branch("fTreeVariableV0CosineOfPointingAngle", &fTreeVariableV0CosineOfPointingAngle	, "fTreeVariableV0CosineOfPointingAngle/D");
  fBkgTree->Branch("fTreeVariablePtV0",               &fTreeVariablePtV0   , "fTreeVariablePtV0/D");
  fBkgTree->Branch("fTreeVariablectau",               &fTreeVariablectau   , "fTreeVariablectau/D");
  fBkgTree->Branch("fTreeVariableInvMassK0s",         &fTreeVariableInvMassK0s, "fTreeVariableInvMassK0s/D");
  fBkgTree->Branch("fTreeVariableInvMassLambda",      &fTreeVariableInvMassLambda, "fTreeVariableInvMassLambda/D");
  fBkgTree->Branch("fTreeVariableInvMassAntiLambda",  &fTreeVariableInvMassAntiLambda, "fTreeVariableInvMassAntiLambda/D");
  fBkgTree->Branch("fTreeVariableEtaV0",              &fTreeVariableEtaV0  , "fTreeVariableEtaV0/D");
  fBkgTree->Branch("fTreeVariablePhiV0",              &fTreeVariablePhiV0, "fTreeVariablePhiV0/D");
  fBkgTree->Branch("fTreeVariablePtArmenteros",       &fTreeVariablePtArmenteros  , "fTreeVariablePtArmenteros/D");
  fBkgTree->Branch("fTreeVariableAlpha",              &fTreeVariableAlpha  , "fTreeVariableAlpha/D");
  fBkgTree->Branch("fTreeVariableDeltaEta",           &fTreeVariableDeltaEta  , "fTreeVariableDeltaEta/D");
  fBkgTree->Branch("fTreeVariableDeltaPhi",           &fTreeVariableDeltaPhi  , "fTreeVariableDeltaPhi/D");
  fBkgTree->Branch("fTreeVariableDeltaTheta",         &fTreeVariableDeltaTheta, "fTreeVariableDeltaTheta/D");
  fBkgTree->Branch("fTreeVariableMultiplicity",       &fTreeVariableMultiplicity , "fTreeVariableMultiplicity/D");
  fBkgTree->Branch("fTreeVariableZvertex",            &fTreeVariableZvertex  , "fTreeVariableZvertex/D");
  fBkgTree->Branch("fTreeVariablePDGCodeTrigger",     &fTreeVariablePDGCodeTrigger  , "fTreeVariablePDGCodeTrigger/I");
  fBkgTree->Branch("fTreeVariablePDGCodeAssoc",       &fTreeVariablePDGCodeAssoc  , "fTreeVariablePDGCodeAssoc/I");
 
  fHistPt = new TH1F("fHistPt", "p_{T} distribution of selected charged tracks in events used for AC", 300, 0, 30); 
  fHistPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 
  fHistDCAxym1 = new TH1F("fHistDCAxym1", "DCAxy method 1 before DCA cuts", 1000, -10, 10); 
  fHistDCAxym1->GetXaxis()->SetTitle("DCAxy method 1 (cm)");
  fHistDCAxym2 = new TH1F("fHistDCAxym2", "DCAxy method 2 before DCA cuts", 1000, -10, 10); 
  fHistDCAxym2->GetXaxis()->SetTitle("DCAxy method 2 (cm)");
  fHistDCAzm1 = new TH1F("fHistDCAzm1", "DCAz method 1 before DCA cuts", 1000, -10, 10); 
  fHistDCAzm1->GetXaxis()->SetTitle("DCAz method 1 (cm)");
  fHistDCAzm2 = new TH1F("fHistDCAzm2", "DCAz method 2 before DCA cuts", 1000, -10, 10); 
  fHistDCAzm2->GetXaxis()->SetTitle("DCAz method 2 (cm)");
 
  fHistPtV0 = new TH1F("fHistPtV0", "p_{T} distribution of selected V0 in events used for AC", 300, 0, 30); 
  fHistPtV0->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPthAssoc = new TH1F("fHistPthAssoc", "p_{T} distribution of selected associated charged particles in events used for AC", 300, 0, 30); 
  fHistPthAssoc->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMaxBefAllCfrDataMC= new TH2F("fHistPtTMaxBefAllCfrDataMC", "p_{T} distribution of trigger particle with Pt maximum in the event", 300, 0, 30, 300, 0,30); 
  fHistPtTMaxBefAllCfrDataMC->GetXaxis()->SetTitle("p^{Trigger, Max}_{T} (reco)(GeV/c)");
  fHistPtTMaxBefAllCfrDataMC->GetYaxis()->SetTitle("p^{Trigger, Max}_{T} (MC)(GeV/c)");

  fHistPtTMinBefAll = new TH1F("fHistPtTMinBefAll", "p_{T} distribution of reco trigger particle with Pt minimum in the event", 300, 0, 30); 
  fHistPtTMinBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMinBefAllMC = new TH1F("fHistPtTMinBefAllMC", "p_{T} distribution of true trigger particle with Pt minimum in the event", 300, 0, 30); 
  fHistPtTMinBefAllMC->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtvsMultBefAll= new TH2F("fHistPtvsMultBefAll", "p_{T} and centrality distribution of charged tracks in events w T>0", 300, 0, 30, 100, 0, 100); 
  fHistPtvsMultBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtvsMultBefAll->GetYaxis()->SetTitle("Centrality");

  fHistPtvsMult= new TH2F("fHistPtvsMult", "p_{T} and centrality distribution of charged tracks in events used for AC", 300, 0, 30, 100, 0, 100); 
  fHistPtvsMult->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtvsMult->GetYaxis()->SetTitle("Centrality");

  fHistPtMaxvsMultBefAll= new TH2F("fHistPtMaxvsMultBefAll", "p_{T} and centrality distribution of charged tracks with maxiumum pt in events w T>0", 600, 0, 30, 100, 0, 100); 
  fHistPtMaxvsMultBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtMaxvsMultBefAll->GetYaxis()->SetTitle("Centrality");

  fHistPtMaxvsMultBefAllGen= new TH2F("fHistPtMaxvsMultBefAllGen", "p_{T} and centrality distribution of charged tracks with maxiumum pt in events w T>0", 600, 0, 30, 100, 0, 100); 
  fHistPtMaxvsMultBefAllGen->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtMaxvsMultBefAllGen->GetYaxis()->SetTitle("Centrality");

  fHistPtMaxvsMultBefAllReco= new TH2F("fHistPtMaxvsMultBefAllReco", "p_{T} and centrality distribution of charged tracks with maxiumum pt in events w T>0", 600, 0, 30, 100, 0, 100); 
  fHistPtMaxvsMultBefAllReco->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtMaxvsMultBefAllReco->GetYaxis()->SetTitle("Centrality");

  fHistPtMaxvsMult= new TH2F("fHistPtMaxvsMult", "p_{T} and centrality distribution of charged tracks with maximum pT in events used for AC)", 600, 0, 30, 100, 0, 100); 
  fHistPtMaxvsMult->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtMaxvsMult->GetYaxis()->SetTitle("Centrality");

  fHistZvertex= new TH1F("fHistZvertex", "Z vertex distribution of selected events used for AC", 40,-20,20);

  fHistFractionSharedTPCClusters = new TH1F ("fHistFractionSharedTPCClusters", "fHistFractionSharedTPCClusters",100, 0,1);

  fHistNumberChargedAllEvents=new TH3F("fHistNumberChargedAllEvents", "fHistNumberChargedAllEvents", 100,0,100, 100,0,100, 60, 0,30);
  fHistNumberChargedAllEvents->GetXaxis()->SetTitle("Multiplicity class");
  fHistNumberChargedAllEvents->GetYaxis()->SetTitle("Number of charged primary particles");
  fHistNumberChargedAllEvents->GetZaxis()->SetTitle("p^{Trigg, Max}_{T} (GeV/c)");

  fHistNumberChargedNoTrigger=new TH3F("fHistNumberChargedNoTrigger", "fHistNumberChargedNoTrigger", 100,0,100, 100,0,100, 60, 0,30);
  fHistNumberChargedNoTrigger->GetXaxis()->SetTitle("Multiplicity class");
  fHistNumberChargedNoTrigger->GetYaxis()->SetTitle("Number of charged primary particles");
  fHistNumberChargedNoTrigger->GetZaxis()->SetTitle("p^{Trigg, Max}_{T} (GeV/c)");

  fHistNumberChargedTrigger=new TH3F("fHistNumberChargedTrigger", "fHistNumberChargedTrigger", 100,0,100, 100,0,100, 60, 0,30);
  fHistNumberChargedTrigger->GetXaxis()->SetTitle("Multiplicity class");
  fHistNumberChargedTrigger->GetYaxis()->SetTitle("Number of charged primary particles");
  fHistNumberChargedTrigger->GetZaxis()->SetTitle("p^{Trigg, Max}_{T} (GeV/c)");

  fHist_eta_phi= new TH2F("fHist_eta_phi", "Distribution of charged tracks in events used for AC", 400, 0,2*TMath::Pi(), 400,-0.8, 0.8);
  fHist_eta_phi->GetYaxis()->SetTitle("Eta");
  fHist_eta_phi->GetXaxis()->SetTitle("Phi (radians)"); 

  fHist_eta_phi_PtMax= new TH2F("fHist_eta_phi_PtMax", "Distribution of charged tracks with maximum Pt in events used for AC", 400, 0,2*TMath::Pi(), 400,-0.8, 0.8);
  fHist_eta_phi_PtMax->GetYaxis()->SetTitle("Eta");
  fHist_eta_phi_PtMax->GetXaxis()->SetTitle("Phi (radians)"); 

  fHist_multiplicity=new TH1F("fHist_multiplicity", "fHist_multiplicity", 100, 0, 100); 
  fHist_multiplicity->SetTitle("Centrality distribution of events used for AC");
  fHist_multiplicity_EvwTrigger= new TH1F("fHist_multiplicity_EvwTrigger", "fHist_multiplicity_EvwTrigger", 100, 0, 100); 
  fHist_multiplicity_EvwTrigger->SetTitle("Centrality distribution of events with NT>0");

  fHistPDG=new TH1F("fHistPDG", "fHistPDG",6400, -3200, 3200);

  fHistTrackBufferOverflow = new TH1F("fHistTrackBufferOverflow","",2,0,2);
  
  fHistEventMult=new TH1F("fHistEventMult", "fHistEventMult", 24, 0.5, 24.5);
  fHistEventMult->GetXaxis()->SetBinLabel(1,"Event cuts");
  fHistEventMult->GetXaxis()->SetBinLabel(2,"Events w/PV and AOD");
  fHistEventMult->GetXaxis()->SetBinLabel(3,"Events w/|Vx|<10 cm"); 
  fHistEventMult->GetXaxis()->SetBinLabel(4,"Events w/ PID"); 
  fHistEventMult->GetXaxis()->SetBinLabel(5,"centrality <= 199"); 
  fHistEventMult->GetXaxis()->SetBinLabel(6,"NO PILE UP"); 
  fHistEventMult->GetXaxis()->SetBinLabel(7,"INT7"); 
  fHistEventMult->GetXaxis()->SetBinLabel(8,"ANY"); 
  fHistEventMult->GetXaxis()->SetBinLabel(9,"isSelected"); 
  fHistEventMult->GetXaxis()->SetBinLabel(10,"Ntrigger>0"); 
  fHistEventMult->GetXaxis()->SetBinLabel(11,"Ntrigger>0 (MC)"); 
  fHistEventMult->GetXaxis()->SetBinLabel(12,"NTrigger>1");
  fHistEventMult->GetXaxis()->SetBinLabel(13,"NTrigger>1 (MC)");
  fHistEventMult->GetXaxis()->SetBinLabel(14,"NFirstPartMC < NFirstReco");
  fHistEventMult->GetXaxis()->SetBinLabel(15,"NFirstPartMC=0 && NFirstPart=1");
  fHistEventMult->GetXaxis()->SetBinLabel(16,"NFirstPartMC!=0 && NFirstPart!=0");
  fHistEventMult->GetXaxis()->SetBinLabel(17,"NSecondPartMC==0 && NSecondPart!=0 (NT>0)");
  fHistEventMult->GetXaxis()->SetBinLabel(18,"NSecondPartMC<NSecondRecoTrue (NT>0)");
  fHistEventMult->GetXaxis()->SetBinLabel(19,"NTrigger>0 && DoubleCounted");
  fHistEventMult->GetXaxis()->SetBinLabel(20,"Common daughterPos");
  fHistEventMult->GetXaxis()->SetBinLabel(21,"Common daughterNeg");
  fHistEventMult->GetXaxis()->SetBinLabel(22,"SelEv (ACEvents)"); //all the events used for angular correlation
  fHistEventMult->GetXaxis()->SetBinLabel(23,"All events");
  fHistEventMult->GetXaxis()->SetBinLabel(24,"AOD event");

  fHistEventV0=new TH1F("fHistEventV0", "fHistEventV0",28, 0.5, 28.5);
  fHistEventV0->SetTitle("Number of V0 which progressively pass the listed selections");
  fHistEventV0->GetXaxis()->SetBinLabel(1,"All V0s");
  fHistEventV0->GetXaxis()->SetBinLabel(2,"V0s OnFly");
  fHistEventV0->GetXaxis()->SetBinLabel(3,"Filterbit daughters"); 
  fHistEventV0->GetXaxis()->SetBinLabel(4,"Chis daughter tracks"); 
  fHistEventV0->GetXaxis()->SetBinLabel(5,"TPC Refit"); 
  fHistEventV0->GetXaxis()->SetBinLabel(6,"NClusters>50"); 
  fHistEventV0->GetXaxis()->SetBinLabel(7,"Reject kink daughter"); 
  fHistEventV0->GetXaxis()->SetBinLabel(8,"Length DTracks>90cm"); 
  fHistEventV0->GetXaxis()->SetBinLabel(9,"CrossedRows/Length>0.8"); 
  fHistEventV0->GetXaxis()->SetBinLabel(10,"TOF or SPD hit for pileup"); 
  fHistEventV0->GetXaxis()->SetBinLabel(11,"PID daughters"); 
  fHistEventV0->GetXaxis()->SetBinLabel(12,"|eta daughters|<0.8"); 
  fHistEventV0->GetXaxis()->SetBinLabel(13,"NCrossedRows>70 && Crossed/Find>0.8"); 
  fHistEventV0->GetXaxis()->SetBinLabel(14,"|eta_K0s|<0.8"); 
  fHistEventV0->GetXaxis()->SetBinLabel(15,"more cuts + 0.45 < Mass < 0.55"); 
  fHistEventV0->GetXaxis()->SetBinLabel(16,"Lrejection"); 
  fHistEventV0->GetXaxis()->SetBinLabel(17,"DCAxy and DCAz Daughter cuts (not applied)"); 
  fHistEventV0->GetXaxis()->SetBinLabel(18,"FilterBit(1) on daughters (not applied)"); 
  fHistEventV0->GetXaxis()->SetBinLabel(19,"0<pT<pTV0max (reco up to here)"); 
  fHistEventV0->GetXaxis()->SetBinLabel(20,"NV0(reco) in ev wNT>0"); 
  fHistEventV0->GetXaxis()->SetBinLabel(21,"NV0(reco true) in ev wNT>0");
  fHistEventV0->GetXaxis()->SetBinLabel(22,"NV0(MC) in ev wNT>0");
  fHistEventV0->GetXaxis()->SetBinLabel(23,"N V0 pT>pTTrig(reco) in ev wNT>0");
  fHistEventV0->GetXaxis()->SetBinLabel(24,"N V0 pT>pTTrig(MC) in ev wNT>0");
  fHistEventV0->GetXaxis()->SetBinLabel(25,"NV0(reco) in SelEv");
  fHistEventV0->GetXaxis()->SetBinLabel(26,"NV0(reco true) in SelEv");
  fHistEventV0->GetXaxis()->SetBinLabel(27,"NV0(MC) in SelEv"); 
  fHistEventV0->GetXaxis()->SetBinLabel(28,"!OOB (Balbino)"); 

  fHistTrack=new TH1F("fHistTrack", "fHistTrack", 19, 0.5, 19.5);
  fHistTrack->GetXaxis()->SetBinLabel(1,"All tracks");
  fHistTrack->GetXaxis()->SetBinLabel(2,"Tracks after filterbit");
  fHistTrack->GetXaxis()->SetBinLabel(3,"Tracks with |eta| < 0.8"); 
  fHistTrack->GetXaxis()->SetBinLabel(4,"Track quality");
  fHistTrack->GetXaxis()->SetBinLabel(5,"TPCCrossedRows>80");
  fHistTrack->GetXaxis()->SetBinLabel(6, " Crossed rows/findable >0.8");
  fHistTrack->GetXaxis()->SetBinLabel(7,"TrackLength>90");
  fHistTrack->GetXaxis()->SetBinLabel(8, "Crossed rows/track length> 0.8");
  fHistTrack->GetXaxis()->SetBinLabel(9,"Charged tracks");
  fHistTrack->GetXaxis()->SetBinLabel(10,"DCAxy < 0.010+0.035/pt**1.1");
  fHistTrack->GetXaxis()->SetBinLabel(11,"DCAz <2");
  fHistTrack->GetXaxis()->SetBinLabel(12,"N.trigger"); //NumberPrimary in all selected events 
  fHistTrack->GetXaxis()->SetBinLabel(13,"N.trigger MC");
  fHistTrack->GetXaxis()->SetBinLabel(14,"N.trigger (NT>0)"); //NumberPrimary in all selected events 
  fHistTrack->GetXaxis()->SetBinLabel(15,"N.trigger MC (NT>0)");
  fHistTrack->GetXaxis()->SetBinLabel(16,"N.trigger (NV0>0)"); //NumberPrimary in events with at least one V0 (one reco V0 for data, one true V0 for MC)
  fHistTrack->GetXaxis()->SetBinLabel(17,"N.trigger MC (NV0>0)");

  fHistTOFBunchCrossing = new TH1F("fHistTOFBunchCrossing", "fHistTOFBunchCrossing", 200, -100, 100);

  fHistTrackAssoc=new TH1F("fHistTrackAssoc", "fHistTrackAssoc", 20, 0.5, 20.5);
  fHistTrackAssoc->GetXaxis()->SetBinLabel(1,"All tracks");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(2,"TrackAssocs after filterbit");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(3,"TrackAssocs with |eta| < 0.8"); 
  fHistTrackAssoc->GetXaxis()->SetBinLabel(4,"TrackAssoc quality");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(5,"TPCCrossedRows>80");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(6,"Crossed rows/findable >0.8");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(7,"Track Length > 90");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(8,"Crossed rows/track length >0.8");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(9,"Charged tracks");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(10,"DCAxy < 0.010+0.035/pt**1.1");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(11,"DCAz <2");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(12,"0<pT<pTV0max (reco up to here)"); 
  fHistTrackAssoc->GetXaxis()->SetBinLabel(13,"NAssoc(reco) in ev wNT>0"); 
  fHistTrackAssoc->GetXaxis()->SetBinLabel(14,"NAssoc(reco true) in ev wNT>0");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(15,"NAssoc(MC) in ev wNT>0");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(16,"pT > pTTrigg (reco) in ev wNT>0 ");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(17,"pT > pTTrigg (MC) in ev wNT>0");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(18,"NAssoc(reco) in SelEv");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(19,"NAssoc(reco true) in SelEv");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(20,"NAssoc(MC) in SelEv"); 

  fHistLengthvsCrossedRows = new TH2F("fHistLengthvsCrossedRows", "fHistLengthvsCrossedRows",  160, 0, 160, 160, 0, 160);
  fHistLengthvsCrossedRows->GetXaxis()->SetTitle("Number of Crossed rows");
  fHistLengthvsCrossedRows->GetYaxis()->SetTitle("Track length");

  fHistLengthvsCrossedRowsAfterSel = new TH2F("fHistLengthvsCrossedRowsAfterSel", "fHistLengthvsCrossedRowsAfterSel",  160, 0, 160, 160, 0, 160);
  fHistLengthvsCrossedRowsAfterSel->GetXaxis()->SetTitle("Number of Crossed rows");
  fHistLengthvsCrossedRowsAfterSel->GetYaxis()->SetTitle("Track length");

  fHistIsCommonParton = new TH1F("fHistIsCommonParton", "fHistIsCommonParton", 2, -0.5,1.5);

  fHistCommonParton = new TH3F("fHistCommonParton", "fHistCommonParton", 50, -25,25, 10, 0, 10, 10, 0, 10);
  fHistCommonParton->GetXaxis()->SetTitle("PdgCode");
  fHistCommonParton->GetYaxis()->SetTitle("TriggerLevel");
  fHistCommonParton->GetZaxis()->SetTitle("K0sLevel");

  fHistCommonPartonTrueCasc = new TH3F("fHistCommonPartonTrueCasc", "fHistCommonPartonTrueCasc", 50, -25,25, 10, 0, 10, 10, 0, 10);
  fHistCommonPartonTrueCasc->GetXaxis()->SetTitle("PdgCode");
  fHistCommonPartonTrueCasc->GetYaxis()->SetTitle("TriggerLevel");
  fHistCommonPartonTrueCasc->GetZaxis()->SetTitle("K0sLevel");

  fHistTriggerComposition=new TH3F("fHistTriggerComposition", "fHistTriggerComposition",10000 , -5000, 5000, 2, 0,2, 300, 0, 30);
  fHistTriggerComposition->GetYaxis()->SetTitle("0=NotPrim, 1=Primary");
  fHistTriggerComposition->GetZaxis()->SetTitle("p_{T} Trigger (GeV/c)");

  fHistTriggerCompositionMCTruth=new TH3F("fHistTriggerCompositionMCTruth", "fHistTriggerCompositionMCTruth",10000 , -5000, 5000, 2, 0,2, 300, 0, 30);
  fHistTriggerCompositionMCTruth->GetYaxis()->SetTitle("0=NotPrim, 1=Primary");
  fHistTriggerCompositionMCTruth->GetZaxis()->SetTitle("p_{T} Trigger (GeV/c)");

  fHistAssocComposition=new TH3F("fHistAssocComposition", "fHistAssocComposition",10000 , -5000, 5000, 2, 0,2, 300, 0, 30);
  fHistAssocComposition->GetYaxis()->SetTitle("0=NotPrim, 1=Primary");
  fHistAssocComposition->GetZaxis()->SetTitle("p_{T} Trigger (GeV/c)");

  fHistAssocCompositionMCTruth=new TH3F("fHistAssocCompositionMCTruth", "fHistAssocCompositionMCTruth",10000 , -5000, 5000, 2, 0,2, 300, 0, 30);
  fHistAssocCompositionMCTruth->GetYaxis()->SetTitle("0=NotPrim, 1=Primary");
  fHistAssocCompositionMCTruth->GetZaxis()->SetTitle("p_{T} Trigger (GeV/c)");

  fMassV0= new TH1F("fMassV0", "Invariant mass of V0 candidates", 100, 0.45, 0.55);
  fMassV0->GetXaxis()->SetTitle("M_{#pi^+ #pi^-}");

  fDCAxyDaughters=new TH2F("fDCAxyDaughters", "DCAxy Neg Daughter vs DCAxy Pos Daughter of V0 candidates", 100, -10, 10, 100, -10, 10);
  fDCAxyDaughters->GetXaxis()->SetTitle("DCAxyPos");
  fDCAxyDaughters->GetYaxis()->SetTitle("DCAxyNeg");
  fDCAzDaughters =new TH2F("fDCAzDaughters", "DCAxy Neg Daughter vs DCAxy Pos Daughter of V0 candidates", 100, -10, 10, 100, -10, 10);
  fDCAzDaughters->GetXaxis()->SetTitle("DCAzPos");
  fDCAzDaughters->GetYaxis()->SetTitle("DCAzNeg");

  fDCAxyDaughtersBis=new TH2F("fDCAxyDaughtersBis", "DCAxy Neg Daughter vs DCAxy Pos Daughter of V0 candidates", 100, -10, 10, 100, -10, 10);
  fDCAxyDaughtersBis->GetXaxis()->SetTitle("DCAxyPos");
  fDCAxyDaughtersBis->GetYaxis()->SetTitle("DCAxyNeg");
  fDCAzDaughtersBis =new TH2F("fDCAzDaughtersBis", "DCAxy Neg Daughter vs DCAxy Pos Daughter of V0 candidates", 100, -10, 10, 100, -10, 10);
  fDCAzDaughtersBis->GetXaxis()->SetTitle("DCAzPos");
  fDCAzDaughtersBis->GetYaxis()->SetTitle("DCAzNeg");

  fDCAxyPos=new TH1F("fDCAxyPos", "DCAxy Pos Daughter of V0 candidates", 100, -5, 5);
  fDCAxyPos->GetXaxis()->SetTitle("DCAxyPos");
  fDCAxyNeg=new TH1F("fDCAxyNeg", "DCAxy Neg Daughter of V0 candidates", 100, -5, 5);
  fDCAxyNeg->GetXaxis()->SetTitle("DCAxyNeg");

  fDCAzPos =new TH1F("fDCAzPos", "DCAz Pos Daughter of V0 candidates", 100, -5, 5);
  fDCAzPos->GetXaxis()->SetTitle("DCAzPos");
  fDCAzNeg =new TH1F("fDCAzNeg", "DCAz Neg Daughter of V0 candidates", 100, -5, 5);
  fDCAzNeg->GetXaxis()->SetTitle("DCAzNeg");

  fHistSecondParticleAll= new TH2F("fHistSecondParticleAll", "Number of V0 MCTrue vs number V0 reco (T>0) ", 60,-0.5,59.5,60,-0.5,59.5);
  fHistSecondParticleAll->GetXaxis()->SetTitle("Number (reco)");
  fHistSecondParticleAll->GetYaxis()->SetTitle("Number (MC)");

  fHistSecondParticleTruthAll= new TH2F("fHistSecondParticleTruthAll", "Number of V0 MCTrue vs number V0 reco (true) (T>0)", 60,-0.5,59.5,60,-0.5,59.5);
  fHistSecondParticleTruthAll->GetXaxis()->SetTitle("Number (reco true)");
  fHistSecondParticleTruthAll->GetYaxis()->SetTitle("Number (MC)");

  fHistSecondParticle= new TH2F("fHistSecondParticle", "Number of V0 MCTrue vs number V0 reco (T>0, V>0)", 60,-0.5,59.5,60,-0.5,59.5);
  fHistSecondParticle->GetXaxis()->SetTitle("Number (reco)");
  fHistSecondParticle->GetYaxis()->SetTitle("Number (MC)");

  fHistSecondParticleTruth= new TH2F("fHistSecondParticleTruth", "Number of V0 MCTrue vs number V0 reco (true) (T>0, V>0)", 60,-0.5,59.5,60,-0.5,59.5);
  fHistSecondParticleTruth->GetXaxis()->SetTitle("Number (reco true)");
  fHistSecondParticleTruth->GetYaxis()->SetTitle("Number (MC)");

  fHistMassvsPt_tagli = new TH2F *[6];

  for(Int_t j=0; j<6; j++){
    fHistMassvsPt_tagli[j] = new TH2F(Form("fHistMassvsPt_" +fV0+ "_tagli_%i",j),Form("fHistMassvsPt_" +fV0+ "_tagli_%i" + " (all selections on V0 applied)",j),400,0.3,0.7,160,0,16);
    fHistMassvsPt_tagli[j]->GetXaxis()->SetTitle("Invariant mass of V0 candidate");
    fHistMassvsPt_tagli[j]->GetYaxis()->SetTitle("p_{T} of V0 candidate");   
  }
  
  fHistMultvsTrigger=new TH2F("fHistMultvsTrigger", "Centrality of selected events (T>0, V>0) vs number of trigger particles", 30, -0.5, 29.5, 100, 0, 100);
  fHistMultvsTrigger->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTrigger->GetYaxis()->SetTitle("Centrality");
  

  fHistMultvsTriggerMCTruth=new TH2F("fHistMultvsTriggerMCTruth", "Centrality of selected events (T>0, V>0) vs number of trigger particles, MC Truth", 30, -0.5, 29.5, 100, 0, 100);
  fHistMultvsTriggerMCTruth->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerMCTruth->GetYaxis()->SetTitle("Centrality");
  

  fHistMultvsTriggerBefAll=new TH2F("fHistMultvsTriggerBefAll", "Centrality of events vs number of trigger particles", 30, -0.5, 29.5, 100, 0, 100);
  fHistMultvsTriggerBefAll->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerBefAll->GetYaxis()->SetTitle("Centrality");
  
 
  fHistMultvsTriggerMCTruthBefAll=new TH2F("fHistMultvsTriggerMCTruthBefAll", "Centrality of events vs number of trigger particles, MC Truth", 30, -0.5, 29.5, 100, 0, 100);
  fHistMultvsTriggerMCTruthBefAll->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerMCTruthBefAll->GetYaxis()->SetTitle("Centrality");
  
  fHistMultvsV0All=new TH2F("fHistMultvsV0All", "Centrality of events w T>0 vs number of reco V0s",60, -0.5, 59.5,100, 0, 100 );
  fHistMultvsV0All->GetXaxis()->SetTitle("Number of V0 particles");
  fHistMultvsV0All->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0AllTruth=new TH2F("fHistMultvsV0AllTruth", "Centrality of events w T>0 vs number of reco true V0s",60, -0.5, 59.5,100, 0, 100 );
  fHistMultvsV0AllTruth->GetXaxis()->SetTitle("Number of V0 reco particles");
  fHistMultvsV0AllTruth->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0MCAll=new TH2F("fHistMultvsV0MCAll", "Centrality of events w T>0 vs number of true V0s",60, -0.5, 59.5,100, 0, 100 );
  fHistMultvsV0MCAll->GetXaxis()->SetTitle("Number of V0 true particles");
  fHistMultvsV0MCAll->GetYaxis()->SetTitle("Centrality");

  fHistTriggerNotLeading=new TH3F("fHistTriggerNotLeading", "Events with trigger not leading in all events with NT>0",60, -0.5, 59.5,100, 0, 100, 60, 0, 30 );
  fHistTriggerNotLeading->GetXaxis()->SetTitle("Number of V0 with p_{T}> p_{T} Trigger");
  fHistTriggerNotLeading->GetYaxis()->SetTitle("Multiplicity class");
  fHistTriggerNotLeading->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");

  fHistTriggerNotLeadingMC=new TH3F("fHistTriggerNotLeadingMC", "Events with trigger not leading in all events with NT>0 (MC Truth)",60, -0.5, 59.5,100, 0, 100, 60, 0, 30 );
  fHistTriggerNotLeadingMC->GetXaxis()->SetTitle("Number of V0 with p_{T}> p_{T} Trigger");
  fHistTriggerNotLeadingMC->GetYaxis()->SetTitle("Multiplicity class");
  fHistTriggerNotLeadingMC->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");

  fHistMassPhoton=new TH1F("fHistMassPhoton", "Inv Mass of two V0 daughters as electrons, after V0 mass selection",  300, 0, 1.5);
  fHistMass2Photon=new TH1F("fHistMass2Photons", "Inv MassSquared of two daughters as electrons",  100, -3, 3);
  
  fHistPtArmvsAlpha=new TH2F("fHistPtArmvsAlpha", "Distribution of V0 candidates before cuts on Arm/Lrej/mass",  80, -1, 1,80, 0, 0.3);
  fHistPtArmvsAlpha->GetXaxis()->SetTitle("Alpha");
  fHistPtArmvsAlpha->GetYaxis()->SetTitle("Pt Armenteros");

  fHistPtArmvsAlphaAfterSelection=new TH2F("fHistPtArmvsAlphaAfterSelection", "Distribution of V0 candidates after applying mass cuts", 80, -1, 1,80, 0, 0.3);

  fHistPtArmvsAlphaAfterPhotonSelection=new TH2F("fHistPtArmvsAlphaAfterPhotonSelection", "Distribution of V0 candidates after photon cuts",  80, -1, 1,80, 0, 0.3);
  fHistPtArmvsAlphaAfterPhotonSelection->GetXaxis()->SetTitle("Alpha");
  fHistPtArmvsAlphaAfterPhotonSelection->GetYaxis()->SetTitle("Pt Armenteros");

  fHistPtArmvsAlphaAfterLambdaRejectionSelection=new TH2F("fHistPtArmvsAlphaAfterLambdaRejectionSelection", "Distribution of V0 candidates after Lrejection",  80, -1, 1,80, 0, 0.3);
  fHistPtArmvsAlphaAfterLambdaRejectionSelection->GetXaxis()->SetTitle("Alpha");
  fHistPtArmvsAlphaAfterLambdaRejectionSelection->GetYaxis()->SetTitle("Pt Armenteros");

  fHistTrigger=new TH1F("fHistTrigger", "Number of reco trigger particle distribution for selected events (also T=0)", 60, -0.5, 59.5); // each entry is an event

  fHistTriggerwV0=new TH1F("fHistTriggerwV0", "Number of reco trigger particle distribution for events used for AC", 60, -0.5, 59.5); // each entry is an event

  fHistTriggerMCTruth=new TH1F("fHistTriggerMCTruth", "Number of true trigger particle distribution for selected events (also T=0)", 60, -0.5, 59.5); // each entry is an event

  fHistTriggerwV0MCTruth=new TH1F("fHistTriggerwV0MCTruth", "Number of true trigger particle distribution for events used for AC", 60, -0.5, 59.5); // each entry is an event

  fHistMultiplicityVsVertexZ=new TH2F("fHistMultiplicityVsVertexZ", "Centrality vs Z vertex of selected events with NT>0 and NV0>0 ",  20, -10, 10,100, 0, 100);
      
  fHistTriggervsMult=new TH1F("fHistTriggervsMult", "Numero di particelle di trigger nei vari intervalli di centralita'", 100, 0, 100);
  fHistTriggervsMult->GetXaxis()->SetTitle("Centrality");

  fHistTriggervsMultMC=new TH1F("fHistTriggervsMultMC", "Numero di particelle di trigger (MCtruth) nei vari intervalli di centralita'", 100, 0, 100);
  fHistTriggervsMultMC->GetXaxis()->SetTitle("Centrality");

  fHistGeneratedTriggerPtPhi=new TH3F("fHistGeneratedTriggerPtPhi", "p_{T} and #phi distribution of generated trigger particles (charged, primary)", 600, 0, 30, 400,0, 2*TMath::Pi(),  100, 0, 100 );
  fHistGeneratedTriggerPtPhi->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedTriggerPtPhi->GetYaxis()->SetTitle("#phi");

  fHistGeneratedTriggerPtEta=new TH3F("fHistGeneratedTriggerPtEta", "p_{T} and #eta distribution of generated trigger particles (primary, charged)", 600, 0, 30, 400,-1.2,1.2,  100, 0, 100 );
  fHistGeneratedTriggerPtEta->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedTriggerPtEta->GetYaxis()->SetTitle("#eta");

  fHistSelectedTriggerPtPhi= new TH3F*[1];
  for(Int_t j=0; j<1; j++){
    fHistSelectedTriggerPtPhi[j]=new TH3F(Form("fHistSelectedTriggerPtPhi_%i",j), "p_{T} and #phi distribution of selected trigger particles (primary)", 600, 0, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
    fHistSelectedTriggerPtPhi[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedTriggerPtPhi[j]->GetYaxis()->SetTitle("#phi");
  }
  fHistSelectedGenTriggerPtPhi= new TH3F*[3];
  for(Int_t j=0; j<3; j++){
    fHistSelectedGenTriggerPtPhi[j]=new TH3F(Form("fHistSelectedGenTriggerPtPhi_%i",j), "p_{T} and #phi distribution of selected trigger particles (primary) (p_{T} and #phi generated))", 600, 0, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
    fHistSelectedGenTriggerPtPhi[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedGenTriggerPtPhi[j]->GetYaxis()->SetTitle("#phi");
  }

  fHistSelectedTriggerPtEta= new TH3F*[1];
  for(Int_t j=0; j<1; j++){
    fHistSelectedTriggerPtEta[j]=new TH3F(Form("fHistSelectedTriggerPtEta_%i",j), "p_{T} and #eta distribution of selected trigger particles (primary)", 600, 0, 30, 400,-1.2, 1.2,  100, 0, 100);
    fHistSelectedTriggerPtEta[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedTriggerPtEta[j]->GetYaxis()->SetTitle("#eta");
  }
  fHistSelectedGenTriggerPtEta= new TH3F*[3];
  for(Int_t j=0; j<3; j++){
    fHistSelectedGenTriggerPtEta[j]=new TH3F(Form("fHistSelectedGenTriggerPtEta_%i",j), "p_{T} and #eta distribution of selected trigger particles (primary) (p_{T} and #eta generated))", 600, 0, 30, 400,-1.2, 1.2,  100, 0, 100);
    fHistSelectedGenTriggerPtEta[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedGenTriggerPtEta[j]->GetYaxis()->SetTitle("#eta");
  }
  
  fHistGeneratedV0PtTMaxPhi=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistGeneratedV0PtTMaxPhi[j]=new TH3F(Form("fHistGeneratedV0PtTMaxPhi_%i",j), "p^{Trigg, Max}_{T} and #phi distribution of generated V0 particles (K0s, primary, events w T>0)", 120, -30, 30, 400,0, 2*TMath::Pi(),  100, 0, 100 );
    fHistGeneratedV0PtTMaxPhi[j]->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
    fHistGeneratedV0PtTMaxPhi[j]->GetYaxis()->SetTitle("#phi");
  }
  
  fHistCPGeneratedV0PtTMaxPhi=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistCPGeneratedV0PtTMaxPhi[j]=new TH3F(Form("fHistCPGeneratedV0PtTMaxPhi_%i",j), "p^{Trigg, Max}_{T} and #phi distribution of generated V0 particles (K0s, primary, events w T>0)", 120, -30, 30, 400,0, 2*TMath::Pi(),  100, 0, 100 );
    fHistCPGeneratedV0PtTMaxPhi[j]->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
    fHistCPGeneratedV0PtTMaxPhi[j]->GetYaxis()->SetTitle("#phi");
  }

  fHistSelectedV0PtTMaxPhi=new TH3F*[1];
  for(Int_t j=0; j<1; j++){
    fHistSelectedV0PtTMaxPhi[j]=new TH3F(Form("fHistSelectedV0PtTMaxPhi_%i",j), "p^{Trigg, Max}_{T} and #phi distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
    fHistSelectedV0PtTMaxPhi[j]->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
    fHistSelectedV0PtTMaxPhi[j]->GetYaxis()->SetTitle("#phi");
  }
  
  fHistGeneratedV0PtTMaxEta=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistGeneratedV0PtTMaxEta[j]=new TH3F(Form("fHistGeneratedV0PtTMaxEta_%i",j), "p^{Trigg, Max}_{T} and #eta distribution of generated V0 particles (K0s, primary, events w T>0)", 120,-30, 30, 400,-1.2,1.2,  100, 0, 100 );
    fHistGeneratedV0PtTMaxEta[j]->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
    fHistGeneratedV0PtTMaxEta[j]->GetYaxis()->SetTitle("#eta");
  }
  
  fHistCPGeneratedV0PtTMaxEta=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistCPGeneratedV0PtTMaxEta[j]=new TH3F(Form("fHistCPGeneratedV0PtTMaxEta_%i",j), "p^{Trigg, Max}_{T} and #eta distribution of generated V0 particles (K0s, primary, events w T>0)", 120, -30, 30, 400,-1.2,1.2,  100, 0, 100 );
    fHistCPGeneratedV0PtTMaxEta[j]->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
    fHistCPGeneratedV0PtTMaxEta[j]->GetYaxis()->SetTitle("#eta");
  }

  fHistSelectedV0PtTMaxEta=new TH3F*[1];
  for(Int_t j=0; j<1; j++){
    fHistSelectedV0PtTMaxEta[j]=new TH3F(Form("fHistSelectedV0PtTMaxEta_%i",j), "p^{Trigg, Max}_{T} and #eta distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,-1.2,1.2,  100, 0, 100 );
    fHistSelectedV0PtTMaxEta[j]->GetXaxis()->SetTitle("p^{Trigg, max}_{T}");
    fHistSelectedV0PtTMaxEta[j]->GetYaxis()->SetTitle("#eta");
  }

  fHistGeneratedV0PtPtTMax=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistGeneratedV0PtPtTMax[j]=new TH3F(Form("fHistGeneratedV0PtPtTMax_%i",j), "p_{T} and p^{Trigg, Max}_{T} distribution of generated V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 120, -30, 30,  100, 0, 100 );
    fHistGeneratedV0PtPtTMax[j]->GetXaxis()->SetTitle("p_{T}");
    fHistGeneratedV0PtPtTMax[j]->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");
  }

  fHistCPGeneratedV0PtPtTMax=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistCPGeneratedV0PtPtTMax[j]=new TH3F(Form("fHistCPGeneratedV0PtPtTMax_%i",j), "p_{T} and p^{Trigg, Max}_{T} distribution of generated V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 120, -30, 30,  100, 0, 100 );
    fHistCPGeneratedV0PtPtTMax[j]->GetXaxis()->SetTitle("p_{T}");
    fHistCPGeneratedV0PtPtTMax[j]->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");
  }
  
  fHistSelectedV0PtPtTMax=new TH3F*[1];
  for(Int_t j=0; j<1; j++){
    fHistSelectedV0PtPtTMax[j]=new TH3F(Form("fHistSelectedV0PtPtTMax_%i",j), "p_{T} and p^{Trigg, Max}_{T} distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 60, 0,30,  100, 0, 100 );
    fHistSelectedV0PtPtTMax[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedV0PtPtTMax[j]->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");
  }
  fHistSelectedGenV0PtPtTMax=new TH3F*[1];
  for(Int_t j=0; j<1; j++){
    fHistSelectedGenV0PtPtTMax[j]=new TH3F(Form("fHistSelectedGenV0PtPtTMax_%i",j), "p_{T} and p^{Trigg, Max}_{T} distribution of selected V0 particles (K0s, primary, events w T>0) (p_{T} generated)", 300, 0, 30, 60, 0,30,  100, 0, 100 );
    fHistSelectedGenV0PtPtTMax[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedGenV0PtPtTMax[j]->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");
  }

  fHistGeneratedV0PtEta=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistGeneratedV0PtEta[j]=new TH3F(Form("fHistGeneratedV0PtEta_%i",j), "p_{T} and #eta distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 480,-1.2,1.2,  100, 0, 100 );
    fHistGeneratedV0PtEta[j]->GetXaxis()->SetTitle("p_{T}");
    fHistGeneratedV0PtEta[j]->GetYaxis()->SetTitle("#eta");
  }

  fHistSelectedV0PtEta=new TH3F*[1];
  for(Int_t j=0; j<1; j++){
    fHistSelectedV0PtEta[j]=new TH3F(Form("fHistSelectedV0PtEta_%i",j), "p_{T} and #eta distribution of selected V0 particles (K0s, primary, events w T>0)", 60, 0, 30, 480,-1.2,1.2,  100, 0, 100 );
    fHistSelectedV0PtEta[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedV0PtEta[j]->GetYaxis()->SetTitle("#eta");
  }

  fHistReconstructedV0PtMass=new TH3F("fHistReconstructedV0PtMass", "p_{T} and mass distribution of reconstructed V0 particles(K0s, primary, event w T>0)", 100, 0.45, 0.55, 160, 0, 16,  100, 0, 100);
  fHistReconstructedV0PtMass->GetYaxis()->SetTitle("p_{T}");
  fHistReconstructedV0PtMass->GetXaxis()->SetTitle("M_{pi^{+} #pi^{-}}");

  fHistSelectedV0PtMass=new TH3F("fHistSelectedV0PtMass", "p_{T} and mass distribution of selected V0 particles (K0s, primary, event w T>0)", 100, 0.45, 0.55,  160, 0, 16, 100, 0, 100);
  fHistSelectedV0PtMass->GetYaxis()->SetTitle("p_{T}");
  fHistSelectedV0PtMass->GetXaxis()->SetTitle("M_{pi^{+} #pi^{-}}");

  fHistTriggerPtRecovsPtGen = new TH2F ("fHistTriggerPtRecovsPtGen", "p_{T, reco} vs p_{T, gen} for primary trigger particles", 600, 0,30,600, 0,30);
  fHistTriggerPtRecovsPtGen->GetXaxis()->SetTitle("p_{T, gen}");
  fHistTriggerPtRecovsPtGen->GetYaxis()->SetTitle("p_{T, reco}");

  fHistTriggerPtRecovsPtGenNotPrim = new TH2F ("fHistTriggerPtRecovsPtGenNotPrim", "p_{T, reco} vs p_{T, gen} for non-primary trigger particles", 600, 0,30,600, 0,30);
  fHistTriggerPtRecovsPtGenNotPrim->GetXaxis()->SetTitle("p_{T, gen}");
  fHistTriggerPtRecovsPtGenNotPrim->GetYaxis()->SetTitle("p_{T, reco}");

  fHistTriggerPtRecovsPtGenPion = new TH2F ("fHistTriggerPtRecovsPtGenPion", "p_{T, reco} vs p_{T, gen} for primary Pion trigger particles", 600, 0,30,600, 0,30);
  fHistTriggerPtRecovsPtGenPion->GetXaxis()->SetTitle("p_{T, gen}");
  fHistTriggerPtRecovsPtGenPion->GetYaxis()->SetTitle("p_{T, reco}");

  fHistTriggerPtRecovsPtGenProton = new TH2F ("fHistTriggerPtRecovsPtGenProton", "p_{T, reco} vs p_{T, gen} for primary Proton trigger particles", 600, 0,30,600, 0,30);
  fHistTriggerPtRecovsPtGenProton->GetXaxis()->SetTitle("p_{T, gen}");
  fHistTriggerPtRecovsPtGenProton->GetYaxis()->SetTitle("p_{T, reco}");

  fHistTriggerPtRecovsPtGenKaon = new TH2F ("fHistTriggerPtRecovsPtGenKaon", "p_{T, reco} vs p_{T, gen} for primary Kaon trigger particles", 600, 0,30,600, 0,30);
  fHistTriggerPtRecovsPtGenKaon->GetXaxis()->SetTitle("p_{T, gen}");
  fHistTriggerPtRecovsPtGenKaon->GetYaxis()->SetTitle("p_{T, reco}");

  fHistAssocPtRecovsPtGen = new TH2F ("fHistAssocPtRecovsPtGen", "p_{T, reco} vs p_{T, gen} for primary associated particles", 600, 0,30,600, 0,30);
  fHistAssocPtRecovsPtGen->GetXaxis()->SetTitle("p_{T, gen}");
  fHistAssocPtRecovsPtGen->GetYaxis()->SetTitle("p_{T, reco}");

  fHistAssocPtRecovsPtGenNotPrim = new TH2F ("fHistAssocPtRecovsPtGenNotPrim", "p_{T, reco} vs p_{T, gen} for non-primary associated particles", 600, 0,30,600, 0,30);
  fHistAssocPtRecovsPtGenNotPrim->GetXaxis()->SetTitle("p_{T, gen}");
  fHistAssocPtRecovsPtGenNotPrim->GetYaxis()->SetTitle("p_{T, reco}");

  fHistAssocPtRecovsPtGenPion = new TH2F ("fHistAssocPtRecovsPtGenPion", "p_{T, reco} vs p_{T, gen} for primary Pion associated particles", 600, 0,30,600, 0,30);
  fHistAssocPtRecovsPtGenPion->GetXaxis()->SetTitle("p_{T, gen}");
  fHistAssocPtRecovsPtGenPion->GetYaxis()->SetTitle("p_{T, reco}");

  fHistAssocPtRecovsPtGenProton = new TH2F ("fHistAssocPtRecovsPtGenProton", "p_{T, reco} vs p_{T, gen} for primary Proton associated particles", 600, 0,30,600, 0,30);
  fHistAssocPtRecovsPtGenProton->GetXaxis()->SetTitle("p_{T, gen}");
  fHistAssocPtRecovsPtGenProton->GetYaxis()->SetTitle("p_{T, reco}");

  fHistAssocPtRecovsPtGenKaon = new TH2F ("fHistAssocPtRecovsPtGenKaon", "p_{T, reco} vs p_{T, gen} for primary Kaon associated particles", 600, 0,30,600, 0,30);
  fHistAssocPtRecovsPtGenKaon->GetXaxis()->SetTitle("p_{T, gen}");
  fHistAssocPtRecovsPtGenKaon->GetYaxis()->SetTitle("p_{T, reco}");

  fHistResolutionTriggerPt=new TH2F("fHistResolutionTriggerPt", "p_{T} resolution of selected trigger particles (primary)", 2000, -10, 10,  60,0,30);
  fHistResolutionTriggerPt->GetXaxis()->SetTitle("p_{T, DATA}-p_{T, MC}");
  fHistResolutionTriggerPt->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionTriggerPhi=new TH2F("fHistResolutionTriggerPhi", "#Phi resolution of selected trigger particles (primary)", 2000, -1, 1,  60,0,30);
  fHistResolutionTriggerPhi->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionTriggerPhi->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionTriggerEta=new TH2F("fHistResolutionTriggerEta", "#Eta resolution of selected trigger particles (primary)", 2000, -1, 1,  60,0,30);
  fHistResolutionTriggerEta->GetXaxis()->SetTitle("Eta_{DATA}-Eta_{MC}");
  fHistResolutionTriggerEta->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionTriggerPhiPhi=new TH2F("fHistResolutionTriggerPhiPhi", "#Phi resolution of selected trigger particles (primary)", 2000, -1, 1,  100, 0, 2*TMath::Pi());
  fHistResolutionTriggerPhiPhi->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionTriggerPhiPhi->GetYaxis()->SetTitle("#phi^{Trigg, max}");

  fHistResolutionTriggerPhiEta=new TH2F("fHistResolutionTriggerPhiEta", "#Phi resolution of selected trigger particles (primary)", 2000, -1, 1,  100, -1,1);
  fHistResolutionTriggerPhiEta->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionTriggerPhiEta->GetYaxis()->SetTitle("#eta^{Trigg, max}");

  fHistResolutionV0Pt=new TH2F("fHistResolutionV0Pt", "p_{T} resolution of selected V0 particles (K0s, primary, event w T>0)", 2000, -10, 10,  60,0,30);
  fHistResolutionV0Pt->GetXaxis()->SetTitle("p_{T, DATA}-p_{T, MC}");
  fHistResolutionV0Pt->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionV0Phi=new TH2F("fHistResolutionV0Phi", "#Phi resolution of selected V0 particles (K0s, primary, event w T>0)", 2000, -1, 1,  60,0,30);
  fHistResolutionV0Phi->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionV0Phi->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionV0PhivsPt=new TH2F("fHistResolutionV0PhivsPt", "#Phi resolution of selected V0 particles (K0s, primary, event w T>0)", 2000, -1, 1,  60,0,30);
  fHistResolutionV0PhivsPt->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionV0PhivsPt->GetYaxis()->SetTitle("p_{T}");

  fHistResolutionV0PtvsPt=new TH2F("fHistResolutionV0PtvsPt", "#Pt resolution of selected V0 particles (K0s, primary, event w T>0) vs p_{T} trigger", 2000, -10, 10, 60,0,30);
  fHistResolutionV0PtvsPt->GetXaxis()->SetTitle("p_{T,DATA}-p_{T,MC}");
  fHistResolutionV0PtvsPt->GetYaxis()->SetTitle("p_{T}");
  
  fHistResolutionV0Eta=new TH2F("fHistResolutionV0Eta", "#Eta resolution of selected V0 particles (K0s, primary, event w T>0)", 2000, -1, 1, 60,0,30);
  fHistResolutionV0Eta->GetXaxis()->SetTitle("Eta_{DATA}-Eta_{MC}");
  fHistResolutionV0Eta->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");

  /*
    fHistResolutionTriggerPt=new TH3F("fHistResolutionTriggerPt", "p_{T} resolution of selected trigger particles (primary)", 2000, -5, 5, 100, 0, 100, 60,0,30);
    fHistResolutionTriggerPt->GetXaxis()->SetTitle("p_{T, DATA}-p_{T, MC}");
    fHistResolutionTriggerPt->GetYaxis()->SetTitle("Centrality");
    fHistResolutionTriggerPt->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");

    fHistResolutionTriggerPhi=new TH3F("fHistResolutionTriggerPhi", "#Phi resolution of selected trigger particles (primary)", 2000, -1, 1, 100, 0, 100, 60,0,30);
    fHistResolutionTriggerPhi->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
    fHistResolutionTriggerPhi->GetYaxis()->SetTitle("Centrality");
    fHistResolutionTriggerPhi->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");

    fHistResolutionTriggerEta=new TH3F("fHistResolutionTriggerEta", "#Eta resolution of selected trigger particles (primary)", 2000, -1, 1, 100, 0, 100, 60,0,30);
    fHistResolutionTriggerEta->GetXaxis()->SetTitle("Eta_{DATA}-Eta_{MC}");
    fHistResolutionTriggerEta->GetYaxis()->SetTitle("Centrality");
    fHistResolutionTriggerEta->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");

    fHistResolutionV0Pt=new TH3F("fHistResolutionV0Pt", "p_{T} resolution of selected V0 particles (K0s, primary, event w T>0)", 2000, -5, 5, 100, 0, 100, 60,0,30);
    fHistResolutionV0Pt->GetXaxis()->SetTitle("p_{T, DATA}-p_{T, MC}");
    fHistResolutionV0Pt->GetYaxis()->SetTitle("Centrality");
    fHistResolutionV0Pt->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");

    fHistResolutionV0Phi=new TH3F("fHistResolutionV0Phi", "#Phi resolution of selected V0 particles (K0s, primary, event w T>0)", 2000, -1, 1, 100, 0, 100, 60,0,30);
    fHistResolutionV0Phi->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
    fHistResolutionV0Phi->GetYaxis()->SetTitle("Centrality");
    fHistResolutionV0Phi->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");

    fHistResolutionV0PhivsPt=new TH3F("fHistResolutionV0PhivsPt", "#Phi resolution of selected V0 particles (K0s, primary, event w T>0)", 4000, -10, 10, 100, 0, 100, 60,0,30);
    fHistResolutionV0PhivsPt->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
    fHistResolutionV0PhivsPt->GetYaxis()->SetTitle("Centrality");
    fHistResolutionV0PhivsPt->GetZaxis()->SetTitle("p_{T}");

    fHistResolutionV0PtvsPt=new TH3F("fHistResolutionV0PtvsPt", "#Pt resolution of selected V0 particles (K0s, primary, event w T>0) vs p_{T} trigger", 2000, -5, 5, 100, 0, 100, 60,0,30);
    fHistResolutionV0PtvsPt->GetXaxis()->SetTitle("p_{T,DATA}-p_{T,MC}");
    fHistResolutionV0PtvsPt->GetYaxis()->SetTitle("Centrality");
    fHistResolutionV0PtvsPt->GetZaxis()->SetTitle("p_{T}");
  
    fHistResolutionV0Eta=new TH3F("fHistResolutionV0Eta", "#Eta resolution of selected V0 particles (K0s, primary, event w T>0)", 2000, -1, 1, 100, 0, 100, 60,0,30);
    fHistResolutionV0Eta->GetXaxis()->SetTitle("Eta_{DATA}-Eta_{MC}");
    fHistResolutionV0Eta->GetYaxis()->SetTitle("Centrality");
    fHistResolutionV0Eta->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");
  */

  fHistPrimaryTrigger= new TH2F **[6];
  for(Int_t j=0; j<6; j++){
    fHistPrimaryTrigger[j]=new TH2F*[1];
    for(Int_t i=0; i<1; i++){
      fHistPrimaryTrigger[j][i]=new TH2F(Form("fHistPrimaryTrigger_%i_cut%i", j,i), "Trigger MC (selected)", 4, 0.5, 4.5, 100, 0,30 );
      fHistPrimaryTrigger[j][i]->GetXaxis()->SetBinLabel(1,"Primary selected triggers");
      fHistPrimaryTrigger[j][i]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected triggers"); 
      fHistPrimaryTrigger[j][i]->GetXaxis()->SetBinLabel(3,"Secondary from material selected triggers"); 
      fHistPrimaryTrigger[j][i]->GetYaxis()->SetTitle("p_{T}"); 
    }
  }

  fHistPrimaryV0= new TH3F**[6];
  for(Int_t j=0; j<6; j++){
    fHistPrimaryV0[j]=new TH3F*[1];
    for(Int_t i=0; i<1; i++){
      fHistPrimaryV0[j][i]=new TH3F(Form("fHistPrimaryV0_%i_cut%i",j,i), "V0 MC (K0s, selected)", 4, 0.5, 4.5, 160, 0, 16, 60, 0, 30);
      fHistPrimaryV0[j][i]->GetXaxis()->SetBinLabel(1,"Primary selected V0s");
      fHistPrimaryV0[j][i]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected V0s"); 
      fHistPrimaryV0[j][i]->GetXaxis()->SetBinLabel(3,"Secondary from material selected V0s"); 
      fHistPrimaryV0[j][i]->GetYaxis()->SetTitle("p_{T}");
      fHistPrimaryV0[j][i]->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");
    } 
  }

  fHistMultiplicityOfMixedEvent=new TH2F("fHistMultiplicityOfMixedEvent", "Distribution of number of events used for the mixing", 100, 0.5, 100.5, 100, 0, 100);

  fEventCuts.AddQAplotsToList(fOutputList);
  
  fOutputList->Add(fHistPDG);
  fOutputList->Add(fHistTrackBufferOverflow);
  fOutputList->Add(fHistEventMult);
  fOutputList->Add(fHistEventV0);
  fOutputList->Add(fHistTrack); 
  fOutputList->Add(fHistTOFBunchCrossing);
  fOutputList->Add(fHistLengthvsCrossedRowsAfterSel);
  fOutputList->Add(fHistLengthvsCrossedRows);
  fOutputList->Add(fHistIsCommonParton);
  fOutputList->Add(fHistCommonParton);
  fOutputList->Add(fHistCommonPartonTrueCasc);
  fOutputList->Add(fHistTriggerComposition); 
  fOutputList->Add(fHistTriggerCompositionMCTruth); 
  fOutputList->Add(fHistAssocComposition); 
  fOutputList->Add(fHistAssocCompositionMCTruth); 
  fOutputList->Add(fHistTrackAssoc); 
  
  //istogrammi riempiti ad ogni evento selezionato (ossia utilizzato per AC) con una entry
  fOutputList->Add(fHistZvertex);
  fOutputList->Add(fHistNumberChargedAllEvents);
  fOutputList->Add(fHistNumberChargedTrigger);
  fOutputList->Add(fHistNumberChargedNoTrigger);
  fOutputList->Add(fHist_multiplicity); 
  fOutputList->Add(fHist_multiplicity_EvwTrigger);
  fOutputList->Add(fHistMultiplicityVsVertexZ);
  fOutputList->Add(fHistFractionSharedTPCClusters);
  fOutputList->Add(fHistMultvsTriggerBefAll);
  fOutputList->Add(fHistMultvsTriggerMCTruthBefAll);
  fOutputList->Add(fHistMultvsTrigger);
  fOutputList->Add(fHistMultvsTriggerMCTruth);
  fOutputList->Add(fHistMultvsV0All);
  fOutputList->Add(fHistMultvsV0AllTruth);
  fOutputList->Add(fHistMultvsV0MCAll);
  fOutputList->Add(fHistTriggerNotLeading);
  fOutputList->Add(fHistTriggerNotLeadingMC);
  fOutputList->Add(fHistSecondParticleAll);
  fOutputList->Add(fHistSecondParticleTruthAll);
  fOutputList->Add(fHistSecondParticle);  
  fOutputList->Add(fHistSecondParticleTruth);
  fOutputList->Add(fHistPt);   
  fOutputList->Add(fHistDCAxym1);     
  fOutputList->Add(fHistDCAxym2);    
  fOutputList->Add(fHistDCAzm1);      
  fOutputList->Add(fHistDCAzm2);     
  fOutputList->Add(fHistPtV0);     
  fOutputList->Add(fHistPthAssoc);     
  fOutputList->Add(fHistPtTMaxBefAllCfrDataMC);
  fOutputList->Add(fHistPtTMinBefAll);     
  fOutputList->Add(fHistPtTMinBefAllMC);     
  fOutputList->Add(fHistPtvsMult);       
  fOutputList->Add(fHistPtvsMultBefAll);       
  fOutputList->Add(fHistPtMaxvsMult);       
  fOutputList->Add(fHistPtMaxvsMultBefAll);       
  fOutputList->Add(fHistPtMaxvsMultBefAllReco);      
  fOutputList->Add(fHistPtMaxvsMultBefAllGen);       
  fOutputList->Add(fHist_eta_phi);
  fOutputList->Add(fHist_eta_phi_PtMax);
  fOutputList->Add(fHistTriggervsMult);
  fOutputList->Add(fHistTriggervsMultMC);
  fOutputList->Add(fHistTrigger);
  fOutputList->Add(fHistTriggerwV0);
  fOutputList->Add(fHistTriggerMCTruth);
  fOutputList->Add(fHistTriggerwV0MCTruth);

  fOutputList->Add(fMassV0);
  fOutputList->Add(fDCAxyDaughters);
  fOutputList->Add(fDCAzDaughters);
  fOutputList->Add(fDCAxyDaughtersBis);
  fOutputList->Add(fDCAzDaughtersBis);
  fOutputList->Add(fDCAxyPos);
  fOutputList->Add(fDCAxyNeg);
  fOutputList->Add(fDCAzPos);
  fOutputList->Add(fDCAzNeg);




  for(Int_t j=0; j < 6; j++){
    fOutputList->Add(fHistMassvsPt_tagli[j]);
    
    for(Int_t i=0; i<1; i++){  
      fOutputList3->Add(fHistPrimaryTrigger[j][i]); 
    }
   
    for(Int_t l=0; l<1; l++){
      fOutputList2->Add(fHistPrimaryV0[j][l]);      
    }
    
  } 
  
  fOutputList->Add(fHistMass2Photon);
  fOutputList->Add(fHistMassPhoton);
  fOutputList->Add(fHistPtArmvsAlpha);
  fOutputList->Add(fHistPtArmvsAlphaAfterSelection);  
  fOutputList->Add(fHistPtArmvsAlphaAfterPhotonSelection);  
  fOutputList->Add(fHistPtArmvsAlphaAfterLambdaRejectionSelection);  
  fOutputList->Add(fHistMultiplicityOfMixedEvent);
   
  fOutputList3->Add(fHistGeneratedTriggerPtPhi);
  fOutputList3->Add(fHistGeneratedTriggerPtEta);
  
  for(Int_t j=0; j < 1; j++){
    fOutputList3->Add(fHistSelectedTriggerPtPhi[j]);
    fOutputList3->Add(fHistSelectedTriggerPtEta[j]);
    fOutputList3->Add(fHistSelectedGenTriggerPtPhi[j]);
    fOutputList3->Add(fHistSelectedGenTriggerPtEta[j]);
  }
 
  for(Int_t j=0; j < 2; j++){
    fOutputList2->Add(fHistGeneratedV0PtTMaxPhi[j]); 
    fOutputList2->Add(fHistGeneratedV0PtTMaxEta[j]); 
    fOutputList2->Add(fHistGeneratedV0PtPtTMax[j]); 
    fOutputList2->Add(fHistGeneratedV0PtEta[j]);
    fOutputList2->Add(fHistCPGeneratedV0PtTMaxPhi[j]);
    fOutputList2->Add(fHistCPGeneratedV0PtTMaxEta[j]);
    fOutputList2->Add(fHistCPGeneratedV0PtPtTMax[j]);
  }
  for(Int_t j=0; j < 1; j++){
    fOutputList2->Add(fHistSelectedV0PtTMaxPhi[j]);
    fOutputList2->Add(fHistSelectedV0PtTMaxEta[j]);
    fOutputList2->Add(fHistSelectedV0PtPtTMax[j]);
    fOutputList2->Add(fHistSelectedV0PtEta[j]);
    fOutputList2->Add(fHistSelectedGenV0PtPtTMax[j]);
  }
 
  fOutputList->Add(fHistReconstructedV0PtMass);
  fOutputList->Add(fHistSelectedV0PtMass);
  fOutputList4->Add(fHistResolutionTriggerPt);
  fOutputList4->Add(fHistResolutionTriggerPhi);
  fOutputList4->Add(fHistResolutionTriggerEta);
  fOutputList4->Add(fHistResolutionTriggerPhiPhi);
  fOutputList4->Add(fHistResolutionTriggerPhiEta);
  fOutputList4->Add(fHistResolutionV0Pt);
  fOutputList4->Add(fHistResolutionV0Phi);
  fOutputList4->Add(fHistResolutionV0PhivsPt);
  fOutputList4->Add(fHistResolutionV0Eta);
  fOutputList4->Add(fHistResolutionV0PtvsPt);
  fOutputList4->Add(fHistTriggerPtRecovsPtGen);
  fOutputList4->Add(fHistAssocPtRecovsPtGen);
  fOutputList4->Add(fHistTriggerPtRecovsPtGenNotPrim);
  fOutputList4->Add(fHistAssocPtRecovsPtGenNotPrim);
  fOutputList4->Add(fHistTriggerPtRecovsPtGenPion);
  fOutputList4->Add(fHistAssocPtRecovsPtGenPion);
  fOutputList4->Add(fHistTriggerPtRecovsPtGenProton);
  fOutputList4->Add(fHistAssocPtRecovsPtGenProton);
  fOutputList4->Add(fHistTriggerPtRecovsPtGenKaon);
  fOutputList4->Add(fHistAssocPtRecovsPtGenKaon);

  PostData(1, fOutputList);  
  PostData(2, fSignalTree);       
  PostData(3, fBkgTree); 
  PostData(4, fOutputList2);     
  PostData(5, fOutputList3);     
  PostData(6, fOutputList4);     
     
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrelationhhK0s::UserExec(Option_t *)
{

  Float_t moltep[6]={0,5,10,30,50,100};  
  Float_t LastzBin;
  Float_t LastcentralityBin;

  fHistEventMult->Fill(23);   
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());  
  if(!fAOD) {
    AliWarning("Error: AOD event not available \n");
    PostData(1, fOutputList);       
    PostData(2, fSignalTree);       
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);     
    return;
  }        
  fHistEventMult->Fill(24);
  
  // fEventCuts.SetManualMode(true);
  // fEventCuts.SetupRun2pp();

  

  /// Use the event cut class to apply the required selections
  if (!fEventCuts.AcceptEvent(fAOD)) {   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3,fBkgTree);
    PostData(4, fOutputList2);     
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);     
    return;
  }


  fHistEventMult->Fill(1);
  
  Int_t iTracks(fAOD->GetNumberOfTracks());         
  Int_t V0Tracks(fAOD->GetNumberOfV0s());           
  Evcounter++;  
  // cout << "\n \n \n ********************************************************* "<< endl;
  // cout << "event number "<<  Evcounter << endl; 
  // cout << "number of tracks before any cut " << iTracks << endl;
  // cout << "number of V0 before any cut " << V0Tracks << endl;
  
  //VERTEX SELECTION AND TRIGGER
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
  const AliAODVertex *lPrimaryBestAODVtx = fAOD->GetPrimaryVertex();
  if (!lPrimaryBestAODVtx){
    AliWarning("No prim. vertex in AOD... return!");
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);  
    PostData(6, fOutputList4);     
    return;
  }
  fHistEventMult->Fill(2);


  AliVVertex *vertexmain =0x0;
  vertexmain = (AliVVertex*) lPrimaryBestAODVtx;
  lPrimaryBestAODVtx->GetXYZ(lBestPrimaryVtxPos);

  if (TMath::Abs(lBestPrimaryVtxPos[2])>10.){
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3, fBkgTree); 
    PostData(4, fOutputList2); 
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);     
    return;
  }
  fHistEventMult->Fill(3);

  //PID RESPONSE
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
  if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
  UInt_t mask = inputHandler->IsEventSelected();

  if (!fPIDResponse){
    AliWarning("cannot get pid response");
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2); 
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);     
    return;
  }
  fHistEventMult->Fill(4);

  Float_t lPercentiles = 0;
 
  
  //This will work for both ESDs and AODs
  AliMultSelection *MultSelection = (AliMultSelection*) fAOD -> FindListObject("MultSelection");
  
  if ( MultSelection ){
  //c cout << "mult sel ok" << endl;
  lPercentiles= MultSelection->GetMultiplicityPercentile("V0M");
  //  cout << lPercentiles << endl;
  }else{
  AliInfo("Didn't find MultSelection!"); 
  }
  
        
  if ( lPercentiles > 199 ){
    PostData(1,fOutputList );
    PostData(2, fSignalTree );
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2); 
    PostData(5, fOutputList3);
    PostData(6, fOutputList4);     
    return;  
  }

  fHistEventMult->Fill(5);

  //event must not be tagged as pileup
  Bool_t isPileUpSpd=kFALSE;
  isPileUpSpd=fAOD->IsPileupFromSPD();
  if(isPileUpSpd){ 
    PostData(1,fOutputList );
    PostData(2, fSignalTree );
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2); 
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);     
    return;
  }
  fHistEventMult->Fill(6);

  Bool_t isSelectedInt7        = kFALSE;
  Bool_t isSelectedAny         = kFALSE;
  Bool_t isSelected            = kFALSE;
  Bool_t isSelectedMB          = kFALSE;
 
  if(fCollidingSystem == "pp"){
    isSelectedInt7        = (mask & AliVEvent::kINT7);
    isSelectedAny         = (mask & AliVEvent::kAnyINT);
    isSelectedMB          = (mask & AliVEvent::kMB);

    
    if(isSelectedInt7 ) isSelected = kTRUE;
    if(fYear == 2010 && isSelectedMB) isSelected = kTRUE;

  }
  else if(fCollidingSystem == "pPb"){
    isSelectedInt7        = (mask & AliVEvent::kINT7);
    isSelectedAny         = (mask & AliVEvent::kAny);

    if(isSelectedInt7 )
      isSelected = kTRUE;
  }
  
  if(isSelectedInt7)
    fHistEventMult->Fill(7);
  else if(isSelectedAny)
    fHistEventMult->Fill(8);
  
  if (isSelected)
    fHistEventMult->Fill(9) ; 
      
  if(!isSelected){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3,fBkgTree);
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);
    // cout << "event does not fulfil centrality selection criteria " << endl;          
    return;
  }
  
  //  cout << "event has passed selection criteria.... first and second particles to be analyzed ...."<< endl;

  Bool_t Generated=kTRUE; //TRUE if generated particles are analyzed
  Int_t labelPrimOrSec=0; 
  Int_t label = 0; //track label
  Float_t TrackLength=-999;
  AliAODTrack* track=0x0;
  AliAODTrack* trackPtTMax=0x0;
  AliAODTrack* globaltrackPtTMax=0x0;
  AliAODTrack *globaltrack = 0x0;
  Bool_t isV0=kFALSE;

  if(fReadMCTruth){
    fMCEvent= MCEvent();
    if (fMCEvent){
      ProcessMCParticles(Generated, track, labelPrimOrSec, lPercentiles, isV0, 0, 0, fIshhCorr, 0,0);
    }
  }
  
  //Collect events according to their centrality and z coordinate of PV
  const Float_t bfield = (InputEvent())->GetMagneticField();
  int fieldsign;
  if (bfield >=0.) fieldsign = 1;
  else fieldsign = -1;
 
  int zBin=0;
  int centralityBin=0;
    
  double zStep=2*10/double(fzVertexBins), zStart=-10.;
    
  for (int i=0; i<fzVertexBins; i++) {
    if ((lBestPrimaryVtxPos[2] > zStart+i*zStep) && (lBestPrimaryVtxPos[2] < zStart+(i+1)*zStep)) {
      zBin=i;
      break;
    }
  } 

  if(lPercentiles < 5.0) centralityBin=19;  // changed <= with < to be consistent with histogram binning, except last bin 
  else if(lPercentiles < 10.0) centralityBin=18;
  else if(lPercentiles < 30.0) centralityBin=17;
  else if(lPercentiles < 50.) centralityBin=16;
  else if(lPercentiles <= 100.) centralityBin=15;

  if (((centralityBin+1) >fnMultBins) || ((zBin+1) > fzVertexBins)){ 
    //c cout<<" ##################  WARNING: I'm going to break bacause of dimensional issues ########################"<<endl;
  }

  fEventColl[zBin][centralityBin]->FifoClear();
  fEvt = fEventColl[zBin][centralityBin]->fEvt;


  //*****************create a global track to get the correct DCA*******************************
  for (Int_t igt = 0; igt < fTrackBufferSize; igt++) farrGT[igt] = -1;
 
  if (fFilterBitValue==128 || fFilterBitValue==512){
    // Read and store global tracks to retrieve PID information for TPC only tracks
    for (Int_t igt = 0; igt < iTracks; igt++) {
      globaltrack = (AliAODTrack*) fAOD->GetTrack(igt);

      if (!globaltrack) continue; 
      if (globaltrack->GetID()<0 ) continue;
      if (!globaltrack->IsOn(AliAODTrack::kTPCrefit)) continue; // there are such tracks with no TPC clusters

      // Check id is not too big for buffer
      if (globaltrack->GetID()>=fTrackBufferSize) {
	printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n",globaltrack->GetID(),fTrackBufferSize);
	fHistTrackBufferOverflow->Fill(1);
	continue;
      }

      //    if ( !(globaltrack->GetFilterMap()) ) { 
      //        cout<<" No filter map for this global track!!  "<<globaltrack->GetFilterMap()<<endl;
      //        continue;
      //    }

      if ( globaltrack->GetTPCNcls()<=0   ) { // such tracks do not have the TPC refit either, filter map is 2 --> ITS constrained
	//      cout<<" No TPC cl for this global track!!  "<<igt<<endl;
	//      if (!globaltrack->IsOn(AliAODTrack::kTPCrefit)) cout<<" ... and also no tpc refit  "<<globaltrack->GetFilterMap()<<endl;
	continue;
      }
      //cout<<" The array will contain "<<igt<<" , it contains "<<farrGT[globaltrack->GetID()]<<endl; 

      // Warn if we overwrite a track
      if (farrGT[globaltrack->GetID()]>=0) { // Two tracks same id  --> checked that it never happens
	//      cout<<" Array already filled "<<farrGT[globaltrack->GetID()]<<endl;
      } else { 
	farrGT[globaltrack->GetID()] = igt;           // solution adopted in the femto framework
	//    cout<<" Set array, now it contains  "<<farrGT[globaltrack->GetID()]<<endl;
      }
    }
  } //end if on filterbit value
  globaltrack = 0x0; 

  //***********************************************************************************
  //-----------------------------------LOOP OVER THE TRACKS
  //***********************************************************************************
  AliVTrack *vtrackg = 0x0;
  AliVTrack *vtrack = 0x0;
  AliVTrack *vtrackgPos = 0x0;
  AliVTrack *vtrackgNeg = 0x0;
  AliExternalTrackParam etp1;
  AliExternalTrackParam etpPos;
  AliExternalTrackParam etpNeg;

  Float_t nTPCCrossedRows=0.;
  float rationCrnFind=0;
  Float_t  FractionSharedTPCCluster=0;
  Float_t nsigmaTOFj=3;
  Float_t nsigmaTPCj=3;
  Int_t charge=0;
  Int_t NumberCharged=0;
  Int_t NumberFirstParticleAll=0;
  Int_t NumberFirstParticle=0;
  Int_t NumberFirstParticleAllPt=0;
  Int_t NumberFirstParticleAllPtMC=0;
  Int_t NumberFirstParticleMC=0;
  Int_t NumberFirstParticle_finale=0;
  Int_t NumberSecondParticleRecoTrue=0;
  Int_t NumberSecondParticle=0;
  Int_t NumberSecondParticleMC=0;
  Int_t NumberSecondParticleAll=0;
  Int_t NumberSecondParticleMCNoAssoc=0;
  Int_t NumberSecondParticleNoAssoc=0;
  Double_t Ptintermediate=0;
  Double_t selectedtrackID=0;
  Int_t pos0or1=0;
  Int_t neg0or1=0;
  Double_t dzglobal[2] = {-999.,-999.};
  Double_t dz[2] = {-999.,-999.};
  Double_t d[2] = {-999.,-999.};
  Double_t dzg[2]= {-999.,-999.}; Double_t covarg[3]={-999.,-999.,-999.};
  Double_t dzgPos[2]= {-999.,-999.}; Double_t covargPos[3]={-999.,-999.,-999.};
  Double_t dzgNeg[2]= {-999.,-999.}; Double_t covargNeg[3]={-999.,-999.,-999.};
  Double_t dzgPtTMax=-999;
  Float_t DCAxy=-999.;
  Float_t DCAz=-999.;
  Float_t labelPtTMax=0;
  Float_t ptTriggerMinimoDati=10000;
  Double_t ptTriggerMassimoDati=0;
  Double_t ptTriggerMassimoAll=0;
  Double_t ptTriggerMassimoDatiBis=0;
  Float_t etaTriggerMassimoDati=0;
  Float_t phiTriggerMassimoDati=0;
  Int_t   PdgCodeTrackPtMax=0;
  Int_t TriggerPdgCode=0;

  //begin loop for trigger particles   
  if (!(fReadMCTruth && !isEfficiency && !isHybridMCTruth)){
  for(Int_t i=0; i < iTracks; i++) {

    track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));        
    fHistTrack->Fill(1);
    if(!track) continue;

    //to know if trigger particle is primary, secondary,...                                   
    if(fReadMCTruth){
      if (fMCEvent){
	TClonesArray* AODMCTrackArrayTrigg =0x0;
	AODMCTrackArrayTrigg = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArrayTrigg == NULL){
	  return;
	  Printf("ERROR: stack not available");
	}

	AliAODMCParticle* particleTrig = static_cast<AliAODMCParticle*>(AODMCTrackArrayTrigg->At(TMath::Abs(track->GetLabel())));
	if (!particleTrig) continue;
	if(particleTrig->IsPhysicalPrimary())     labelPrimOrSec=1;
	else if(particleTrig->IsSecondaryFromWeakDecay())      labelPrimOrSec=2;
	else if(particleTrig->IsSecondaryFromMaterial())      labelPrimOrSec=3;
	else labelPrimOrSec=4;
	TriggerPdgCode = particleTrig->GetPdgCode();

      }
    }

    if(!track->TestFilterBit(fFilterBitValue)) continue; 
    fHistTrack->Fill(2);
    
    if(TMath::Abs(track->Eta())>fEtaTrigger)  continue;
    fHistTrack->Fill(3);

    NumberCharged++; //all charged particles
    
    if(track->Chi2perNDF()>4.)continue;
    fHistTrack->Fill(4);
    
    nTPCCrossedRows=track->GetTPCNCrossedRows();
    if(nTPCCrossedRows<80) continue;
    fHistTrack->Fill(5); 
        
    rationCrnFind=nTPCCrossedRows/track->GetTPCNclsF();
    if(rationCrnFind<0.8)  continue;
    fHistTrack->Fill(6); 

    Float_t lTrackLength = -1;
    Float_t lTrackLengthBis = -1;
    lTrackLength = GetLengthInActiveZone( track, /*1,*/ 2.0, 220.0, fAOD->GetMagneticField());
    lTrackLengthBis = GetLengthInActiveZone( track, /*1,*/ 0., 250.0, fAOD->GetMagneticField()); //test value

    Float_t     CrossedRowsOverLengthBis=    track->GetTPCClusterInfo(2,1)/lTrackLength;
    Float_t     CrossedRowsOverLength=(Float_t)nTPCCrossedRows/lTrackLength; //->both work in the same way

    fHistLengthvsCrossedRows ->Fill(  (Float_t)nTPCCrossedRows, lTrackLength );

    if (lTrackLength<90) continue;
    fHistTrack->Fill(7);  

    if (CrossedRowsOverLength< 0.8) continue;
    fHistTrack->Fill(8); 
    fHistLengthvsCrossedRowsAfterSel ->Fill(  (Float_t)nTPCCrossedRows, lTrackLength);

    //cluster selection (not applied)
    Int_t nClustersTPC = -1;
    //if(fCutRequireTPCStandAlone) {
    // nClustersTPC = track->GetTPCNclsIter1();
    //  }
    //  else {
    nClustersTPC = track->GetTPCncls(0);
    //}
    Int_t nClustersTPCShared = track->GetTPCnclsS();
    Float_t fracClustersTPCShared = -1.;
    if (nClustersTPC!=0) {
      fracClustersTPCShared = Float_t(nClustersTPCShared)/Float_t(nClustersTPC);
    }
    fHistFractionSharedTPCClusters->Fill(fracClustersTPCShared);
    //  if (fracClustersTPCShared > 0.2) continue;
    //    fHistTrack->Fill(9);


    //Golden cut (not applied)
    //  if ((track->GetChi2TPCConstrainedVsGlobal())>36) continue;
    //    fHistTrack->Fill(7);

    if((track->Charge())==0) continue;
    fHistTrack->Fill(9);
	  
    if (fFilterBitValue==128 || fFilterBitValue==512){ 
      //!!!!!!!!!information on DCA is correct only if taken in this way for tracks with FB128!!!!!!!!
      if (-track->GetID()-1 >= fTrackBufferSize) {
	printf ("Exceeding buffer size!!");
	continue;
      }

      vtrackg = fAOD->GetTrack(farrGT[-track->GetID()-1]);
      if (!vtrackg) {
	printf ("No global info! iTrack %d, ID %d\n",i,track->GetID());
	continue;
      }
      if (farrGT[-track->GetID()-1]>=iTracks || farrGT[-track->GetID()-1]<0) { /*cout<<"This index is out of range!!"<<farrGT[-track->GetID()-1]<<endl;*/ continue;}
      globaltrack = dynamic_cast<AliAODTrack*>(vtrackg); 
      if(!globaltrack) AliFatal("Not a standard AOD");

      dzg[0]= -999.; dzg[1] = -999.;
      etp1.CopyFromVTrack(vtrackg); 
      etp1.PropagateToDCA(vertexmain,(InputEvent())->GetMagneticField(),100.,dzg,covarg);    
      //      cout << "Propagate to DCA method : DCAxy " << dzg[0] << " DCAz " << dzg[1] << endl;

      /* does not seem to work properly 
      vtrackg->GetImpactParameters(&DCAxy, &DCAz);
      dz[0] = DCAxy;   
      dz[1] = DCAz; 
      cout << "GetImpactParamter method : DCAxy " << dz[0] << " DCAz " << dz[1] << endl;
      */
	
      /* The following method gives equal results to propagate to DCA     
      Float_t fPosition[2]={-999.};
      static_cast<AliAODTrack*>(fAOD->GetTrack(farrGT[-track->GetID()-1]))->GetPosition(fPosition);
      cout << "z position pos " << fPosition[2] << " z position primary vertex " << lBestPrimaryVtxPos[2] << endl;
      cout << "y position pos " << fPosition[1] << " y position primary vertex " << lBestPrimaryVtxPos[1] << endl;
      cout << "x position pos " << fPosition[0] << " x position primary vertex " << lBestPrimaryVtxPos[0] << endl;
      cout << fPosition[2]- lBestPrimaryVtxPos[2] << endl; //this value is = to dzg[1]
      cout << sqrt(pow (fPosition[0]- lBestPrimaryVtxPos[0],2) + pow(fPosition[1]- lBestPrimaryVtxPos[1],2)) << endl; //this value is = to dzg[0]
      */
    } //end if on filterbit value

    //    cout << "DCAxy trigger particle " << dzg[0] << " DCAz trigger particle " << dzg[1] << endl;	  

    if (fFilterBitValue!=128 && fFilterBitValue!=512){

      /* does not seem to work properly
      DCAxy=-999.;      DCAz =-999.;
      track->GetImpactParameters(&DCAxy, &DCAz);
      dz[0] = DCAxy;   
      dz[1] = DCAz; 
      cout << "GetImpactParamter method : DCAxy " << dz[0] << " DCAz " << dz[1] << endl;
      */

      //It seems that for FB 256 the two following methods (my method and propagate t DCA) give the same results (which seem reasonable: both DCAz and DCA yx are < 0.2 and have the correct shape), while the previous method (GetImpactParameter) give smaller results. 
      dzg[0]= -999.; dzg[1] = -999.;
      etp1.CopyFromVTrack(track); 
      etp1.PropagateToDCA(vertexmain,(InputEvent())->GetMagneticField(),100.,dzg,covarg);    
      //      cout << "Propagate to DCA method : DCAxy " << dzg[0] << " DCAz " << dzg[1] << endl;

      /*      Float_t fPosition[2]={-999.};
      track->GetPosition(fPosition);
      cout << "z position pos " << fPosition[2] << " z position primary vertex " << lBestPrimaryVtxPos[2] << endl;
      cout << "y position pos " << fPosition[1] << " y position primary vertex " << lBestPrimaryVtxPos[1] << endl;
      cout << "x position pos " << fPosition[0] << " x position primary vertex " << lBestPrimaryVtxPos[0] << endl;
      cout << "DCAz my method " << fPosition[2]- lBestPrimaryVtxPos[2] << endl; //this value is = to dzg[1]
      cout << "|DCAxy| my method " << sqrt(pow (fPosition[0]- lBestPrimaryVtxPos[0],2) + pow(fPosition[1]- lBestPrimaryVtxPos[1],2)) << endl; //this value is = to dzg[0]
      */
    }

    //    fHistDCAxym1->Fill(dz[0]);
    //    fHistDCAzm1->Fill(dz[1]);

    fHistDCAxym2->Fill(dzg[0]);
    fHistDCAzm2 ->Fill(dzg[1]);

    //it seems that Propagate to DCA give the best results for all FB:
    //    if (fFilterBitValue==128){
      dzglobal[0]=dzg[0];
      dzglobal[1]=dzg[1];
      /*
    }
    else{
      dzglobal[0]=dz[0];
      dzglobal[1]=dz[1];
    }
      */

    if(TMath::Abs(dzglobal[0])> (0.0105 + 0.0350/pow(track->Pt(),1.1))) continue;
    fHistTrack->Fill(10);
    if(TMath::Abs(dzglobal[1])> 2.) continue;
    fHistTrack->Fill(11);

    NumberFirstParticleAllPt++; 

    if(track->Pt()< ptTriggerMinimoDati) ptTriggerMinimoDati=track->Pt(); 
    if(track->Pt()> ptTriggerMassimoDati){
      ptTriggerMassimoDati =track->Pt(); 
      etaTriggerMassimoDati=track->Eta(); 
      phiTriggerMassimoDati=track->Phi(); 
      labelPtTMax=track->GetLabel();
      dzgPtTMax = dzglobal[1];
      trackPtTMax = static_cast<AliAODTrack*>(fAOD->GetTrack(i));        
      globaltrackPtTMax = dynamic_cast<AliAODTrack*>(vtrackg); 
      ptTriggerMassimoDatiBis=	trackPtTMax->Pt();
      PdgCodeTrackPtMax = labelPrimOrSec;
    }

    if (track->Pt()>Ptintermediate){
      Ptintermediate=track->Pt();
      selectedtrackID= track->GetID();
      NumberFirstParticle_finale= NumberFirstParticle;     
    }

    if(track->Pt()> fminPtj && track->Pt()<fmaxPtj){ //Here I select trigger particle with pT>fminPtj and save their characteristics
      NumberFirstParticle++;   
      if((!fReadMCTruth)|| (fReadMCTruth && isEfficiency) || (fReadMCTruth && isHybridMCTruth) ){
	//save first particle information (leading particle)
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fCharge       = track->Charge();
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fPt           = track->Pt();
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fEta          = track->Eta();
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fPhi          = track->Phi();
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fTheta        = track->Theta();
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fDCAz         = dzglobal[1];
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fDCAxy        = dzglobal[0];
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fMultiplicity = lPercentiles;
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fZvertex      = lBestPrimaryVtxPos[2];
	fEvt->fReconstructedFirst[NumberFirstParticle-1].isP           = labelPrimOrSec;
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fPDGcode      = TriggerPdgCode;
      }
      fHistPtvsMultBefAll->Fill(track->Pt(), lPercentiles);
    }

  }//end loop for trigger particles
}//done only if not MC truth

  TClonesArray* AODMCTrackArray =0x0;  
  Float_t ptTriggerMinimoMC=10000;
  Double_t ptTriggerMassimoMC=0;
  Double_t etaTriggerMassimoMC=0;
  Double_t phiTriggerMassimoMC=0;
  TriggerPdgCode=0;
  
  //begin loop for trigger particles (MC truth analysis)
  if(fReadMCTruth){
    if (fMCEvent){
      //      cout << " loop for trigger particles.. MC truth " << endl;
      AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      if (AODMCTrackArray == NULL){
	return;
	Printf("ERROR: stack not available");
      }
      for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
       
	AliAODMCParticle* trParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
       
	if (!trParticle) continue;       
	if((trParticle->Charge())==0)continue;
	if(TMath::Abs(trParticle->Eta())>fEtaTrigger)continue; //I need to select particles within this eta range!
	if (!(trParticle->IsPhysicalPrimary()))continue; 
	if (!(TMath::Abs(trParticle ->GetPdgCode()) == 211 || TMath::Abs(trParticle ->GetPdgCode()) == 321|| TMath::Abs(trParticle ->GetPdgCode()) ==2212 || TMath::Abs(trParticle ->GetPdgCode()) == 11|| TMath::Abs(trParticle ->GetPdgCode()) == 13)) continue;
	NumberFirstParticleAllPtMC++; 
	if(trParticle->Pt()<= fminPtj || trParticle->Pt()>=fmaxPtj)continue;
	if(trParticle->Pt()< ptTriggerMinimoMC) ptTriggerMinimoMC=trParticle->Pt();
	if(trParticle->Pt()> ptTriggerMassimoMC){
	  ptTriggerMassimoMC=trParticle->Pt();
	  etaTriggerMassimoMC=trParticle->Eta();
	  phiTriggerMassimoMC=trParticle->Phi();
	  TriggerPdgCode =trParticle ->GetPdgCode();
	}
	NumberFirstParticleMC++;
	if (isEfficiency || isHybridMCTruth) continue;	  
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fCharge       = trParticle->Charge();
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fPt           = trParticle->Pt();
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fEta          = trParticle->Eta();
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fPhi          = trParticle->Phi();
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fTheta        = trParticle->Theta();
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fDCAz         = 0;
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fDCAxy        = 0;
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fMultiplicity = lPercentiles;
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fZvertex      = lBestPrimaryVtxPos[2];
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].isP           = 1;
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fPDGcode     = trParticle->GetPdgCode();

      }
      //      cout << "I found " << 	NumberFirstParticleMC << " trigger particles (pT > 3 GeV/c, not only the highest pt one)" << endl;
      //      cout << " pt of trigger particle " << 	  ptTriggerMassimoMC << endl;
    }
  }  //end loop for trigger particles (MC truth analysis)
  //  cout << " end of loop for trigger particles (MCtruth)" << endl;

  if ((!fReadMCTruth || (fReadMCTruth &&isEfficiency) || (fReadMCTruth && isHybridMCTruth))&&  ptTriggerMassimoDati!=0 )  fHistPtMaxvsMultBefAll->Fill(ptTriggerMassimoDati, lPercentiles);
  if (fReadMCTruth && !isEfficiency && !isHybridMCTruth &&  ptTriggerMassimoMC!=0)  fHistPtMaxvsMultBefAll->Fill(ptTriggerMassimoMC, lPercentiles);  

  if(fReadMCTruth){ //to determine "event loss"
    if(ptTriggerMassimoDati!=0)  fHistPtMaxvsMultBefAllReco->Fill(ptTriggerMassimoDati, lPercentiles);  
    if(ptTriggerMassimoMC!=0)    fHistPtMaxvsMultBefAllGen ->Fill(ptTriggerMassimoMC, lPercentiles);  
  }

  fHistPtTMinBefAll->Fill(ptTriggerMinimoDati);
  fHistPtTMinBefAllMC->Fill(ptTriggerMinimoMC);

  fHistPtTMaxBefAllCfrDataMC->Fill(ptTriggerMassimoDati,ptTriggerMassimoMC);
  fHistTrigger->Fill(NumberFirstParticle);
  fHistTriggerMCTruth->Fill(NumberFirstParticleMC);
  fHistMultvsTriggerBefAll->Fill(NumberFirstParticle, lPercentiles);
  fHistMultvsTriggerMCTruthBefAll->Fill(NumberFirstParticleMC, lPercentiles);

  fHistTrack->AddBinContent(12, NumberFirstParticle);
  fHistTrack->AddBinContent(13, NumberFirstParticleMC);
  if (NumberFirstParticle>0) fHistEventMult->Fill(10);   
  if (NumberFirstParticleMC>0) fHistEventMult->Fill(11);   
  if(NumberFirstParticle>1)fHistEventMult->Fill(12); 
  if(NumberFirstParticleMC>1)fHistEventMult->Fill(13); 
  if(NumberFirstParticleMC<NumberFirstParticle) fHistEventMult->Fill(14); 
  if(NumberFirstParticleMC==0 && NumberFirstParticle==1) fHistEventMult->Fill(15); 
  if(NumberFirstParticleMC!=0 && NumberFirstParticle!=0) fHistEventMult->Fill(16); 

  
  if((fReadMCTruth && isEfficiency) || (!fReadMCTruth) || (fReadMCTruth && isHybridMCTruth)) {NumberFirstParticleAll=NumberFirstParticle; ptTriggerMassimoAll=ptTriggerMassimoDati;}
  else if (fReadMCTruth && !isEfficiency && !isHybridMCTruth) {NumberFirstParticleAll=NumberFirstParticleMC; ptTriggerMassimoAll=ptTriggerMassimoMC; NumberFirstParticleAllPt = NumberFirstParticleAllPtMC;}
  if(NumberFirstParticleAll==0) fHistNumberChargedNoTrigger->Fill(lPercentiles, NumberCharged,0 );
  if(NumberFirstParticleAll!=0) fHistNumberChargedTrigger->Fill(lPercentiles, NumberCharged,ptTriggerMassimoAll );
  fHistNumberChargedAllEvents->Fill(lPercentiles, NumberCharged,ptTriggerMassimoAll );

  if(NumberFirstParticleAllPt==0){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3,fBkgTree);
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);
    // cout  << "event does not have Trigger particles " << endl;          
    return;
  }
 
  //************Filling selected histograms for trigger particle efficiency calculation  ***************
  //only the highest pT trigger particle in each event is used, since this is the definition of trigger particle used in the angular correlation
  //! The selected histograms are filled also by trigger particles with pT< fminPtj, a cut on pT is therefore required in post processing phase
  Generated = kFALSE;     
  isV0=kFALSE;
  Int_t VPdgTrig[50]={0};
  Int_t VParticleTrigLabel[50]={0};
  if(fReadMCTruth && isEfficiency){
    if(fMCEvent){
      //      ProcessMCParticles(Generated, trackPtTMax, labelPrimOrSec, lPercentiles, isV0, dzgPtTMax,0, fIshhCorr);
      ProcessMCParticles(Generated, trackPtTMax, labelPrimOrSec, lPercentiles, isV0, dzgPtTMax,0, fIshhCorr, VPdgTrig, VParticleTrigLabel);
    }
  }
  //***********************************************************************************************************

  if(NumberFirstParticleAll==0){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3,fBkgTree);
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);
    //    cout  << "event does not have Trigger particles with pT> 3 GeV/c " << endl;          
    return;
  }
  fHist_multiplicity_EvwTrigger->Fill(lPercentiles);
  fHistTrack->AddBinContent(14, NumberFirstParticle);
  fHistTrack->AddBinContent(15, NumberFirstParticleMC);
  if (fReadMCTruth && !isEfficiency)  fHistTriggerCompositionMCTruth->Fill(TriggerPdgCode,1,ptTriggerMassimoMC);

  Int_t V0PDGCode=0;    
  Int_t hAssocPDGCode=0;
  Int_t labelPrimOrSecV0=0;
  //*****************************************************************************************************
  //***********************this part is executed only if hh correlation analysis is chosen***************
  if (fIshhCorr){
    //    cout << "\n\n**********\n**********\n************\nhh correlation analysis chosen " << endl;
    //MC generated charged particles associated
    isV0=kTRUE;
    Generated=kTRUE;

    if(fReadMCTruth){
      if (fMCEvent){
	ProcessMCParticles(Generated, trackPtTMax, labelPrimOrSec, lPercentiles, isV0, 0, ptTriggerMassimoDati, fIshhCorr, VPdgTrig, VParticleTrigLabel);
      }
    }

    //LOOP for charged particles as associated particles
    bool skipAssoc=kFALSE; 
    if (!(fReadMCTruth && !isEfficiency)){
    for(Int_t i=0; i < iTracks; i++) {
      track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));        
      fHistTrackAssoc->Fill(1);
      if(!track) continue;

      if(!track->TestFilterBit(fFilterBitValue)) continue;
      fHistTrackAssoc->Fill(2);
    
      if(TMath::Abs(track->Eta())>fEtahAssoc)  continue;
      fHistTrackAssoc->Fill(3);
    
      if(track->Chi2perNDF()>4.)continue;
      fHistTrackAssoc->Fill(4);
    
      nTPCCrossedRows=track->GetTPCNCrossedRows();
      if(nTPCCrossedRows<80) continue;
      fHistTrackAssoc->Fill(5);
        
      rationCrnFind=nTPCCrossedRows/track->GetTPCNclsF();
      if(rationCrnFind<0.8)  continue;
      fHistTrackAssoc->Fill(6);
   
      Float_t lTrackLengthAssoc = GetLengthInActiveZone( track, /*1,*/ 2.0, 220.0, fAOD->GetMagneticField());
      if (lTrackLengthAssoc<90) continue;
      fHistTrackAssoc->Fill(7);  
      
      Float_t      CrossedRowsOverLengthAssoc=(Float_t)nTPCCrossedRows/lTrackLengthAssoc;
      if (CrossedRowsOverLengthAssoc< 0.8) continue;
      fHistTrackAssoc->Fill(8); 

      if((track->Charge())==0) continue;
      fHistTrackAssoc->Fill(9);

      if (fFilterBitValue==128 || fFilterBitValue==512){
	//	cout << "track ID (FB 128 selected" << 	track->GetID()<< endl;;
	
	if (-track->GetID()-1 >= fTrackBufferSize) {
	  printf ("Exceeding buffer size!!");
	  continue;
	}
	vtrackg = fAOD->GetTrack(farrGT[-track->GetID()-1]);
	if (!vtrackg) {
	  printf ("No global info! iTrack %d, ID %d\n",i,track->GetID());
	  continue;
	}
	if (farrGT[-track->GetID()-1]>=iTracks || farrGT[-track->GetID()-1]<0) { /*cout<<"This index is out of range!!"<<farrGT[-track->GetID()-1]<<endl;*/  continue;}
	globaltrack = dynamic_cast<AliAODTrack*>(vtrackg); 
	if(!globaltrack) AliFatal("Not a standard AOD");

	dzg[0]= -999.; dzg[1] = -999.;
	etp1.CopyFromVTrack(vtrackg); 
	etp1.PropagateToDCA(vertexmain,(InputEvent())->GetMagneticField(),100.,dzg,covarg);    

      } //end if on filterbit value 

      if (fFilterBitValue!=128 && fFilterBitValue!=512){
	dzg[0]= -999.; dzg[1] = -999.;
	etp1.CopyFromVTrack(vtrackg); 
	etp1.PropagateToDCA(vertexmain,(InputEvent())->GetMagneticField(),100.,dzg,covarg);    

	/* does not seem to work properly 
	DCAxy=-999.;      DCAz =-999.;
	track->GetImpactParameters(&DCAxy, &DCAz);
	dz[0] = DCAxy;   
	dz[1] = DCAz; 
	*/
      }

	dzglobal[0]=dzg[0];
	dzglobal[1]=dzg[1];

      if(TMath::Abs(dzglobal[0])> (0.0105 + 0.0350/pow(track->Pt(),1.1))) continue;
      fHistTrackAssoc->Fill(10);
      if(TMath::Abs(dzglobal[1])> 2.) continue;
      fHistTrackAssoc->Fill(11);

      if(!(track->Pt()> fminPthAssoc && track->Pt()<fmaxPthAssoc)) continue;
      skipAssoc=kFALSE;
      fHistTrackAssoc->Fill(12);

      if (track->Pt()>=ptTriggerMassimoDati) skipAssoc=kTRUE;
      if (skipAssoc) continue;
      NumberSecondParticle++;

      //I fill histograms with associated charged particles which have been selected****************
      labelPrimOrSecV0=0;


      //MC analysis
      AliAODMCParticle *VParticlehAssoc[50]={0};
      Int_t VPdghAssoc[50]={0};
      Int_t VParticlehAssocLabel[50]={0};
      Bool_t IsCommonParton =0;
      Int_t TrigLabel=0;
      Int_t hAssocLabel=0;
      Int_t PdgCodeCommon=-10;

      if(fReadMCTruth){
	if (fMCEvent){
	  AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	  if (AODMCTrackArray == NULL){
	    return;
	    Printf("ERROR: stack not available");
	  }
	
	  AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(track->GetLabel())));

	  //Common Parton -> ****************************
	  VParticlehAssoc[0]=particle; 
	  VPdghAssoc[0]=       VParticlehAssoc[0]->GetPdgCode();
	  VParticlehAssocLabel[0]=       VParticlehAssoc[0]->GetLabel();

	  //	  cout << "assoc: " ;
	  for (Int_t i=0; i<50; i++){
	    VParticlehAssocLabel[i+1]=VParticlehAssoc[i]->GetMother();
	    VParticlehAssoc[i+1] = static_cast<AliAODMCParticle*>(AODMCTrackArray-> At(TMath::Abs(VParticlehAssocLabel[i+1])));
	    VPdghAssoc[i+1] = VParticlehAssoc[i+1]->GetPdgCode();
	    //        cout << VPdghAssoc[i] << " (" << VParticlehAssocLabel[i] << ") " << "<-" ;                               
	    if ((VParticlehAssocLabel[i] ==1 || VParticlehAssocLabel[i] ==-1) && (VPdghAssoc[i]==2212) && (VParticlehAssocLabel[i] ==  VParticlehAssocLabel[i+1] )) {
	      //cout << endl;                                                                                               
	      break;
	    }
	  }

	  for (Int_t i=1; i<50; i++){ //I start from one since last element cannot be a parton but is a hadron              
	    if (IsCommonParton==1) break;
	    for (Int_t j=1; j<50; j++){
	      if ((VParticlehAssocLabel[i] == VParticleTrigLabel[j] ) &&  VParticleTrigLabel[j]!=0 && ( TMath::Abs(VPdghAssoc[i]) <=8 ||  TMath::Abs(VPdghAssoc[i]) ==21)) { //both Xi and Trigger particle have a common ancestor which has to be a quark or a gluon-> therefore te cascade comes form the jet defined by the trigger particle                                         
		IsCommonParton =1;
		TrigLabel=j;
		hAssocLabel=i;
		PdgCodeCommon =VPdghAssoc[hAssocLabel];
		break;
	      }
	    }
	  }
	  if (IsCommonParton)      {
	    fHistCommonParton->Fill(PdgCodeCommon, TrigLabel, hAssocLabel);
	    fHistIsCommonParton->Fill(1);
	  }
	  else {
	    fHistIsCommonParton->Fill(0);
	  }

	  //Common Parton <- ****************************

	  hAssocPDGCode= particle->GetPdgCode();
	  if(! (particle->IsPhysicalPrimary())){
	    fHistAssocPtRecovsPtGenNotPrim->Fill(particle->Pt(), track->Pt());
	  }
	  if(particle->IsPhysicalPrimary()){
	    fHistAssocComposition->Fill(particle->GetPdgCode(),1, ptTriggerMassimoDati);
	    fHistAssocPtRecovsPtGen->Fill(particle->Pt(), track->Pt());
	    if (particle->GetPdgCode()==211 || particle->GetPdgCode()==-211)      fHistAssocPtRecovsPtGenPion->Fill(particle->Pt(), track->Pt());
	    if (particle->GetPdgCode()==2212 || particle->GetPdgCode()==-2212)      fHistAssocPtRecovsPtGenProton->Fill(particle->Pt(), track->Pt());
	    if (particle->GetPdgCode()==321 || particle->GetPdgCode()==-321)      fHistAssocPtRecovsPtGenKaon->Fill(particle->Pt(), track->Pt());
	    /* 3D histos

	       fHistResolutionV0Pt->Fill(track->Pt()- particle->Pt(), lPercentiles, ptTriggerMassimoDati);
	       fHistResolutionV0Phi->Fill(track->Phi()- particle->Phi(), lPercentiles, ptTriggerMassimoDati);
	       fHistResolutionV0PhivsPt->Fill(track->Phi()- particle->Phi(), lPercentiles, track->Pt());
	       fHistResolutionV0PtvsPt->Fill(track->Pt()- particle->Pt(), lPercentiles, track->Pt());
	       fHistResolutionV0Eta->Fill(track->Eta()- particle->Eta(), lPercentiles, ptTriggerMassimoDati);
	    */
	    // 2D histos
	    fHistResolutionV0Pt->Fill(track->Pt()- particle->Pt(), ptTriggerMassimoDati);
	    fHistResolutionV0Phi->Fill(track->Phi()- particle->Phi(), ptTriggerMassimoDati);
	    fHistResolutionV0PhivsPt->Fill(track->Phi()- particle->Phi(), track->Pt());
	    fHistResolutionV0PtvsPt->Fill(track->Pt()- particle->Pt(), track->Pt());
	    fHistResolutionV0Eta->Fill(track->Eta()- particle->Eta(), ptTriggerMassimoDati);
	    //
	    if(  (TMath::Abs(dzglobal[1]) < 1.)) {
	      fHistSelectedV0PtTMaxPhi[0]->Fill(ptTriggerMassimoDati, track->Phi(), lPercentiles);
	      fHistSelectedV0PtTMaxEta[0]->Fill(ptTriggerMassimoDati, track->Eta(), lPercentiles);
	      fHistSelectedV0PtPtTMax[0]->Fill(track->Pt(),ptTriggerMassimoDati , lPercentiles);
	      fHistSelectedV0PtEta[0]->Fill(track->Pt(), track->Eta(), lPercentiles);
	      fHistSelectedGenV0PtPtTMax[0]->Fill(particle->Pt(),ptTriggerMassimoDati , lPercentiles);
	    } 

	    labelPrimOrSecV0=1;
	    NumberSecondParticleRecoTrue++;
	  }
	  else if(particle->IsSecondaryFromWeakDecay())      labelPrimOrSecV0=2;
	  else if(particle->IsSecondaryFromMaterial())      labelPrimOrSecV0=3;
	  else labelPrimOrSecV0=4;
	    
	  if (!(particle->IsPhysicalPrimary())) 	    fHistAssocComposition->Fill(particle->GetPdgCode(),0.1, ptTriggerMassimoDati);

	  for (Int_t m =0; m<5;m++){
	    if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	      for(Int_t p=1; p<=4; p++){
		if (labelPrimOrSecV0==p) {
		  if(  (TMath::Abs(dzglobal[1]) < 1.)) {
		    fHistPrimaryV0[m][0]->Fill(p, track->Pt(), ptTriggerMassimoDati);   
		  }
		}
	      }
	    }
	  }
	  for(Int_t p=1; p<=4; p++){
	    if (labelPrimOrSecV0==p){
	      if(  (TMath::Abs(dzglobal[1]) < 1.)) {
		fHistPrimaryV0[5][0]->Fill(p, track->Pt(), ptTriggerMassimoDati);   
	      }
	    }
	  }
	}
      }
      //save second particle information
      if((!fReadMCTruth)|| (fReadMCTruth && isEfficiency)){
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sCharge       = track->Charge();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sPt           = track->Pt();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sEta          = track->Eta();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sPhi          = track->Phi();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sTheta        = track->Theta();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sDCAz         = dzglobal[1];
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sDCAxy        = dzglobal[0];
	fEvt->fReconstructedSecond[NumberSecondParticle-1].isP           = labelPrimOrSecV0;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sPDGcode      = hAssocPDGCode;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sAssocOrNot     = skipAssoc;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sIsCommonParton      = IsCommonParton;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sPdgCommonParton     = PdgCodeCommon;

	fHistPthAssoc->Fill(track->Pt());
      }

    }  //end loop for charged particles as associated particles
    } //done only if it is not MC truth

    //begin MC loop for charged particles as associated particles
    if(fReadMCTruth){
      if (fMCEvent){
	AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArray == NULL){
	  return;
	  Printf("ERROR: stack not available");
	}
	for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
       
	  AliAODMCParticle* trParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
       
	  if (!trParticle) continue;     
	  if((trParticle->Charge())==0)continue;
	  if(TMath::Abs(trParticle->Eta())>fEtahAssoc)continue; //I need to select particles within this eta range!
	  if (!(trParticle->IsPhysicalPrimary()))continue; 
	  if(trParticle->Pt()<= fminPthAssoc || trParticle->Pt()>=fmaxPthAssoc)continue;
	  bool	skipAssoc_MC=kFALSE;
	  if ((trParticle->Pt())>=ptTriggerMassimoMC) {
	    skipAssoc_MC=kTRUE;
	  }
	  if (skipAssoc_MC)      continue;
	  //qui potri introdurre istogramma charged particles selezionate in funzione di pT trigger max, pt h associata e centralità
	  //	  fHistPthAssoc->Fill(trParticle->Pt());
	  NumberSecondParticleMC++;
	  if(isEfficiency) continue;

	  if (fReadMCTruth && !isEfficiency)  fHistAssocCompositionMCTruth->Fill(trParticle->GetPdgCode(),1, ptTriggerMassimoMC);
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sCharge       = trParticle->Charge();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPt           = trParticle->Pt();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sEta          = trParticle->Eta();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPhi          = trParticle->Phi();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sTheta        = trParticle->Theta();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDCAz         = dzglobal[1];
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDCAxy        = dzglobal[0];
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].isP           = 1;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPDGcode        = trParticle->GetPdgCode();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sAssocOrNot     = skipAssoc_MC;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sIsCommonParton        = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPdgCommonParton       = 0;

	}
      }
    } //end MC truth loop for charged particles as associated
  }


  //*****************************end hh section*************************************************

  if (!fIshhCorr){
    //    cout << "\n\n**********\n**********\n************\nhV0 correlation analysis chosen " << endl;
    //    cout << "hV0 correlation analysis chosen" << endl;

    Float_t ycut[2]={0.5, 0.5};
    Float_t PosEtaCut[2]={0.8, 0.8};
    Float_t NegEtaCut[2]={0.8, 0.8};
    Float_t nsigmaTPC=3;
    Float_t nsigmaTOF=3;
    Float_t kMaxTPCSigmaPion=3;
    Double_t rapidityV0=0;
    Double_t EV0[2]={0};
    bool goodPiPlus=kFALSE;
    bool goodPiMinus=kFALSE;
    bool goodPiPlusTPC=kFALSE;
    bool goodPiMinusTPC=kFALSE;
    bool goodPiPlusTOF=kFALSE;
    bool goodPiMinusTOF=kFALSE;
    Float_t kMaxTOFSigmaPion=3;
    Float_t kTOFLow[2]={1000, 1000}; 
    Float_t kMaxDCA[2]={0.6, 1000}; //was 0.5 for K0s
    Float_t kMinCosAngle[2]={0.97, 0.997}; //was 0.995 for K0s
    Float_t kMaxDCADaughters[2]={1, 1};
    Float_t kMinDCAPrimary[2]={0.05, 0.06}; //was 0.06 for K0s
    Float_t kMinDL[2]={0.5, 0.5};
    Float_t kMinRadius[2]={0.5, 0.5};
    Float_t kctauval[2]={40, 30}; //was 20 for K0s
    Float_t kctau[2]={0, 0};
    Float_t Mass[2]={0.497611, 1.115683};
    Int_t ParticleType=0;

    Int_t labelPos=0;
    Int_t labelNeg=0;
    AliAODMCParticle* particlePos;
    AliAODMCParticle* particleNeg;
    Int_t PdgPos=0;
    Int_t PdgNeg=0;
    Int_t labelMotherPos=0;
    Int_t labelMotherNeg=0;
    AliAODMCParticle* MotherPos;
    AliAODMCParticle* MotherNeg;
    Int_t PdgMotherPos=0;
    Int_t PdgMotherNeg=0;
    TLorentzVector vPos;
    TLorentzVector vNeg;
    TLorentzVector vPhoton;
    Double_t MassPhoton;
    Double_t MassPhoton2;
    Double_t MassElectron= 0.0005109989461;
    Double_t massLambda = 1.115683;
    Int_t isaK0s=0;
    if(fV0=="K0s") ParticleType =0;
    if(fV0=="Lambda") ParticleType=1;

    //MC generated V0
    isV0=kTRUE;
    Generated=kTRUE;

    if(fReadMCTruth){
      if (fMCEvent){
	ProcessMCParticles(Generated, track, labelPrimOrSec, lPercentiles, isV0, 0, ptTriggerMassimoAll, fIshhCorr, VPdgTrig, VParticleTrigLabel);
      }
    }
  
    // cout << "\n \n  here I start the loop on v0s " << endl;
    if (!(fReadMCTruth && !isEfficiency)){
    for(Int_t i(0); i < V0Tracks; i++) {       
      isaK0s=0; //it will be put to 1 for true K0s in MC
      rapidityV0=0;
      EV0[2]={0};
      kctau[2]={0};
      goodPiPlus=kFALSE;
      goodPiMinus=kFALSE;

      AliAODv0* v0 = fAOD->GetV0(i);
      if(!v0) continue;       
      fHistEventV0->Fill(1);

      if(v0->GetOnFlyStatus()) continue;
      fHistEventV0->Fill(2);

      AliAODTrack* tempTrack=(AliAODTrack*)v0->GetDaughter(0);
      if(tempTrack->Charge()>0) {pos0or1=0; neg0or1=1;}
      else {pos0or1=1; neg0or1=0;}
      AliAODTrack* prongTrackPos=(AliAODTrack*)v0->GetDaughter(pos0or1);
      AliAODTrack* prongTrackNeg=(AliAODTrack*)v0->GetDaughter(neg0or1);

      if (!prongTrackPos || !prongTrackNeg) {
	Printf("ERROR: Could not retreive one of the daughter track");
	continue;
      }

      // Filter like-sign V0 (next: add counter and distribution)
      if ( prongTrackPos->GetSign() == prongTrackNeg->GetSign()) {
	continue;
      }

      labelPos=prongTrackPos->GetLabel();
      labelNeg=prongTrackNeg->GetLabel();
      
      AliAODMCParticle *VParticleV0[50]={0};
      Int_t VPdgV0[50]={0};
      Int_t VParticleV0Label[50]={0};
      Bool_t IsCommonParton =0;
      Int_t TrigLabel=0;
      Int_t V0Label=0;
      Int_t PdgCodeCommon=-10;

      if (fReadMCTruth){
	if (fMCEvent){
	  AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	  if (AODMCTrackArray == NULL){
	    return;
	    Printf("ERROR: stack not available");
	  }
	  particlePos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelPos)));
	  particleNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelNeg)));
	  if(labelPos>=0)	PdgPos = particlePos->GetPdgCode();
	  if(labelNeg>=0)	PdgNeg = particleNeg->GetPdgCode();
 
	  labelMotherPos=particlePos->GetMother();
	  labelMotherNeg=particleNeg->GetMother();
	  MotherPos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelMotherPos));
	  MotherNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelMotherNeg));
	
	  if (labelMotherPos>=0) PdgMotherPos = MotherPos->GetPdgCode();
	  if (labelMotherNeg>=0)	PdgMotherNeg = MotherNeg->GetPdgCode();
	

	  //Common Parton -> ****************************
	  VParticleV0[0]=MotherPos; 
	  VPdgV0[0]=       VParticleV0[0]->GetPdgCode();
	  VParticleV0Label[0]=       VParticleV0[0]->GetLabel();

	  //	  cout << "assoc: " ;
	  for (Int_t i=0; i<50; i++){
	    VParticleV0Label[i+1]=VParticleV0[i]->GetMother();
	    VParticleV0[i+1] = static_cast<AliAODMCParticle*>(AODMCTrackArray-> At(TMath::Abs(VParticleV0Label[i+1])));
	    VPdgV0[i+1] = VParticleV0[i+1]->GetPdgCode();
	    //        cout << VPdgV0[i] << " (" << VParticleV0Label[i] << ") " << "<-" ;                               
	    if ((VParticleV0Label[i] ==1 || VParticleV0Label[i] ==-1) && (VPdgV0[i]==2212) && (VParticleV0Label[i] ==  VParticleV0Label[i+1] )) {
	      //cout << endl;                                                                                               
	      break;
	    }
	  }

	  for (Int_t i=1; i<50; i++){ //I start from one since last element cannot be a parton but is a hadron              
	    if (IsCommonParton==1) break;
	    for (Int_t j=1; j<50; j++){
	      if ((VParticleV0Label[i] == VParticleTrigLabel[j] ) &&  VParticleTrigLabel[j]!=0 && ( TMath::Abs(VPdgV0[i]) <=8 ||  TMath::Abs(VPdgV0[i]) ==21)) { //both Xi and Trigger particle have a common ancestor which has to be a quark or a gluon-> therefore te cascade comes form the jet defined by the trigger particle                                         
		IsCommonParton =1;
		TrigLabel=j;
		V0Label=i;
		PdgCodeCommon =VPdgV0[V0Label];
		break;
	      }
	    }
	  }

	  if (IsCommonParton)      {
	    fHistCommonParton->Fill(PdgCodeCommon, TrigLabel, V0Label);
	    fHistIsCommonParton->Fill(1);
	  }
	  else {
	    fHistIsCommonParton->Fill(0);
	  }

	  //Common Parton <- ****************************

	}
      }

      if(ParticleType==0) EV0[ParticleType]= TMath::Sqrt(pow(v0->MassK0Short(),2)+ pow( v0->Px(),2)+ pow( v0->Py(),2)+ pow( v0->Pz(),2));
      if(ParticleType==1) EV0[ParticleType]= TMath::Sqrt(pow(v0->MassLambda(),2)+ pow( v0->Px(),2)+ pow( v0->Py(),2)+ pow( v0->Pz(),2));

      rapidityV0  = 0.5*TMath::Log(( EV0[ParticleType] + v0->Pz()) / (EV0[ParticleType]   - v0->Pz()) );
    
      fTreeVariablePtV0=TMath::Sqrt( pow( v0->Px(),2)+ pow( v0->Py(),2));
      fTreeVariableInvMassK0s=v0->MassK0Short();    
      fTreeVariableInvMassLambda=v0->MassLambda();    
      fTreeVariableInvMassAntiLambda=v0->MassAntiLambda();    
      kctau[ParticleType]=Mass[ParticleType]*v0->DecayLengthV0(lBestPrimaryVtxPos)/TMath::Sqrt( pow( v0->Px(),2)+ pow( v0->Py(),2) + pow(v0->Pz(),2));
      vPos.SetPxPyPzE(prongTrackPos->Px(),prongTrackPos->Py(),prongTrackPos->Pz(),sqrt(pow(prongTrackPos->P(),2) + pow(MassElectron,2)));
      vNeg.SetPxPyPzE(prongTrackNeg->Px(),prongTrackNeg->Py(),prongTrackNeg->Pz(),sqrt(pow(prongTrackNeg->P(),2) + pow(MassElectron,2)));
      vPhoton = vPos+vNeg;
      MassPhoton2=vPhoton.M2();
      MassPhoton=vPhoton.M();
      fHistMass2Photon->Fill(MassPhoton2);
   

      //-------------------------daughter cuts---------------

      Float_t nTPCCrossedRowspos=prongTrackPos ->GetTPCNCrossedRows();
      Float_t nTPCCrossedRowsneg=prongTrackNeg ->GetTPCNCrossedRows();
      Float_t rationCrnFindpos=nTPCCrossedRowspos/prongTrackPos->GetTPCNclsF();
      Float_t rationCrnFindneg=nTPCCrossedRowsneg/prongTrackNeg->GetTPCNclsF();

//      if(!prongTrackPos->TestFilterBit(1)) continue; //I do not use the Filterbit for V0 daughter tracks because this suppresses the efficiency
//      if(!prongTrackNeg->TestFilterBit(1)) continue;
      fHistEventV0->Fill(3);
      if(prongTrackPos->Chi2perNDF()>4.)continue;
      if(prongTrackNeg->Chi2perNDF()>4.)continue;
      fHistEventV0->Fill(4);

      Float_t       lV0Radius			= v0->RadiusV0();
      //      Float_t       lV0RadiusBis			= TMath::Sqrt(pow(v0->DecayVertexV0X(),2) + pow(v0->DecayVertexV0Y(),2)); //this works exactly as lV0Radius


      //TPC Refit for daugther tracks
      ULong_t pStatus    = prongTrackPos->GetStatus();
      ULong_t nStatus    = prongTrackNeg->GetStatus();

      if ((pStatus&AliAODTrack::kTPCrefit)    == 0) {
	AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!");
	continue;
      }
      if ((nStatus&AliAODTrack::kTPCrefit)    == 0) {
	AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!");
	continue;
      }
      fHistEventV0->Fill(5);

      //TPCClusters
      Int_t      lPosTPCClusters   = prongTrackPos->GetTPCNcls();
      Int_t      lNegTPCClusters   = prongTrackNeg->GetTPCNcls();
      if(lPosTPCClusters  < 50) {
	AliWarning("Pb / V0 Pos. track has less than 50 TPC clusters ... continue!");
	continue;
      }
      
      if(lNegTPCClusters  < 50) {
	AliWarning("Pb / V0 Neg. track has less than 50 TPC clusters ... continue!");
	continue;
      }
      fHistEventV0->Fill(6);


      //GetKinkIndex condition                                                                                                        
      Bool_t V0VarPosIsKink=kFALSE;
      Bool_t V0VarNegIsKink=kFALSE;
      if( prongTrackPos->GetKinkIndex(0)>0 ) V0VarPosIsKink = kTRUE;
      if( prongTrackNeg->GetKinkIndex(0)>0 ) V0VarNegIsKink = kTRUE;

      if (V0VarPosIsKink || V0VarNegIsKink) continue;
      fHistEventV0->Fill(7);

      //Tracklength selection
      Float_t lTrackLengthpos = -1;
      Float_t lTrackLengthneg = -1;
      lTrackLengthpos = GetLengthInActiveZone( prongTrackPos, 2.0, 220.0, fAOD->GetMagneticField());
      lTrackLengthneg = GetLengthInActiveZone( prongTrackNeg, 2.0, 220.0, fAOD->GetMagneticField());
      Float_t lPosTrackNcrOverLength = (Float_t)prongTrackPos->GetTPCClusterInfo(2,1)/(lTrackLengthpos-TMath::Max(lV0Radius-85.,0.));
      Float_t lNegTrackNcrOverLength = (Float_t)prongTrackNeg->GetTPCClusterInfo(2,1)/(lTrackLengthneg-TMath::Max(lV0Radius-85.,0.));
      
      if (lTrackLengthpos<90 || lTrackLengthneg<90) continue;
      fHistEventV0->Fill(8);

      if (lPosTrackNcrOverLength< 0.8 || lNegTrackNcrOverLength< 0.8) continue;
      fHistEventV0->Fill(9); 

     //TPC PID
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( prongTrackPos, (AliPID::EParticleType)2))< 3.) goodPiPlusTPC=kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( prongTrackNeg, (AliPID::EParticleType)2))< 3.) goodPiMinusTPC=kTRUE;

      //A. Balbino way of doing oob pileup:
      Float_t      fV0_NegTOFBunchCrossing = prongTrackNeg->GetTOFBunchCrossing(fAOD->GetMagneticField());
      Float_t      fV0_PosTOFBunchCrossing = prongTrackPos->GetTOFBunchCrossing(fAOD->GetMagneticField());
      Bool_t IsOOBPileUpBalbino=kTRUE;
      Float_t cutval_V0=-95;
      //      if (!(fV0_NegTrackStatus & AliESDtrack::kITSrefit) && !(fV0_PosTrackStatus & AliESDtrack::kITSrefit) && (fV0_NegTOFBunchCrossing < cutval_V0[kV0_TOFBunchCrossing]) && (fV0_PosTOFBunchCrossing < cutval_V0[kV0_TOFBunchCrossing])) continue;
      IsOOBPileUpBalbino = !((nStatus & AliESDtrack::kITSrefit) || (pStatus & AliESDtrack::kITSrefit) || (fV0_NegTOFBunchCrossing > cutval_V0) || (fV0_PosTOFBunchCrossing > cutval_V0));
      if (!IsOOBPileUpBalbino)        fHistEventV0->Fill(28); 
      fHistTOFBunchCrossing->Fill(fV0_PosTOFBunchCrossing);

      AliPIDResponse::EDetPidStatus statusTOFPos;
      AliPIDResponse::EDetPidStatus statusTOFNeg;
      bool isTOFPIDok = kFALSE;
      statusTOFPos = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,prongTrackPos);
      statusTOFNeg = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,prongTrackNeg);

      //TOF PID (not used)
      if ( (statusTOFPos ==  AliPIDResponse::kDetPidOk) && (statusTOFNeg ==  AliPIDResponse::kDetPidOk) ) {  
	if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF( prongTrackPos, (AliPID::EParticleType)2))< 3.) goodPiPlusTOF=kTRUE; 
	if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF( prongTrackNeg, (AliPID::EParticleType)2))< 3.) goodPiMinusTOF=kTRUE;
	//	isTOFPIDok=kTRUE;
      }

      Bool_t HasPointOnSPDPos=kFALSE;
      Bool_t HasPointOnSPDNeg=kFALSE;
      if (prongTrackPos->HasPointOnITSLayer(0) || prongTrackPos->HasPointOnITSLayer(1) ) HasPointOnSPDPos=kTRUE;
      if (prongTrackNeg->HasPointOnITSLayer(0) || prongTrackNeg->HasPointOnITSLayer(1) ) HasPointOnSPDNeg=kTRUE;

      if ( !(statusTOFPos ==  AliPIDResponse::kDetPidOk) && !(statusTOFNeg ==  AliPIDResponse::kDetPidOk) && !(HasPointOnSPDPos) && !(HasPointOnSPDNeg)) continue; //out of bunch pile up (a V0 daughter track is required to have at least a hit in the TOF or in one of the first two ITS layers (SPD))
      fHistEventV0->Fill(10);

      if(isTOFPIDok){
	if(goodPiPlusTPC && goodPiPlusTOF){
	  goodPiPlus=kTRUE;
	}      
	if(goodPiMinusTPC && goodPiMinusTOF){
	  goodPiMinus=kTRUE;
	}      
      } else {
	if(goodPiPlusTPC){
	  goodPiPlus=kTRUE;
	}      
	if(goodPiMinusTPC){
	  goodPiMinus=kTRUE;
	}
      }


      if(!goodPiMinus || !goodPiPlus )             continue;
      fHistEventV0->Fill(11);

      if(TMath::Abs(v0->EtaProng(pos0or1))>0.8) continue;
      if(TMath::Abs(v0->EtaProng(neg0or1))>0.8) continue;
      fHistEventV0->Fill(12);
    
      //DAUGHTER TRACK QUALITY 
      if(nTPCCrossedRowspos<80) continue;
      if(nTPCCrossedRowsneg<80) continue;
      if(rationCrnFindpos<0.8)  continue;
      if(rationCrnFindneg<0.8)  continue;
      fHistEventV0->Fill(13);    
    
      //      if(TMath::Abs(v0->Eta()) >  fEtaV0Assoc)	   	continue; // to be done offline
      fHistEventV0->Fill(14);    
    
      if(fReadMCTruth){
	if (fMCEvent){
	  //cout << "\n this particle has passed all but pt cuts: let's fill the mass Pt histo for true reco K0s "<< endl;
	  if(PdgPos==211 && PdgNeg==-211 && PdgMotherPos == 310 &&  PdgMotherNeg == 310 && labelMotherPos==labelMotherNeg  && MotherPos->IsPhysicalPrimary()){
	    fHistReconstructedV0PtMass->Fill(v0->MassK0Short(),v0->Pt(), lPercentiles);
	    fHistCommonPartonTrueCasc->Fill(PdgCodeCommon, TrigLabel, V0Label);
	  }
	}
      }
    
      if(v0->CosPointingAngle(lBestPrimaryVtxPos) < kMinCosAngle[ParticleType]) 	continue;
      //    if(v0->DecayLengthV0(lBestPrimaryVtxPos) > kMaxDL)          	continue;
      if(v0->DecayLengthV0(lBestPrimaryVtxPos) < kMinDL[ParticleType])	   	continue;

      double v0Dca = v0->DcaV0ToPrimVertex();
      if(v0->DcaNegToPrimVertex() < kMinDCAPrimary[ParticleType])      continue;
      if(v0->DcaPosToPrimVertex() < kMinDCAPrimary[ParticleType])      continue;
      if(v0->DcaV0Daughters() > kMaxDCADaughters[ParticleType])        continue;
      if(v0Dca > kMaxDCA[ParticleType]) 	        	       continue;
      if(kctau[ParticleType]>kctauval[ParticleType])                   continue;
   
      //if(v0->PtArmV0()< 0.2*TMath::Abs(v0->AlphaV0()))                    continue;
      //fHistPtArmvsAlphaAfterSelection->Fill(v0->AlphaV0(), v0->PtArmV0());    
      fHistPtArmvsAlpha->Fill(v0->AlphaV0(),v0->PtArmV0());       

      if(v0->MassK0Short()< 0.55 && v0->MassK0Short()> 0.45) {
	fHistPtArmvsAlphaAfterSelection->Fill(v0->AlphaV0(),v0->PtArmV0());    
      }

      fHistMassPhoton->Fill(MassPhoton);   
      if(MassPhoton>0.1){
	fHistPtArmvsAlphaAfterPhotonSelection->Fill(v0->AlphaV0(),v0->PtArmV0());    
      }

      if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.005 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.005) {
	fHistPtArmvsAlphaAfterLambdaRejectionSelection->Fill(v0->AlphaV0(),v0->PtArmV0());    
      }

      if(v0->MassK0Short()> 0.55 || v0->MassK0Short()< 0.45) continue;
      fHistEventV0->Fill(15);    
      // if(TMath::Abs((v0->MassLambda() - massLambda))< 0.005) continue;
      // if(TMath::Abs((v0->MassAntiLambda() - massLambda))< 0.005) continue;
      fHistEventV0->Fill(16);    

      //V0 daugthter DCA cut (just a check)
      Float_t DCAxyPos=-999.; Float_t DCAxyNeg=-999.;
      Float_t DCAzPos=-999.; Float_t DCAzNeg=-999.;
      Float_t DCAxyPosGlobal=-999.; Float_t DCAxyNegGlobal=-999.;
      Float_t DCAzPosGlobal=-999.; Float_t DCAzNegGlobal=-999.;

      Bool_t  DCAxyDaughters=kFALSE;
      //      cout << " prongtrackneg id " << prongTrackNeg->GetID() << " prong track pos id " << prongTrackPos->GetID() << endl;

      //      cout << " \ngetting track parameters for positive and negative V0 daughter " << endl;
      prongTrackPos->GetImpactParameters(&DCAxyPos, &DCAzPos);
      prongTrackNeg->GetImpactParameters(&DCAxyNeg, &DCAzNeg);
      //      cout << " DCAxyPos: "  << DCAxyPos <<  " DCAxyNeg: "  << DCAxyNeg <<  " DCAzPos: "  << DCAzPos <<  " DCAzNeg: "  << DCAzNeg <<endl;
      
      //      cout << "\n******* getting track parameters for positive V0 daughter (globaltrack)" << endl;
      etpPos.CopyFromVTrack(prongTrackPos); 
      etpPos.PropagateToDCA(vertexmain,(InputEvent())->GetMagneticField(),100.,dzgPos,covargPos);    
      DCAxyPosGlobal = dzgPos[0];
      DCAzPosGlobal = dzgPos[1];
      //      cout << " from prongtrackneg: DCAxyPosGlobal: "  << DCAxyPosGlobal <<  " DCAxyNegGlobal: "  << DCAxyNegGlobal <<  " DCAzPosGlobal: "  << DCAzPosGlobal <<  " DCAzNegGlobal: "  << DCAzNegGlobal <<endl;

      //      cout << "\n****** getting track parameters for negative V0 daughter (globaltrack)" << endl;
      etpNeg.CopyFromVTrack(prongTrackNeg); 
      etpNeg.PropagateToDCA(vertexmain,(InputEvent())->GetMagneticField(),100.,dzgNeg,covargNeg);    
      DCAxyNegGlobal = dzgNeg[0];
      DCAzNegGlobal = dzgNeg[1];
      //      cout << " from prongtrackneg: DCAxyPosGlobal: "  << DCAxyPosGlobal <<  " DCAxyNegGlobal: "  << DCAxyNegGlobal <<  " DCAzPosGlobal: "  << DCAzPosGlobal <<  " DCAzNegGlobal: "  << DCAzNegGlobal <<endl;

      /*
      Float_t fPositionPos[2]={-999.};
      prongTrackPos->GetPosition(fPositionPos);
      cout << "z position pos " << fPositionPos[2] << " z position primary vertex " << lBestPrimaryVtxPos[2] << endl;
      cout << "y position pos " << fPositionPos[1] << " y position primary vertex " << lBestPrimaryVtxPos[1] << endl;
      cout << "x position pos " << fPositionPos[0] << " x position primary vertex " << lBestPrimaryVtxPos[0] << endl;
      */
      if (TMath::Abs(DCAxyPosGlobal)<=2.4 && TMath::Abs(DCAxyNegGlobal)<=2.4 ) {
	//	  fHistEventV0->Fill(17);
	  DCAxyDaughters=kTRUE;
      }
      if (TMath::Abs(DCAzPosGlobal)<=3.2 && TMath::Abs(DCAzNegGlobal)<=3.2 && DCAxyDaughters ) {
	fHistEventV0->Fill(17);
      }
      if(prongTrackPos->TestFilterBit(1) && prongTrackNeg->TestFilterBit(1)) 	fHistEventV0->Fill(18);

      bool skipV0=kFALSE;
      if(fReadMCTruth){
	if (fMCEvent){
	  if(PdgPos==211 && PdgNeg==-211 && PdgMotherPos == 310 &&  PdgMotherNeg == 310 && labelMotherPos==labelMotherNeg){
	    if(MotherPos->IsPhysicalPrimary()){
	      fHistSelectedV0PtMass->Fill(fTreeVariableInvMassK0s,fTreeVariablePtV0, lPercentiles);
	    }	    
	  }
	}
      }
    
    
      if(!(v0->Pt()> fminPtV0 && v0->Pt()<fmaxPtV0) )continue;
      fHistEventV0->Fill(19);     
      if (v0->Pt()>=ptTriggerMassimoDati){
	skipV0=kTRUE;
	NumberSecondParticleNoAssoc++;
      }
      //      if (skipV0)	continue;
     

      fMassV0->Fill(v0->MassK0Short());    
      fDCAxyDaughters->Fill(DCAxyPosGlobal, DCAxyNegGlobal);
      fDCAzDaughters->Fill(DCAzPosGlobal, DCAzNegGlobal);
      fDCAxyDaughtersBis->Fill(DCAxyPos, DCAxyNeg);
      fDCAzDaughtersBis->Fill(DCAzPos, DCAzNeg);
      fDCAxyPos->Fill(DCAxyPosGlobal);
      fDCAxyNeg->Fill(DCAxyNegGlobal);
      fDCAzPos->Fill(DCAzPosGlobal);
      fDCAzNeg->Fill(DCAzNegGlobal);

      for (Int_t m =0; m<5;m++){
	if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	  fHistMassvsPt_tagli[m]->Fill(v0->MassK0Short(),v0->Pt());      
	}
      }
      fHistMassvsPt_tagli[5]->Fill(v0->MassK0Short(),v0->Pt());      

      NumberSecondParticle++;

      labelPrimOrSecV0=0;

      if(fReadMCTruth){
	if (fMCEvent){
	  V0PDGCode=0;
	  if(PdgPos==211 && PdgNeg==-211 && PdgMotherPos == 310 &&  PdgMotherNeg == 310 && labelMotherPos==labelMotherNeg){
	    V0PDGCode=PdgMotherPos;
	    if(MotherPos->IsPhysicalPrimary()){

	      isaK0s=1;
	      //c cout <<"selected v0 Pt " <<  v0->Pt() << endl;
	      /* 3D histos
		 fHistResolutionV0Pt->Fill(v0->Pt()- MotherPos->Pt(), lPercentiles, ptTriggerMassimoDati);
		 fHistResolutionV0Phi->Fill(v0->Phi()- MotherPos->Phi(), lPercentiles, ptTriggerMassimoDati);
		 fHistResolutionV0PhivsPt->Fill(v0->Phi()- MotherPos->Phi(), lPercentiles, v0->Pt());
		 fHistResolutionV0PtvsPt->Fill(v0->Pt()- MotherPos->Pt(), lPercentiles, v0->Pt());
		 fHistResolutionV0Eta->Fill(v0->Eta()- MotherPos->Eta(), lPercentiles, ptTriggerMassimoDati);
	      */
	      //2D histos
	      fHistAssocPtRecovsPtGen->Fill(MotherPos->Pt(), v0->Pt());
	      if(!skipV0){
	      fHistResolutionV0Pt->Fill(v0->Pt()- MotherPos->Pt(), ptTriggerMassimoDati);
	      fHistResolutionV0Phi->Fill(v0->Phi()- MotherPos->Phi(), ptTriggerMassimoDati);
	      fHistResolutionV0PhivsPt->Fill(v0->Phi()- MotherPos->Phi() , v0->Pt());
	      fHistResolutionV0PtvsPt->Fill(v0->Pt()- MotherPos->Pt(), v0->Pt());
	      fHistResolutionV0Eta->Fill(v0->Eta()- MotherPos->Eta(),  ptTriggerMassimoDati);
	      //
	      }

	      if (TMath::Abs(rapidityV0)< 0.5){
	      fHistSelectedV0PtTMaxPhi[0]->Fill(ptTriggerMassimoDati, v0->Phi(), lPercentiles);
	      fHistSelectedV0PtTMaxEta[0]->Fill(ptTriggerMassimoDati, v0->Eta(), lPercentiles);
	      fHistSelectedV0PtPtTMax[0]->Fill(v0->Pt(),ptTriggerMassimoDati , lPercentiles);
	      fHistSelectedV0PtEta[0]->Fill(v0->Pt(), v0->Eta() , lPercentiles);
	      fHistSelectedGenV0PtPtTMax[0]->Fill(MotherPos->Pt(),ptTriggerMassimoDati , lPercentiles);
	      }

	      labelPrimOrSecV0=1;
	      NumberSecondParticleRecoTrue++;
	    }
	    else if(MotherPos->IsSecondaryFromWeakDecay())      labelPrimOrSecV0=2;
	    else if(MotherPos->IsSecondaryFromMaterial())      labelPrimOrSecV0=3;
	    else labelPrimOrSecV0=4;
	    
	    if(!MotherPos->IsPhysicalPrimary())      fHistAssocPtRecovsPtGenNotPrim->Fill(MotherPos->Pt(), v0->Pt());

	    for (Int_t m =0; m<5;m++){
	      if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
		for(Int_t p=1; p<=4; p++){
		  if (labelPrimOrSecV0==p) {
		    if (TMath::Abs(rapidityV0)< 0.5)		    fHistPrimaryV0[m][0]->Fill(p, v0->Pt(), ptTriggerMassimoDati);   
		  }
		}
	      }
	    }
	    for(Int_t p=1; p<=4; p++){
	      if (labelPrimOrSecV0==p){
		if (TMath::Abs(rapidityV0)< 0.5)		fHistPrimaryV0[5][0]->Fill(p, v0->Pt(), ptTriggerMassimoDati);
	      }
	    }
	  }
	}
      }

      //save second particle information (V0)
      //cout << "save second particle information (V0) "<< endl;
      if((fReadMCTruth && isEfficiency) || (!fReadMCTruth)){
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sLabelMotherPos = labelMotherPos;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sLabelPos       = labelPos;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sLabelNeg       = labelNeg;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sPDGcode        = V0PDGCode;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].isP             = labelPrimOrSecV0;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sDcaPosV0       = v0->DcaPosToPrimVertex();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sDcaNegV0       = v0->DcaNegToPrimVertex();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sPtArmV0        = v0->PtArmV0();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sAlphaV0        = v0->AlphaV0();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sInvMassK0s     = v0->MassK0Short();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sInvMassLambda   = v0->MassLambda();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sInvMassAntiLambda = v0->MassAntiLambda();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sCosPointingAngle  = v0->CosPointingAngle(lBestPrimaryVtxPos);
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sDcaV0ToPV    = v0->DcaV0ToPrimVertex();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sRap          = rapidityV0;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sPt           = v0->Pt();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sctau         = kctau[ParticleType];
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sEta          = v0->Eta();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sTheta        = v0->Theta();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sPhi          = v0->Phi();
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sCharge       = 0;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sDCAz         = 0;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sDCAxy        = 0;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sAssocOrNot   = skipV0;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sIsCommonParton        = IsCommonParton;
	fEvt->fReconstructedSecond[NumberSecondParticle-1].sPdgCommonParton       = PdgCodeCommon;


      
	fHistPtV0->Fill(v0->Pt());
      }
    } //end loop for K0s particles as associated
    } // done only if not MC truth  

    cout << " beginning loop for associated particles (MCtruth) " << endl;
    //begin MC truth loop for K0s particles as associated 
    if(fReadMCTruth){
      if (fMCEvent){
	cout << " MC truth for associated particles " << endl;
	AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArray == NULL){
	  return;
	  Printf("ERROR: stack not available");
	}
	for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
	  Bool_t skipV0_MC=kFALSE;
	  AliAODMCParticle* particleV0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
	  if (!particleV0) continue;
	  if((particleV0->GetPdgCode())!=310) continue;
	  //	  if(TMath::Abs(particleV0->Eta())> fEtaV0Assoc)continue;
	  if (!(particleV0->IsPhysicalPrimary()))continue; 
	  if(!(particleV0->Pt()> fminPtV0 && particleV0->Pt()<fmaxPtV0) )continue;

	  //question: should I take all daughter tracks or only those with |eta| < 0.8 as in the data?
	  AliAODMCParticle* particlePos;
	  AliAODMCParticle* particleNeg;
	  
	  Int_t labelPos = particleV0->GetDaughterLabel(0);
	  Int_t labelNeg = particleV0->GetDaughterLabel(1);
	  particlePos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelPos)));
	  particleNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelNeg)));
	  
	  //	  labelPos =  (particleV0->GetDaughterFirst());
	  //	  labelNeg =  (particleV0->GetDaughterLast());
	  //cout << "label pos " << labelPos << " labelNeg " << labelNeg << " q+ " << particlePos->Charge() << " q- " << particleNeg->Charge() << endl;
	  //cout << " PDG+ " << particlePos->GetPdgCode() << " PDG- " << particleNeg->GetPdgCode() << endl;
	  //cout << "eta + " << particlePos->Eta() << endl;
	  //cout << "eta - " << particleNeg->Eta() << endl;
	  //	  EtaPos = particlePos->Eta();
	  //	  EtaNeg = particleNeg->Eta();

	  if (!isHybridMCTruth){
	    if ((particleV0->Pt())>=ptTriggerMassimoMC) {
	      skipV0_MC=kTRUE;
	      NumberSecondParticleMCNoAssoc++;
	    }
	  }
	  else {
	    if ((particleV0->Pt())>=ptTriggerMassimoDati) {
	      skipV0_MC=kTRUE;
	      NumberSecondParticleMCNoAssoc++;
	    }
	  }
	  //	  if (skipV0_MC)      continue;

	  NumberSecondParticleMC++;
	  if(isEfficiency) continue;

	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPDGcode      = particleV0->GetPdgCode();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].isP           = 1;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDcaPosV0     = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDcaNegV0     = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPtArmV0      = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sAlphaV0      = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sInvMassK0s   = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sInvMassLambda   = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sInvMassAntiLambda = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sCosPointingAngle  = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDcaV0ToPV    = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sRap          = particleV0->Y();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPt           = particleV0->Pt();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sctau         = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sEta          = particleV0->Eta();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sTheta        = particleV0->Theta();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPhi          = particleV0->Phi();
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sCharge       = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDCAz         = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDCAxy        = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sAssocOrNot   = skipV0_MC;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sIsCommonParton        = 0;
	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPdgCommonParton       = 0;
	  //	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sEtaPos       = particlePos->Eta();
	  //	  fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sEtaNeg       = particleNeg->Eta();


	}
	//	  cout << "I found " << 	NumberSecondParticleMC - 	    NumberSecondParticleMCNoAssoc<< " associated particles (pT < pt,Trig)" << endl;
	//	  cout << "I found " << 	NumberSecondParticleMC << " associated particles (all pt)" << endl;
	  
      }
    } //end MC truth loop for K0s particles as associated 
    //    cout << " end of loop for assoc particles (MCtruth)" << endl;
  }
  //************************end hV0 section*************************************************

  if (!fIshhCorr){
    fHistEventV0->AddBinContent(20, NumberSecondParticle);    
    fHistEventV0->AddBinContent(21, NumberSecondParticleRecoTrue);    
    fHistEventV0->AddBinContent(22, NumberSecondParticleMC);  
    fHistEventV0->AddBinContent(23, NumberSecondParticleNoAssoc);    
    fHistEventV0->AddBinContent(24, NumberSecondParticleMCNoAssoc);      
  }
  if (fIshhCorr){
    fHistTrackAssoc->AddBinContent(13, NumberSecondParticle);    
    fHistTrackAssoc->AddBinContent(14, NumberSecondParticleRecoTrue);    
    fHistTrackAssoc->AddBinContent(15, NumberSecondParticleMC);  
    fHistTrackAssoc->AddBinContent(16, NumberSecondParticleNoAssoc);    
    fHistTrackAssoc->AddBinContent(17, NumberSecondParticleMCNoAssoc);      
  }
  fHistMultvsV0All->Fill(NumberSecondParticle,lPercentiles);
  fHistMultvsV0AllTruth->Fill(NumberSecondParticleRecoTrue,lPercentiles);
  fHistMultvsV0MCAll->Fill(NumberSecondParticleMC,lPercentiles);
  fHistSecondParticleAll->Fill(NumberSecondParticle, NumberSecondParticleMC);
  fHistSecondParticleTruthAll->Fill(NumberSecondParticleRecoTrue, NumberSecondParticleMC);
  if(NumberSecondParticleMC==0 && NumberSecondParticle!=0) fHistEventMult->Fill(17); 
  if(NumberSecondParticleMC < NumberSecondParticleRecoTrue) fHistEventMult->Fill(18); 
  fHistTriggerNotLeading->Fill(NumberSecondParticleNoAssoc,lPercentiles, ptTriggerMassimoDati);
  fHistTriggerNotLeadingMC->Fill(NumberSecondParticleMCNoAssoc,lPercentiles, ptTriggerMassimoMC);

  Bool_t DoubleCounted=kFALSE;
  Bool_t CommonDaughtPos=kFALSE;
  Bool_t CommonDaughtNeg=kFALSE;
  if(!fIshhCorr){
    if((fReadMCTruth && isEfficiency) || (!fReadMCTruth)){
      for(Int_t j=0; j < NumberSecondParticle; j++){
	if(  fEvt->fReconstructedSecond[j].isP           ==1){
	  // fHistEventMult->Fill(20);
	  for(Int_t i=0; i< NumberSecondParticle; i++){
	    if(  fEvt->fReconstructedSecond[i].isP           ==1  && i!=j){
	      if(    fEvt->fReconstructedSecond[j].sLabelMotherPos ==     fEvt->fReconstructedSecond[i].sLabelMotherPos) {
		DoubleCounted=kTRUE;

		if(    fEvt->fReconstructedSecond[j].sLabelPos ==     fEvt->fReconstructedSecond[i].sLabelPos) CommonDaughtPos=kTRUE;
		if(    fEvt->fReconstructedSecond[j].sLabelNeg ==     fEvt->fReconstructedSecond[i].sLabelNeg) CommonDaughtNeg=kTRUE;
	      }
	    }
	  }
	}
      }
    }
    if (DoubleCounted)  fHistEventMult->Fill(19);
    if (CommonDaughtPos)  fHistEventMult->Fill(20);
    if (CommonDaughtNeg)  fHistEventMult->Fill(21);
  }

  if((fReadMCTruth && isEfficiency) || (!fReadMCTruth)) NumberSecondParticleAll=NumberSecondParticle;
  else if (fReadMCTruth && !isEfficiency) NumberSecondParticleAll=NumberSecondParticleMC;

  if(NumberSecondParticleAll==0){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3,fBkgTree);
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);
    //    cout << "event has trigger particle but no V0 " << endl;     
    return;
  }

  fHistSecondParticle->Fill(NumberSecondParticle, NumberSecondParticleMC);
  fHistSecondParticleTruth->Fill(NumberSecondParticleRecoTrue, NumberSecondParticleMC);

  if(!fIshhCorr){
    fHistEventV0->AddBinContent(25, NumberSecondParticle);    
    fHistEventV0->AddBinContent(26, NumberSecondParticleRecoTrue);    
    fHistEventV0->AddBinContent(27, NumberSecondParticleMC);    
  }
  if(fIshhCorr){
    fHistTrackAssoc->AddBinContent(18, NumberSecondParticle);    
    fHistTrackAssoc->AddBinContent(19, NumberSecondParticleRecoTrue);    
    fHistTrackAssoc->AddBinContent(20, NumberSecondParticleMC);    
  }
  // FifoShiftok=kTRUE;
  // LastzBin=zBin;
  // LastcentralityBin=centralityBin;
  
  fHistEventMult->Fill(22);  
  if(!fReadMCTruth || (fReadMCTruth && isEfficiency)){
    fEvt->fNumberCandidateFirst = NumberFirstParticle;
    fEvt->fNumberCandidateSecond = NumberSecondParticle;
  }
  else if (fReadMCTruth && !isEfficiency && !isHybridMCTruth) {
    fEvt->fNumberCandidateFirst = NumberFirstParticleMC;
    fEvt->fNumberCandidateSecond = NumberSecondParticleMC;
  }
  else { //for hybrid
    fEvt->fNumberCandidateFirst = NumberFirstParticle;
    fEvt->fNumberCandidateSecond = NumberSecondParticleMC;
  }

  //Fill histos for selected events (events with which I perform the same-event and the mixed-event angular correlation, that is having NT> 0 and NV0>0
  fHistTriggerwV0->Fill(NumberFirstParticle);
  fHistTriggerwV0MCTruth->Fill(NumberFirstParticleMC);
  fHist_multiplicity->Fill(lPercentiles);
  fHistZvertex->Fill(lBestPrimaryVtxPos[2]);
  fHistTrack->AddBinContent(16, NumberFirstParticle);
  fHistTrack->AddBinContent(17, NumberFirstParticleMC);
  fHistMultvsTrigger->Fill(NumberFirstParticle, lPercentiles);
  fHistMultvsTriggerMCTruth->Fill(NumberFirstParticleMC, lPercentiles);
  fHistMultiplicityVsVertexZ->Fill(lBestPrimaryVtxPos[2], lPercentiles);
 
  for(Int_t i=0; i< 100; i++){
    if(( lPercentiles <i+1) && (lPercentiles >= i) ) fHistTriggervsMult->AddBinContent(i+1,NumberFirstParticle);  
  }
  for(Int_t i=0; i< 100; i++){
    if(( lPercentiles <i+1) && (lPercentiles >= i) ) fHistTriggervsMultMC->AddBinContent(i+1,NumberFirstParticleMC);  
  }


  for(Int_t l=0; l< NumberFirstParticleAll; l++){
    fHistPt->Fill(fEvt->fReconstructedFirst[l].fPt);
    fHistPtvsMult->Fill(fEvt->fReconstructedFirst[l].fPt, lPercentiles);
    fHist_eta_phi->Fill(fEvt->fReconstructedFirst[l].fPhi, fEvt->fReconstructedFirst[l].fEta);
 
  }

  if (!fReadMCTruth || (fReadMCTruth && isEfficiency) || (isHybridMCTruth)){
    fHistPtMaxvsMult->Fill(ptTriggerMassimoDati, lPercentiles);
    fHist_eta_phi_PtMax->Fill(phiTriggerMassimoDati, etaTriggerMassimoDati);
  }
  if (fReadMCTruth && !isEfficiency && !isHybridMCTruth) {
    fHistPtMaxvsMult->Fill(ptTriggerMassimoMC, lPercentiles);
    fHist_eta_phi_PtMax->Fill(phiTriggerMassimoMC, etaTriggerMassimoMC);
  }

  /* this part is implemented above in a different way 
    if(fFirstpart == fSecondpart){ 
    DoPairshh(lPercentiles, fieldsign);  
    else{
    //Remove candidates that are at the same time a ptc1 and ptc2
    for (int i=0; i < fEvt->fNumberCandidateFirst; i++) {
    for (int j=0; j<fEvt->fNumberCandidateSecond; j++) {
    if (fEvt->fReconstructedFirst[i].index == fEvt->fReconstructedSecond[j].index) {
    //cout<<"the track can be both tracks!"<<endl;
    fEvt->fReconstructedFirst[i].doSkipOver = kTRUE;
    fEvt->fReconstructedSecond[j].doSkipOver = kTRUE;
    }
    }
    }
  */
  
  //--------------------------------------------------------------
  Double_t ptTriggerMassimo=0;
  if (!fReadMCTruth || (fReadMCTruth && isEfficiency) || (fReadMCTruth && isHybridMCTruth) ) ptTriggerMassimo=ptTriggerMassimoDati;
  if (fReadMCTruth && !isEfficiency && !isHybridMCTruth) ptTriggerMassimo=ptTriggerMassimoMC;

  /*
  cout << "\npttrigger massimo dati (per hybrid, data, MC reco) " << ptTriggerMassimoDati << endl;
  cout << "pttrigger massimo mc " << ptTriggerMassimoMC << endl;
  cout << " number candidates first " << fEvt->fNumberCandidateFirst << endl;
  cout << " number first particle " <<   NumberFirstParticleAll << " mc " << NumberFirstParticleMC << " data " << NumberFirstParticle << endl;
  cout << " number candidates second " << fEvt->fNumberCandidateSecond << endl;
  cout << " number first particle " <<   NumberSecondParticleAll << " mc " << NumberSecondParticleMC << " data " << NumberSecondParticle << endl;
  */

  DoPairsh1h2((Int_t)lPercentiles, fieldsign, lBestPrimaryVtxPos[2], ptTriggerMassimo);  

  PostData(1, fOutputList);     
  PostData(2, fSignalTree);
  PostData(3, fBkgTree);
  PostData(4, fOutputList2);  
  PostData(5, fOutputList3);     
  PostData(6, fOutputList4);     

  fEventColl[zBin][centralityBin]->FifoShift();       
}


//----------------------------------------------------------------------------------------------------

void AliAnalysisTaskCorrelationhhK0s::DoPairsh1h2 ( const Float_t lPercentiles, int fieldsign, Double_t lBestPrimaryVtxPos, Double_t ptTriggerMassimo)  {

  //-----------
  double DCAxyP1  = -999. ;
  double DCAzP1   = -999. ;  
  double DCAxyP2  = -999. ; 
  double DCAzP2   = -999. ;  

  double  ptP1 = -999.;
  double  ptP2 = -999.;

  bool isP1  = kFALSE;
  bool isaP1 = kFALSE;
  bool isP2  = kFALSE;
  bool isaP2 = kFALSE;

  Int_t  SignP1 = -999;
  Int_t  SignP2 = -999;

  double dphi  = -999.;
  double dphis = -999.;
  double deta  = -999.;
  double dtheta = -999.;
  
  bool isMC1 = kFALSE;
  bool isMC2 = kFALSE;

  bool isMCvector = kFALSE;
  bool sameMother = kFALSE;
  bool sameGrandMother = kFALSE;
  
  Int_t mcMotherLabelP1 = -999;
  Int_t mcMotherLabelP2 = -999;

  Int_t mcGrandMotherLabelP1 = -999;
  Int_t mcGrandMotherLabelP2 = -999;
 
  Int_t typeP1 = -999;
  Int_t typeP2 = -999;
 
  Int_t mcPDGMotherP1 = 0;
  Int_t mcPDGMotherP2 = 0;
  //  Int_t mcMotherBin = 0;
 
  Int_t mcPDGGrandMother = 0;
  Int_t mcGrandMotherBin = 0;

  Int_t mcPDGcodeP1 = 0;
  Int_t mcPDGcodeP2 = 0;

  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
 
  double pairKstar   = 0.;
  double pairKstarMC = 0.;
  double pairMass  = 0.;
  double pairMassE = 0.;
  double pairKt    = 0.;
  double pTmax = 0.;

  //  cout << " I'm doing the trigger-assoc association " << endl;
  for (int i=0; i<fEvt->fNumberCandidateFirst; i++) {
    //    cout << "pt " << fEvt->fReconstructedFirst[i].fPt << endl;
    //I select as trigger particle only the highest-pT one
    if ( fEvt->fReconstructedFirst[i].fPt <ptTriggerMassimo ) continue;
    //    cout << " I got the trigger particle ! "<< endl;
    for (int eventNumber=0; eventNumber<fnEventsToMix+1; eventNumber++) { 
      //if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateFirst)!=0.) evmultmixed++; 
      if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateSecond)!=0.) {
	evmultmixed++; 
      }
      Int_t count=0;
      for (int j=0; j<(fEvt+eventNumber)->fNumberCandidateSecond; j++) {
	//	if (j>0) continue; //new and wrong!!	
	//c 	if ((fEvt+eventNumber)->fReconstructedSecond[j].doSkipOver) continue;
	//	cout << "pt of the associated " << (fEvt+eventNumber)->fReconstructedFirst[j].fPt<< endl;
	deta   = CalculateDeltaEta(fEvt->fReconstructedFirst[i].fEta, (fEvt+eventNumber)->fReconstructedSecond[j].sEta);
	dtheta = CalculateDeltaTheta(fEvt->fReconstructedFirst[i].fTheta, (fEvt+eventNumber)->fReconstructedSecond[j].sTheta);
	//dphi   = CalculateDeltaPhi(fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi);
	dphi = CalculateDPhiStar(fEvt->fReconstructedFirst[i].fCharge, (fEvt+eventNumber)->fReconstructedSecond[j].sCharge, fieldsign,fEvt->fReconstructedFirst[i].fPt , (fEvt+eventNumber)->fReconstructedSecond[j].sPt, fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi,0); 
	
	if (eventNumber==0) {//Same event pair histogramming
	  count++;
	  //	  cout << " counter..."  << endl;
	  fTreeVariablePtTrigger              = fEvt->fReconstructedFirst[i].fPt;		       
	  fTreeVariableChargeTrigger          = fEvt->fReconstructedFirst[i].fCharge;		       
	  fTreeVariableEtaTrigger             = fEvt->fReconstructedFirst[i].fEta; 		       
	  fTreeVariablePhiTrigger	      =	fEvt->fReconstructedFirst[i].fPhi;             
	  fTreeVariableDCAz		      =	fEvt->fReconstructedFirst[i].fDCAz;             
	  fTreeVariableDCAxy		      =	fEvt->fReconstructedFirst[i].fDCAxy;   
	  fTreeVariablePDGCodeTrigger         = fEvt->fReconstructedFirst[i].fPDGcode ;
	  fTreeVariableChargeAssoc            = fEvt->fReconstructedSecond[j].sCharge;		       
	  fTreeVariableAssocDCAz	      =	fEvt->fReconstructedSecond[j].sDCAz;             
	  fTreeVariableAssocDCAxy	      =	fEvt->fReconstructedSecond[j].sDCAxy;              
	  fTreeVariableRapK0Short	      =	fEvt->fReconstructedSecond[j].sRap;     	      
	  fTreeVariableDcaV0ToPrimVertex      = fEvt->fReconstructedSecond[j].sDcaV0ToPV;	       	      
	  fTreeVariableDcaPosToPrimVertex     = fEvt->fReconstructedSecond[j].sDcaPosV0; 	       	      
	  fTreeVariableDcaNegToPrimVertex     = fEvt->fReconstructedSecond[j].sDcaNegV0;  	       	      
	  fTreeVariableV0CosineOfPointingAngle= fEvt->fReconstructedSecond[j].sCosPointingAngle;      
	  fTreeVariablePtV0		      = fEvt->fReconstructedSecond[j].sPt;     
	  fTreeVariablectau		      = fEvt->fReconstructedSecond[j].sctau;     
	  fTreeVariableInvMassK0s	      =	fEvt->fReconstructedSecond[j].sInvMassK0s;   
	  fTreeVariableInvMassLambda	      =	fEvt->fReconstructedSecond[j].sInvMassLambda;   
	  fTreeVariableInvMassAntiLambda      =	fEvt->fReconstructedSecond[j].sInvMassAntiLambda;   
	  fTreeVariableEtaV0		      =	fEvt->fReconstructedSecond[j].sEta;    
	  fTreeVariablePhiV0		      =	fEvt->fReconstructedSecond[j].sPhi;   
	  fTreeVariablePtArmenteros           = fEvt->fReconstructedSecond[j].sPtArmV0;     
	  fTreeVariableAlpha	              = fEvt->fReconstructedSecond[j].sAlphaV0;  
	  fTreeVariableDeltaEta	      	      = deta;  
	  fTreeVariableDeltaPhi		      = dphi;
	  fTreeVariableDeltaTheta             = dtheta;
	  fTreeVariableMultiplicity	      = lPercentiles;
	  fTreeVariableZvertex                = lBestPrimaryVtxPos;
	  fTreeVariablePDGCodeAssoc           = fEvt->fReconstructedSecond[j].sPDGcode;
	  fTreeVariableisPrimaryTrigger       =  fEvt->fReconstructedFirst[i].isP ;
	  fTreeVariableisPrimaryV0            =  fEvt->fReconstructedSecond[j].isP ;
	  fTreeVariableSkipAssoc              =  fEvt->fReconstructedSecond[j].sAssocOrNot ;
	  fTreeVariableIsCommonParton         =  fEvt->fReconstructedSecond[j].sIsCommonParton       ;
          fTreeVariablePdgCommonParton        =  fEvt->fReconstructedSecond[j].sPdgCommonParton      ;
	  //	  fTreeVariableEtaDaughterPos         =  fEvt->fReconstructedSecond[j].sEtaPos ;
	  //	  fTreeVariableEtaDaughterNeg         =  fEvt->fReconstructedSecond[j].sEtaNeg ;

	  fSignalTree->Fill();  

	}

	else {//Mixed-event pair histogramming
	  fTreeVariablePtTrigger              = fEvt->fReconstructedFirst[i].fPt;		       
	  fTreeVariableChargeTrigger          = fEvt->fReconstructedFirst[i].fCharge;		       
	  fTreeVariableEtaTrigger             = fEvt->fReconstructedFirst[i].fEta; 		       
	  fTreeVariablePhiTrigger	      =	fEvt->fReconstructedFirst[i].fPhi;             
	  fTreeVariableDCAz		      =	fEvt->fReconstructedFirst[i].fDCAz;             
	  fTreeVariableDCAxy		      =	fEvt->fReconstructedFirst[i].fDCAxy;   
	  fTreeVariablePDGCodeTrigger         = fEvt->fReconstructedFirst[i].fPDGcode ;
	  fTreeVariableisPrimaryTrigger       = fEvt->fReconstructedFirst[i].isP;

	  fTreeVariableChargeAssoc            = (fEvt+eventNumber)->fReconstructedSecond[j].sCharge;		       
	  fTreeVariableAssocDCAz	      =	(fEvt+eventNumber)->fReconstructedSecond[j].sDCAz;             
	  fTreeVariableAssocDCAxy	      =	(fEvt+eventNumber)->fReconstructedSecond[j].sDCAxy;                         
	  fTreeVariableRapK0Short	      =	(fEvt+eventNumber)->fReconstructedSecond[j].sRap;     	      
	  fTreeVariableDcaV0ToPrimVertex      = (fEvt+eventNumber)->fReconstructedSecond[j].sDcaV0ToPV;	       	      
	  fTreeVariableDcaPosToPrimVertex     = (fEvt+eventNumber)->fReconstructedSecond[j].sDcaPosV0; 	       	      
	  fTreeVariableDcaNegToPrimVertex     = (fEvt+eventNumber)->fReconstructedSecond[j].sDcaNegV0;  	       	      
	  fTreeVariableV0CosineOfPointingAngle= (fEvt+eventNumber)->fReconstructedSecond[j].sCosPointingAngle;      
	  fTreeVariablePtV0		      = (fEvt+eventNumber)->fReconstructedSecond[j].sPt;     
	  fTreeVariablectau		      = (fEvt+eventNumber)->fReconstructedSecond[j].sctau;     
	  fTreeVariableInvMassK0s	      =	(fEvt+eventNumber)->fReconstructedSecond[j].sInvMassK0s;   
	  fTreeVariableInvMassAntiLambda      =	(fEvt+eventNumber)->fReconstructedSecond[j].sInvMassAntiLambda;   
	  fTreeVariableInvMassLambda	      =	(fEvt+eventNumber)->fReconstructedSecond[j].sInvMassLambda;   
	  fTreeVariableEtaV0		      =	(fEvt+eventNumber)->fReconstructedSecond[j].sEta;    
	  fTreeVariablePhiV0		      =	(fEvt+eventNumber)->fReconstructedSecond[j].sPhi;   
	  fTreeVariablePtArmenteros           = (fEvt+eventNumber)->fReconstructedSecond[j].sPtArmV0;     
	  fTreeVariableAlpha	              = (fEvt+eventNumber)->fReconstructedSecond[j].sAlphaV0;  
	  fTreeVariablePDGCodeAssoc           = (fEvt+eventNumber)->fReconstructedSecond[j].sPDGcode;
	  fTreeVariableisPrimaryV0            =  (fEvt+eventNumber)->fReconstructedSecond[j].isP;
	  fTreeVariableSkipAssoc              =  (fEvt+eventNumber)->fReconstructedSecond[j].sAssocOrNot ;
	  //	  fTreeVariableEtaDaughterPos         =  (fEvt+eventNumber)->fReconstructedSecond[j].sEtaPos ;
	  //	  fTreeVariableEtaDaughterNeg         =  (fEvt+eventNumber)->fReconstructedSecond[j].sEtaNeg ;
	  fTreeVariableIsCommonParton         =    (fEvt+eventNumber)->fReconstructedSecond[j].sIsCommonParton       ;
          fTreeVariablePdgCommonParton        =    (fEvt+eventNumber)->fReconstructedSecond[j].sPdgCommonParton      ;

	  fTreeVariableDeltaEta	       	      =deta;  
	  fTreeVariableDeltaPhi		      =dphi;
	  fTreeVariableDeltaTheta             =dtheta;      
	  fTreeVariableMultiplicity	      = lPercentiles;
	  fTreeVariableZvertex                = lBestPrimaryVtxPos;

       
	  fBkgTree->Fill();  
	  
	} //mixed
	
      } //end second particle loop

    }//end event loop

    if (evmultmixed!=0) multmixedcounted = kTRUE;
    
  } // end first particle loop
  
  //  cout << " I've done the trigger-assoc association " << endl;
  if  (multmixedcounted) 
    fHistMultiplicityOfMixedEvent->Fill(evmultmixed, lPercentiles); //tells me with how many events the mixed-event was done
  
}

/*
//----------------------------------------------------------------------------------------------
//void AliAnalysisTaskKPFemto::DoPairshh (const Float_t lPercentiles, int fieldsign) {
void AliAnalysisTaskKPFemto::DoPairshh (const Int_t lPercentiles, int fieldsign) {
return;
}
*/


//-----------------------------------------------------------------------------------------------
double AliAnalysisTaskCorrelationhhK0s::CalculateDPhiStar(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2,Double_t rad) { //AliFemtoUser/AliFemtoPairCutDetaDphi.h

  const Double_t unit_factor = 0.299792458 / 2.0;
  const Double_t b_field = 0.5006670488586 * magSign;

  Double_t  shift1 = TMath::ASin(unit_factor * chg1 * b_field * rad / ptv1);
  Double_t  shift2 = TMath::ASin(unit_factor * chg2 * b_field * rad / ptv2);

  double dps = (phi1 + shift1) - (phi2 + shift2);
  
  //  dps = TVector2::Phi_mpi_pi(dps); //to be checked

  return dps; //deltaphi;
  
}

//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskCorrelationhhK0s::CalculateDeltaEta( Double_t eta1, Double_t eta2 )  {   //AliFemtoUser/AliFemtoPairCutDetaDphi.h
  const double deta = eta2 - eta1;
  return deta;
}
//_______________________________________________________________
Double_t AliAnalysisTaskCorrelationhhK0s::CalculateDeltaTheta( Double_t theta1, Double_t theta2 )  {  
  const double dtheta = theta2 - theta1;
  return dtheta;
}
//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskCorrelationhhK0s::CalculateDeltaPhi( Double_t phi1, Double_t phi2 )  {   //AliFemtoUser/AliFemtoPairCutDetaDphi.h
  const double dphi = phi2 - phi1;
  return dphi;
}

//_____________________________________________________________________________
void AliAnalysisTaskCorrelationhhK0s::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
Float_t AliAnalysisTaskCorrelationhhK0s::GetLengthInActiveZone( AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b ){
  // Input parameters:
  //   deltaY - user defined "dead region" in cm
  //   deltaZ - user defined "active region" in cm (250 cm drift lenght - 14 cm L1 delay
  //   b     - magnetic field 
  AliESDtrack esdTrack( gt );
  esdTrack.SetESDEvent((AliESDEvent*) gt->GetEvent() );
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(gt);
  esdTrack.ResetTrackParamIp(&etp);
  return esdTrack.GetLengthInActiveZone(1, deltaY, deltaZ, b);
}
//_____________________________________________________________________________
