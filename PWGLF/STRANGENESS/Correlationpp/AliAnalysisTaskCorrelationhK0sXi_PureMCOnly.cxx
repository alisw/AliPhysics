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
#include "AliAnalysisTaskCorrelationhK0sXi_PureMCOnly.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliAnalysisCorrelationEventCollection.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"

#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskCorrelationhK0sXi_PureMCOnly;   
using namespace std;          
ClassImp(AliAnalysisTaskCorrelationhK0sXi_PureMCOnly) 
const Float_t kctauval[2]={40, 30}; 
const Float_t Mass[2]={0.497611, 1.115683};

AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::AliAnalysisTaskCorrelationhK0sXi_PureMCOnly() :AliAnalysisTaskSE(), 
  fAnalysisType("AOD"), 
  fCollidingSystem("pp"), 
  fOutputList(0), 
  fSignalTree(0), 
  fBkgTree(0), 
  fOutputList2(0),
  fOutputList3(0),  
  fOutputList4(0),  
  fMCEvent(0), 
  fIshhCorr(0),
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
  fisInclusiveINELgt0(0),
  lPercentilesMin(0),
  lPercentilesMax(0),
  fisHM(0),
  fHistInelCorr(0),
  fHistEvtNoTrigger(0),
  fHistPt(0), 
  fHistPtV0(0), 
  fHistPtTMaxBefAllCfrDataMC(0),
  fHistPtTMinBefAll(0),
  fHistPtTMinBefAllMC(0),
  fHistPtvsMult(0), 
  fHistPtvsMultBefAll(0), 
  fHistPtMaxvsMult(0), 
  fHistPtMaxvsMultBefAll(0), 
  fHistZvertex(0),  
  fHistMultForwardvsMidRap(0),
  fHistMultForwardvsMidRapvsPt(0),
  fHistFractionSharedTPCClusters(0),
  fHistNumberChargedAllEvents(0),
  fHistNumberChargedNoTrigger(0),
  fHistNumberChargedTrigger(0),
  fHist_eta_phi(0),  
  fHist_eta_phi_PtMax(0),  
  fHist_multiplicityAllSelEvents(0),
  fHist_multiplicity(0),
  fHist_multiplicity_EvwTrigger(0),
  fHistEventMult(0), 
  fHistEventV0(0), 
  fHistEventV0Pt(0), 
  fHistTrack(0),
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
  fHistMultvsV0All(0),
  fHistMultvsV0AllTruth(0),
  fHistMultvsV0MCAll(0),
  fHistTriggerNotLeading(0),
  fHistTriggerNotLeadingMC(0),
  fHistMultvsTriggerBefAll(0),
  fHistMultvsTriggerMCTruthBefAll(0),
  fHistMultvsTrigger(0),
  fHistMultvsTriggerMCTruth(0),
  fHistTrigger(0),
  fHistTriggerMCTruth(0),
  fHistTriggerwV0(0),
  fHistTriggerwV0MCTruth(0),
  fHistMultiplicityVsVertexZ(0),
  fHistTriggervsMult(0),
  fHistTriggervsMultMC(0),
  fHistMultiplicityOfMixedEvent(0),
  fHistGeneratedV0Pt(0),
  fHistGeneratedTriggerPtPhi(0),
  fHistGeneratedV0PtTMaxPhi(0),
  fHistGeneratedTriggerPtEta(0),
  fHistGeneratedV0PtTMaxEta(0),
  fHistGeneratedV0PtPtTMax(0),
  fHistGeneratedV0PtEta(0),
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
  fTreeVariableMultiplicityV0A(0),                   
  fTreeVariableMultiplicityV0C(0),                   
  fTreeVariableZvertex(0),
  fTreeVariablePDGCodeTrigger(0),
  fTreeVariablePDGCodeAssoc(0),
  fTreeVariableLabelTrigger(0),
  fTreeVariableLabelPos(0),
  fTreeVariableLabelNeg(0),
  FifoShiftok(kFALSE)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purpos, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::AliAnalysisTaskCorrelationhK0sXi_PureMCOnly(const char* name) : AliAnalysisTaskSE(name),
										     fAnalysisType("AOD"), 
										     fCollidingSystem("pp"), 
										     fOutputList(0), 
										     fSignalTree(0), 
										     fBkgTree(0), 
										     fOutputList2(0), 
										     fOutputList3(0), 
										     fOutputList4(0), 
										     fMCEvent(0), 
										     fIshhCorr(0),
  fEventColl(0x0), 
  fEvt(0x0), 
  fzVertexBins(10), 
  fnMultBins(20),	 
  fMaxFirstMult(50),
  fMaxSecondMult(150),
  fnEventsToMix(50),
  fEtaTrigger(0.8),
  fEtahAssoc(0.8),								      fEtaV0Assoc(0.8),
  fFilterBitValue(128),
  fYear(2016),
  fisInclusiveINELgt0(0),
  lPercentilesMin(0),
  lPercentilesMax(0),
  fisHM(0),
  fHistInelCorr(0),
  fHistEvtNoTrigger(0),
  fHistPt(0), 
  fHistPtV0(0), 
  fHistPtTMaxBefAllCfrDataMC(0), 
  fHistPtTMinBefAll(0),
  fHistPtTMinBefAllMC(0),
  fHistPtvsMult(0), 
  fHistPtvsMultBefAll(0), 
  fHistPtMaxvsMult(0), 
  fHistPtMaxvsMultBefAll(0), 
  fHistZvertex(0),  
  fHistMultForwardvsMidRap(0),
  fHistMultForwardvsMidRapvsPt(0),
  fHistFractionSharedTPCClusters(0),
  fHistNumberChargedAllEvents(0),
  fHistNumberChargedNoTrigger(0),
  fHistNumberChargedTrigger(0),
  fHist_eta_phi(0),  
  fHist_eta_phi_PtMax(0),  
  fHist_multiplicityAllSelEvents(0),
  fHist_multiplicity(0),
  fHist_multiplicity_EvwTrigger(0),
  fHistEventMult(0), 
  fHistEventV0(0), 
  fHistEventV0Pt(0), 
  fHistTrack(0),
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
  fHistMultvsV0All(0),
  fHistMultvsV0AllTruth(0),
  fHistMultvsV0MCAll(0), 
  fHistTriggerNotLeading(0),
  fHistTriggerNotLeadingMC(0),
  fHistMultvsTriggerBefAll(0),
  fHistMultvsTriggerMCTruthBefAll(0),
  fHistMultvsTrigger(0),
  fHistMultvsTriggerMCTruth(0),
  fHistTrigger(0),
  fHistTriggerMCTruth(0),
  fHistTriggerwV0(0),
  fHistTriggerwV0MCTruth(0),
  fHistMultiplicityVsVertexZ(0),
  fHistTriggervsMult(0),
  fHistTriggervsMultMC(0),
  fHistMultiplicityOfMixedEvent(0),
  fHistGeneratedV0Pt(0),
  fHistGeneratedTriggerPtPhi(0),
  fHistGeneratedV0PtTMaxPhi(0),
  fHistGeneratedTriggerPtEta(0),
  fHistGeneratedV0PtTMaxEta(0),
  fHistGeneratedV0PtPtTMax(0),
  fHistGeneratedV0PtEta(0),
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
  fTreeVariableMultiplicityV0A(0),                   
  fTreeVariableMultiplicityV0C(0),                   
  fTreeVariableZvertex(0),
  fTreeVariablePDGCodeTrigger(0),
  fTreeVariablePDGCodeAssoc(0),
  fTreeVariableLabelTrigger(0),
  fTreeVariableLabelPos(0),
  fTreeVariableLabelNeg(0),
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
AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::~AliAnalysisTaskCorrelationhK0sXi_PureMCOnly()
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
  
//_____________________________________________________________________________
void AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::UserCreateOutputObjects()
{
  Int_t NumBinsMult=100;
  if (fisHM) NumBinsMult=100;
  Float_t UpperLimitMult =300;
  if (fisHM) UpperLimitMult = 300;

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
  
  const char* nameoutputSignalTree = GetOutputSlot(2)->GetContainer()->GetName();
  fSignalTree= new TTree(nameoutputSignalTree,"fSignalTree");
  fSignalTree->Branch("fTreeVariablePtTrigger",          &fTreeVariablePtTrigger   , "fTreeVariablePtTrigger/D");
  fSignalTree->Branch("fTreeVariableChargeTrigger",      &fTreeVariableChargeTrigger, "fTreeVariableChargeTrigger/I");
  fSignalTree->Branch("fTreeVariableEtaTrigger",         &fTreeVariableEtaTrigger  , "fTreeVariableEtaTrigger/D");
  fSignalTree->Branch("fTreeVariablePhiTrigger",         &fTreeVariablePhiTrigger, "fTreeVariablePhiTrigger/D");
  fSignalTree->Branch("fTreeVariableDCAz",               &fTreeVariableDCAz  , "fTreeVariableDCAz/D");
  fSignalTree->Branch("fTreeVariableDCAxy",              &fTreeVariableDCAxy  , "fTreeVariableDCAxy/D");
  fSignalTree->Branch("fTreeVariableLabelTrigger",       &fTreeVariableLabelTrigger  , "fTreeVariableLabelTrigger/I");
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
  fSignalTree->Branch("fTreeVariableMultiplicityV0A",    &fTreeVariableMultiplicityV0A , "fTreeVariableMultiplicityV0A/D");
  fSignalTree->Branch("fTreeVariableMultiplicityV0C",    &fTreeVariableMultiplicityV0C , "fTreeVariableMultiplicityV0C/D");
  fSignalTree->Branch("fTreeVariableZvertex",            &fTreeVariableZvertex  , "fTreeVariableZvertex/D");
  fSignalTree->Branch("fTreeVariablePDGCodeTrigger",     &fTreeVariablePDGCodeTrigger  , "fTreeVariablePDGCodeTrigger/I");
  fSignalTree->Branch("fTreeVariablePDGCodeAssoc",       &fTreeVariablePDGCodeAssoc  , "fTreeVariablePDGCodeAssoc/I");
  fSignalTree->Branch("fTreeVariableLabelPos",           &fTreeVariableLabelPos  , "fTreeVariableLabelPos/I");
  fSignalTree->Branch("fTreeVariableLabelNeg",           &fTreeVariableLabelNeg  , "fTreeVariableLabelNeg/I");

  const char* nameoutputBkgTree = GetOutputSlot(3)->GetContainer()->GetName();
  fBkgTree= new TTree(nameoutputBkgTree,"fBkgTree");
  fBkgTree->Branch("fTreeVariablePtTrigger",          &fTreeVariablePtTrigger   , "fTreeVariablePtTrigger/D");
  fBkgTree->Branch("fTreeVariableChargeTrigger",      &fTreeVariableChargeTrigger, "fTreeVariableChargeTrigger/I");
  fBkgTree->Branch("fTreeVariableEtaTrigger",         &fTreeVariableEtaTrigger  , "fTreeVariableEtaTrigger/D");
  fBkgTree->Branch("fTreeVariablePhiTrigger",         &fTreeVariablePhiTrigger, "fTreeVariablePhiTrigger/D");
  fBkgTree->Branch("fTreeVariableDCAz",               &fTreeVariableDCAz  , "fTreeVariableDCAz/D");
  fBkgTree->Branch("fTreeVariableDCAxy",              &fTreeVariableDCAxy  , "fTreeVariableDCAxy/D");
  fBkgTree->Branch("fTreeVariableChargeAssoc",        &fTreeVariableChargeAssoc  , "fTreeVariableChargeAssoc/I");
  fBkgTree->Branch("fTreeVariableLabelTrigger",       &fTreeVariableLabelTrigger  , "fTreeVariableLabelTrigger/I");
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
  fBkgTree->Branch("fTreeVariableMultiplicityV0A",    &fTreeVariableMultiplicityV0A , "fTreeVariableMultiplicityV0A/D");
  fBkgTree->Branch("fTreeVariableMultiplicityV0C",    &fTreeVariableMultiplicityV0C , "fTreeVariableMultiplicityV0C/D");
  fBkgTree->Branch("fTreeVariableZvertex",            &fTreeVariableZvertex  , "fTreeVariableZvertex/D");
  fBkgTree->Branch("fTreeVariablePDGCodeTrigger",     &fTreeVariablePDGCodeTrigger  , "fTreeVariablePDGCodeTrigger/I");
  fBkgTree->Branch("fTreeVariablePDGCodeAssoc",       &fTreeVariablePDGCodeAssoc  , "fTreeVariablePDGCodeAssoc/I");
  fBkgTree->Branch("fTreeVariableLabelPos",           &fTreeVariableLabelPos  , "fTreeVariableLabelPos/I");
  fBkgTree->Branch("fTreeVariableLabelNeg",           &fTreeVariableLabelNeg  , "fTreeVariableLabelNeg/I");
 
  fHistPt = new TH1F("fHistPt", "p_{T} distribution of selected charged tracks in events used for AC", 300, 0, 30); 
  fHistPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  
  fHistPtV0 = new TH1F("fHistPtV0", "p_{T} distribution of selected V0 in events used for AC", 300, 0, 30); 
  fHistPtV0->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMaxBefAllCfrDataMC= new TH2F("fHistPtTMaxBefAllCfrDataMC", "p_{T} distribution of trigger particle with Pt maximum in the event", 300, 0, 30, 300, 0,30); 
  fHistPtTMaxBefAllCfrDataMC->GetXaxis()->SetTitle("p^{Trigger, Max}_{T} (reco)(GeV/c)");
  fHistPtTMaxBefAllCfrDataMC->GetYaxis()->SetTitle("p^{Trigger, Max}_{T} (MC)(GeV/c)");

  fHistPtTMinBefAll = new TH1F("fHistPtTMinBefAll", "p_{T} distribution of reco trigger particle with Pt minimum in the event", 300, 0, 30); 
  fHistPtTMinBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMinBefAllMC = new TH1F("fHistPtTMinBefAllMC", "p_{T} distribution of true trigger particle with Pt minimum in the event", 300, 0, 30); 
  fHistPtTMinBefAllMC->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtvsMultBefAll= new TH2F("fHistPtvsMultBefAll", "p_{T} and centrality distribution of charged tracks in events w T>0", 300, 0, 30, NumBinsMult, 0, UpperLimitMult); 
  fHistPtvsMultBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtvsMultBefAll->GetYaxis()->SetTitle("Centrality");

  fHistPtvsMult= new TH2F("fHistPtvsMult", "p_{T} and centrality distribution of charged tracks in events used for AC", 300, 0, 30, NumBinsMult, 0, UpperLimitMult); 
  fHistPtvsMult->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtvsMult->GetYaxis()->SetTitle("Centrality");

  fHistPtMaxvsMultBefAll= new TH2F("fHistPtMaxvsMultBefAll", "p_{T} and centrality distribution of charged tracks with maxiumum pt in events w T>0", 600, 0, 30, NumBinsMult, 0, UpperLimitMult); 
  fHistPtMaxvsMultBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtMaxvsMultBefAll->GetYaxis()->SetTitle("Centrality");

  fHistPtMaxvsMult= new TH2F("fHistPtMaxvsMult", "p_{T} and centrality distribution of charged tracks with maximum pT in events used for AC)", 600, 0, 30, NumBinsMult, 0, UpperLimitMult); 
  fHistPtMaxvsMult->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtMaxvsMult->GetYaxis()->SetTitle("Centrality");

  fHistZvertex= new TH1F("fHistZvertex", "Z vertex distribution of selected events used for AC", 40,-20,20);

  fHistMultForwardvsMidRap = new TH2F("fHistMultForwardvsMidRap", "Charged particles in V0 acceptance vs dN/deta at midrapidity in events with a trigger particle", NumBinsMult,0,UpperLimitMult,300, 0, 300);
  fHistMultForwardvsMidRap->GetXaxis()->SetTitle("N charged particles at midrapidity (|#eta|<0.5)");
  fHistMultForwardvsMidRap->GetYaxis()->SetTitle("N charged particles in V0 acceptance");

  fHistMultForwardvsMidRapvsPt = new TH3F("fHistMultForwardvsMidRapvsPt", "Charged particles in V0 acceptance vs dN/deta at midrapidity in events with a trigger particle", NumBinsMult,0,UpperLimitMult,300, 0, 300, 300, 0, 30);
  fHistMultForwardvsMidRapvsPt->GetXaxis()->SetTitle("N charged particles at midrapidity (|#eta|<0.5)");
  fHistMultForwardvsMidRapvsPt->GetYaxis()->SetTitle("N charged particles in V0 acceptance");
  fHistMultForwardvsMidRapvsPt->GetZaxis()->SetTitle("p_{T} trigger particle");

  fHistFractionSharedTPCClusters = new TH1F ("fHistFractionSharedTPCClusters", "fHistFractionSharedTPCClusters",100, 0,1);

  fHistNumberChargedAllEvents=new TH3F("fHistNumberChargedAllEvents", "fHistNumberChargedAllEvents", NumBinsMult,0,UpperLimitMult, 100,0,100, 60, 0,30);
  fHistNumberChargedAllEvents->GetXaxis()->SetTitle("Multiplicity class");
  fHistNumberChargedAllEvents->GetYaxis()->SetTitle("Number of charged primary particles");
  fHistNumberChargedAllEvents->GetZaxis()->SetTitle("p^{Trigg, Max}_{T} (GeV/c)");

  fHistNumberChargedNoTrigger=new TH3F("fHistNumberChargedNoTrigger", "fHistNumberChargedNoTrigger", NumBinsMult,0,UpperLimitMult, 100,0,100, 60, 0,30);
  fHistNumberChargedNoTrigger->GetXaxis()->SetTitle("Multiplicity class");
  fHistNumberChargedNoTrigger->GetYaxis()->SetTitle("Number of charged primary particles");
  fHistNumberChargedNoTrigger->GetZaxis()->SetTitle("p^{Trigg, Max}_{T} (GeV/c)");

  fHistNumberChargedTrigger=new TH3F("fHistNumberChargedTrigger", "fHistNumberChargedTrigger", NumBinsMult,0,UpperLimitMult, 100,0,100, 60, 0,30);
  fHistNumberChargedTrigger->GetXaxis()->SetTitle("Multiplicity class");
  fHistNumberChargedTrigger->GetYaxis()->SetTitle("Number of charged primary particles");
  fHistNumberChargedTrigger->GetZaxis()->SetTitle("p^{Trigg, Max}_{T} (GeV/c)");

  fHist_eta_phi= new TH2F("fHist_eta_phi", "Distribution of charged tracks in events used for AC", 400, 0,2*TMath::Pi(), 400,-0.8, 0.8);
  fHist_eta_phi->GetYaxis()->SetTitle("Eta");
  fHist_eta_phi->GetXaxis()->SetTitle("Phi (radians)"); 

  fHist_eta_phi_PtMax= new TH2F("fHist_eta_phi_PtMax", "Distribution of charged tracks with maximum Pt in events used for AC", 400, 0,2*TMath::Pi(), 400,-0.8, 0.8);
  fHist_eta_phi_PtMax->GetYaxis()->SetTitle("Eta");
  fHist_eta_phi_PtMax->GetXaxis()->SetTitle("Phi (radians)"); 

  fHist_multiplicityAllSelEvents=new TH1F("fHist_multiplicityAllSelEvents", "fHist_multiplicityAllSelEvents", NumBinsMult, 0, UpperLimitMult); 
  fHist_multiplicityAllSelEvents->SetTitle("Centrality distribution of selected INT7/HM events");
  fHist_multiplicity=new TH1F("fHist_multiplicity", "fHist_multiplicity", NumBinsMult, 0, UpperLimitMult); 
  fHist_multiplicity->SetTitle("Centrality distribution of events used for AC");
  fHist_multiplicity_EvwTrigger= new TH1F("fHist_multiplicity_EvwTrigger", "fHist_multiplicity_EvwTrigger", NumBinsMult, 0, UpperLimitMult); 
  fHist_multiplicity_EvwTrigger->SetTitle("Centrality distribution of events with NT>0");

  fHistPDG=new TH1F("fHistPDG", "fHistPDG",6400, -3200, 3200);

  fHistTrackBufferOverflow = new TH1F("fHistTrackBufferOverflow","",2,0,2);
  
  fHistEventMult=new TH1F("fHistEventMult", "fHistEventMult", 24, 0.5, 24.5);
  fHistEventMult->GetXaxis()->SetBinLabel(1,"AllEvents");
  fHistEventMult->GetXaxis()->SetBinLabel(2,"Events w/|Vx|<10 cm"); 
  fHistEventMult->GetXaxis()->SetBinLabel(3,"MCStack"); 
  fHistEventMult->GetXaxis()->SetBinLabel(4,"INEL>0"); 
  fHistEventMult->GetXaxis()->SetBinLabel(5,"isSelected"); 
  fHistEventMult->GetXaxis()->SetBinLabel(6,"Ntrigger>0"); 
  fHistEventMult->GetXaxis()->SetBinLabel(7,"NTrigger>1");
  fHistEventMult->GetXaxis()->SetBinLabel(8,"SelEv (ACEvents)"); //all the events used for angular correlation

  fHistEventV0=new TH1F("fHistEventV0", "fHistEventV0",31, 0.5, 31.5);
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
  fHistEventV0->GetXaxis()->SetBinLabel(29,"All hybrid"); 
  fHistEventV0->GetXaxis()->SetBinLabel(30,"Not mother of trigger"); 
  fHistEventV0->GetXaxis()->SetBinLabel(31,"Trigger not from K0s"); 

  fHistEventV0Pt=new TH2F("fHistEventV0Pt", "fHistEventV0Pt",18, 0.5, 18.5, 300, 0 ,30);
  fHistEventV0Pt->GetXaxis()->SetBinLabel(1,"All hybrid"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(2,"Mother of trigger"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(3,"Trigger from K0s"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(4,"2 && !3"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(5,"!2 && 3"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(6,"2 && 3"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(7,"All gen (eff)"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(8,"Mother of trigger"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(9,"Trigger from K0s"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(10,"2 && !3"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(11,"!2 && 3"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(12,"2 && 3"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(13,"All reco (MC)"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(14,"Mother of trigger"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(15,"Trigger from K0s"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(16,"2 && !3"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(17,"!2 && 3"); 
  fHistEventV0Pt->GetXaxis()->SetBinLabel(18,"2 && 3"); 

  fHistTrack=new TH1F("fHistTrack", "fHistTrack", 19, 0.5, 19.5);
  fHistTrack->GetXaxis()->SetBinLabel(1,"All tracks");
  fHistTrack->GetXaxis()->SetBinLabel(2,"Charged tracks"); 
  fHistTrack->GetXaxis()->SetBinLabel(3,"Tracks with |eta| < 0.8"); 
  fHistTrack->GetXaxis()->SetBinLabel(4,"Primary"); 
  fHistTrack->GetXaxis()->SetBinLabel(5,"pKpiemu"); 
  fHistTrack->GetXaxis()->SetBinLabel(6,"pt selected"); 
  fHistTrack->GetXaxis()->SetBinLabel(7,"N.trigger MC");
  fHistTrack->GetXaxis()->SetBinLabel(8,"N.trigger MC (NT>0)");
  fHistTrack->GetXaxis()->SetBinLabel(9,"N.trigger MC (NV0>0)");

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
  
  fHistMultvsTrigger=new TH2F("fHistMultvsTrigger", "Centrality of selected events (T>0, V>0) vs number of trigger particles", 30, -0.5, 29.5, NumBinsMult, 0, UpperLimitMult);
  fHistMultvsTrigger->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTrigger->GetYaxis()->SetTitle("Centrality");
  

  fHistMultvsTriggerMCTruth=new TH2F("fHistMultvsTriggerMCTruth", "Centrality of selected events (T>0, V>0) vs number of trigger particles, MC Truth", 30, -0.5, 29.5, NumBinsMult, 0, UpperLimitMult);
  fHistMultvsTriggerMCTruth->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerMCTruth->GetYaxis()->SetTitle("Centrality");
  

  fHistMultvsTriggerBefAll=new TH2F("fHistMultvsTriggerBefAll", "Centrality of events vs number of trigger particles", 30, -0.5, 29.5, NumBinsMult, 0, UpperLimitMult);
  fHistMultvsTriggerBefAll->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerBefAll->GetYaxis()->SetTitle("Centrality");
  
 
  fHistMultvsTriggerMCTruthBefAll=new TH2F("fHistMultvsTriggerMCTruthBefAll", "Centrality of events vs number of trigger particles, MC Truth", 30, -0.5, 29.5, NumBinsMult, 0, UpperLimitMult);
  fHistMultvsTriggerMCTruthBefAll->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerMCTruthBefAll->GetYaxis()->SetTitle("Centrality");
  
  fHistMultvsV0All=new TH2F("fHistMultvsV0All", "Centrality of events w T>0 vs number of reco V0s",60, -0.5, 59.5,NumBinsMult, 0, UpperLimitMult );
  fHistMultvsV0All->GetXaxis()->SetTitle("Number of V0 particles");
  fHistMultvsV0All->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0AllTruth=new TH2F("fHistMultvsV0AllTruth", "Centrality of events w T>0 vs number of reco true V0s",60, -0.5, 59.5,NumBinsMult, 0, UpperLimitMult );
  fHistMultvsV0AllTruth->GetXaxis()->SetTitle("Number of V0 reco particles");
  fHistMultvsV0AllTruth->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0MCAll=new TH2F("fHistMultvsV0MCAll", "Centrality of events w T>0 vs number of true V0s",60, -0.5, 59.5,NumBinsMult, 0, UpperLimitMult );
  fHistMultvsV0MCAll->GetXaxis()->SetTitle("Number of V0 true particles");
  fHistMultvsV0MCAll->GetYaxis()->SetTitle("Centrality");

  fHistTriggerNotLeading=new TH3F("fHistTriggerNotLeading", "Events with trigger not leading in all events with NT>0",60, -0.5, 59.5,NumBinsMult, 0, UpperLimitMult, 60, 0, 30 );
  fHistTriggerNotLeading->GetXaxis()->SetTitle("Number of V0 with p_{T}> p_{T} Trigger");
  fHistTriggerNotLeading->GetYaxis()->SetTitle("Multiplicity class");
  fHistTriggerNotLeading->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");

  fHistTriggerNotLeadingMC=new TH3F("fHistTriggerNotLeadingMC", "Events with trigger not leading in all events with NT>0 (MC Truth)",60, -0.5, 59.5,NumBinsMult, 0, UpperLimitMult, 60, 0, 30 );
  fHistTriggerNotLeadingMC->GetXaxis()->SetTitle("Number of V0 with p_{T}> p_{T} Trigger");
  fHistTriggerNotLeadingMC->GetYaxis()->SetTitle("Multiplicity class");
  fHistTriggerNotLeadingMC->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");

  fHistInelCorr=new TH2F("fHistInelCorr", "Correlation between #tracklets and #tracks", 200, 0, 1000, 200, 0, 1000);
  fHistInelCorr->GetXaxis()->SetTitle("Tracks");
  fHistInelCorr->GetYaxis()->SetTitle("Tracklets");

  fHistEvtNoTrigger=new TH2F("fHistEvtNoTrigger", "#Tracks and multiplicity of events w/o trigger", 200, 0, 1000, NumBinsMult, 0, UpperLimitMult);
  fHistEvtNoTrigger->GetXaxis()->SetTitle("#Tracks");
  fHistEvtNoTrigger->GetYaxis()->SetTitle("Percentile");
  
  fHistTrigger=new TH1F("fHistTrigger", "Number of reco trigger particle distribution for selected events (also T=0)", 60, -0.5, 59.5); // each entry is an event

  fHistTriggerwV0=new TH1F("fHistTriggerwV0", "Number of reco trigger particle distribution for events used for AC", 60, -0.5, 59.5); // each entry is an event

  fHistTriggerMCTruth=new TH1F("fHistTriggerMCTruth", "Number of true trigger particle distribution for selected events (also T=0)", 60, -0.5, 59.5); // each entry is an event

  fHistTriggerwV0MCTruth=new TH1F("fHistTriggerwV0MCTruth", "Number of true trigger particle distribution for events used for AC", 60, -0.5, 59.5); // each entry is an event

  fHistMultiplicityVsVertexZ=new TH2F("fHistMultiplicityVsVertexZ", "Centrality vs Z vertex of selected events with NT>0 and NV0>0 ",  20, -10, 10,NumBinsMult, 0, UpperLimitMult);
      
  fHistTriggervsMult=new TH1F("fHistTriggervsMult", "Numero di particelle di trigger nei vari intervalli di centralita'", NumBinsMult, 0, UpperLimitMult);
  fHistTriggervsMult->GetXaxis()->SetTitle("Centrality");

  fHistTriggervsMultMC=new TH1F("fHistTriggervsMultMC", "Numero di particelle di trigger (MCtruth) nei vari intervalli di centralita'", NumBinsMult, 0, UpperLimitMult);
  fHistTriggervsMultMC->GetXaxis()->SetTitle("Centrality");

  fHistGeneratedTriggerPtPhi=new TH3F("fHistGeneratedTriggerPtPhi", "p_{T} and #phi distribution of generated trigger particles (charged, primary)", 600, 0, 30, 400,0, 2*TMath::Pi(),  NumBinsMult, 0, UpperLimitMult );
  fHistGeneratedTriggerPtPhi->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedTriggerPtPhi->GetYaxis()->SetTitle("#phi");

  fHistGeneratedTriggerPtEta=new TH3F("fHistGeneratedTriggerPtEta", "p_{T} and #eta distribution of generated trigger particles (primary, charged)", 600, 0, 30, 400,-1.2,1.2,  NumBinsMult, 0, UpperLimitMult );
  fHistGeneratedTriggerPtEta->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedTriggerPtEta->GetYaxis()->SetTitle("#eta");

  fHistGeneratedV0Pt= new TH2F("fHistGeneratedV0Pt", "p_{T} distribution of generated trigger particles from pileup (charged, primary)", 600, 0, 30, NumBinsMult, 0, UpperLimitMult );
  
  fHistGeneratedV0PtTMaxPhi=new TH3F("fHistGeneratedV0PtTMaxPhi", "p^{Trigg, Max}_{T} and #phi distribution of generated V0 particles (K0s, primary, events w T>0)", 120, -30, 30, 400,0, 2*TMath::Pi(),  NumBinsMult, 0, UpperLimitMult);
  fHistGeneratedV0PtTMaxPhi->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
  fHistGeneratedV0PtTMaxPhi->GetYaxis()->SetTitle("#phi");
  
  fHistGeneratedV0PtTMaxEta=new TH3F("fHistGeneratedV0PtTMaxEta", "p^{Trigg, Max}_{T} and #eta distribution of generated V0 particles (K0s, primary, events w T>0)", 120,-30, 30, 240,-1.2,1.2,  NumBinsMult, 0, UpperLimitMult );
  fHistGeneratedV0PtTMaxEta->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
  fHistGeneratedV0PtTMaxEta->GetYaxis()->SetTitle("#eta");
 
  fHistGeneratedV0PtPtTMax=new TH3F("fHistGeneratedV0PtPtTMax", "p_{T} and p^{Trigg, Max}_{T} distribution of generated V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 120, -30, 30,  NumBinsMult, 0, UpperLimitMult );
  fHistGeneratedV0PtPtTMax->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedV0PtPtTMax->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");
  
  fHistGeneratedV0PtEta=new TH3F("fHistGeneratedV0PtEta", "p_{T} and #eta distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 240,-1.2,1.2,  NumBinsMult, 0, UpperLimitMult );
  fHistGeneratedV0PtEta->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedV0PtEta->GetYaxis()->SetTitle("#eta");

  fHistMultiplicityOfMixedEvent=new TH2F("fHistMultiplicityOfMixedEvent", "Distribution of number of events used for the mixing", 100, 0.5, 100.5, NumBinsMult, 0, UpperLimitMult);

  fOutputList->Add(fHistPDG);
  fOutputList->Add(fHistTrackBufferOverflow);
  fOutputList->Add(fHistEventMult);
  fOutputList->Add(fHistEventV0);
  fOutputList->Add(fHistEventV0Pt);
  fOutputList->Add(fHistTrack); 
  fOutputList->Add(fHistIsCommonParton);
  fOutputList->Add(fHistCommonParton);
  fOutputList->Add(fHistCommonPartonTrueCasc);
  fOutputList->Add(fHistTriggerComposition); 
  fOutputList->Add(fHistTriggerCompositionMCTruth); 
  fOutputList->Add(fHistAssocComposition); 
  fOutputList->Add(fHistAssocCompositionMCTruth); 
  
  //istogrammi riempiti ad ogni evento selezionato (ossia utilizzato per AC) con una entry
  fOutputList->Add(fHistZvertex);
  fOutputList->Add(fHistMultForwardvsMidRap);
  fOutputList->Add(fHistMultForwardvsMidRapvsPt);
  fOutputList->Add(fHistNumberChargedAllEvents);
  fOutputList->Add(fHistNumberChargedTrigger);
  fOutputList->Add(fHistNumberChargedNoTrigger);
  fOutputList->Add(fHist_multiplicityAllSelEvents); 
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
  fOutputList->Add(fHistInelCorr);   
  fOutputList->Add(fHistEvtNoTrigger);   
  fOutputList->Add(fHistPt);   
  fOutputList->Add(fHistPtV0);     
  fOutputList->Add(fHistPtTMaxBefAllCfrDataMC);
  fOutputList->Add(fHistPtTMinBefAll);     
  fOutputList->Add(fHistPtTMinBefAllMC);   
  fOutputList->Add(fHistPtvsMult);       
  fOutputList->Add(fHistPtvsMultBefAll);       
  fOutputList->Add(fHistPtMaxvsMult);       
  fOutputList->Add(fHistPtMaxvsMultBefAll);       
  fOutputList->Add(fHist_eta_phi);
  fOutputList->Add(fHist_eta_phi_PtMax);
  fOutputList->Add(fHistTriggervsMult);
  fOutputList->Add(fHistTriggervsMultMC);
  fOutputList->Add(fHistTrigger);
  fOutputList->Add(fHistTriggerwV0);
  fOutputList->Add(fHistTriggerMCTruth);
  fOutputList->Add(fHistTriggerwV0MCTruth);

  fOutputList->Add(fHistMultiplicityOfMixedEvent);
  fOutputList3->Add(fHistGeneratedTriggerPtPhi);
  fOutputList3->Add(fHistGeneratedTriggerPtEta);

  fOutputList2->Add(fHistGeneratedV0Pt);
  fOutputList2->Add(fHistGeneratedV0PtTMaxPhi); 
  fOutputList2->Add(fHistGeneratedV0PtTMaxEta); 
  fOutputList2->Add(fHistGeneratedV0PtPtTMax); 
  fOutputList2->Add(fHistGeneratedV0PtEta);

  PostData(1, fOutputList);  
  PostData(2, fSignalTree);       
  PostData(3, fBkgTree); 
  PostData(4, fOutputList2);     
  PostData(5, fOutputList3);     
  PostData(6, fOutputList4);     
     
}
//_____________________________________________________________________________
void AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::UserExec(Option_t *)
{

  Float_t LastzBin;
  Float_t LastcentralityBin;
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};

  fHistEventMult->Fill(1);   

  //PVz
  fMCEvent= MCEvent();
  if (!fMCEvent) {
    Printf("ERROR: Could not retrieve MC event \n");
    return;
  }

  TArrayF mcPrimaryVtx;
  AliGenEventHeader* mcHeader=fMCEvent->GenEventHeader();
  if(!mcHeader){
    Printf("ERROR: Could not retrieve MC header \n");
    return;
  }
  mcHeader->PrimaryVertex(mcPrimaryVtx);
  lBestPrimaryVtxPos[2]= mcPrimaryVtx.At(2);
  //  cout << "Z vertex position " << lBestPrimaryVtxPos[2] << endl;

  if (TMath::Abs(lBestPrimaryVtxPos[2]) > 10.){
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3, fBkgTree);
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);
    return;
  }
  fHistEventMult->Fill(2);   

  //stack

  AliStack    *lMCstack  = 0x0;
  lMCstack = fMCEvent->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3, fBkgTree);
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);
    return;
  }

  fHistEventMult->Fill(3);   

  // Multiplicity Information Acquistion    
  Float_t lNchEta5   = 0;
  Float_t lNchEta8   = 0;
  Float_t lNchEta8to15   = 0;
  Float_t lNchEta10  = 0;
  Float_t lNchVZEROA = 0;
  Float_t lNchVZEROC = 0;
  Bool_t lEvSel_INELgtZEROStackPrimaries=kFALSE;
    
  //----- Loop on Stack ----------------------------------------------------------------
  //  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++) //stack 
  //  cout << "Loop on particles: " << fMCEvent->GetNumberOfTracks() << endl;
  //  cout << "Particles in stack: " << lMCstack->GetNtrack() << endl;
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < fMCEvent->GetNumberOfTracks(); iCurrentLabelStack++)
    {   // This is the begining of the loop on tracks
      //      TParticle* particleOne = lMCstack->Particle(iCurrentLabelStack); //stack
      TParticle* particleOne = static_cast<TParticle*>(fMCEvent->Particle(iCurrentLabelStack));
      if(!particleOne) continue;
      if(!particleOne->GetPDG()) continue;
      Double_t lThisCharge = particleOne->GetPDG()->Charge()/3.;
      if(TMath::Abs(lThisCharge)<0.001) continue;
      if(! (fMCEvent->IsPhysicalPrimary(iCurrentLabelStack)) ) continue;

      //Double_t gpt = particleOne -> Pt();
      Double_t geta = particleOne -> Eta();
      if( TMath::Abs(geta) < 0.5 ) lNchEta5++;
      if( TMath::Abs(geta) < 0.8 ) lNchEta8++;
      if( (TMath::Abs(geta) > 0.8) && (TMath::Abs(geta) < 1.5) ) lNchEta8to15++;
      if( TMath::Abs(geta) < 1.0 ) lNchEta10++;
      if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZEROStackPrimaries = kTRUE;
      if( 2.8 < geta && geta < 5.1 ) lNchVZEROA++;
      if(-3.7 < geta && geta <-1.7 ) lNchVZEROC++;
    }//End of loop on tracks

  if (lEvSel_INELgtZEROStackPrimaries)   fHistEventMult->Fill(4);   
  if (fisInclusiveINELgt0){
    if (!lEvSel_INELgtZEROStackPrimaries) {
      PostData(1, fOutputList);
      PostData(2, fSignalTree );
      PostData(3, fBkgTree);
      PostData(4, fOutputList2);  
      PostData(5, fOutputList3);     
      PostData(6, fOutputList4);
      return;
    }
  }

  Float_t  lPercentiles = lNchVZEROA + lNchVZEROC;

  if (fisHM) {
    if (lPercentiles < 120) {
      PostData(1, fOutputList);
      PostData(2, fSignalTree );
      PostData(3, fBkgTree);
      PostData(4, fOutputList2);  
      PostData(5, fOutputList3);     
      PostData(6, fOutputList4);
      return;
    }
  }
  fHistEventMult->Fill(5);   

  fHist_multiplicityAllSelEvents->Fill(lPercentiles);

  Bool_t Generated=kTRUE; //TRUE if generated particles are analyzed
  Int_t labelPrimOrSec=0; 
  Int_t label = 0; //track label
  Float_t TrackLength=-999;
  Bool_t isV0=kFALSE;
  int fieldsign =0;

  //Collect events according to their centrality and z coordinate of PV
  int zBin=0;
  int centralityBin=0;    
  double zStep=2*10/double(fzVertexBins), zStart=-10.;
    
  for (int i=0; i<fzVertexBins; i++) {
    if ((lBestPrimaryVtxPos[2] > zStart+i*zStep) && (lBestPrimaryVtxPos[2] < zStart+(i+1)*zStep)) {
      zBin=i;
      break;
    }
  } 

  if (fisHM){
    if(lPercentiles > 300) centralityBin=19; 
    else if(lPercentiles > 240) centralityBin=18;
    else if(lPercentiles > 220) centralityBin=17;
    else if(lPercentiles > 200) centralityBin=16;
    else if(lPercentiles > 180) centralityBin=15;
    else if(lPercentiles > 160) centralityBin=14;
    else if(lPercentiles > 140) centralityBin=13;
    else if(lPercentiles > 120) centralityBin=12;
    else centralityBin = 11;
  }
  else {
    if(lPercentiles > 250) centralityBin=19; 
    else if(lPercentiles > 225) centralityBin=18;
    else if(lPercentiles > 200) centralityBin=17;
    else if(lPercentiles > 175) centralityBin=16;
    else if(lPercentiles > 150) centralityBin=15;
    else if(lPercentiles > 120) centralityBin=14;
    else if(lPercentiles > 105) centralityBin=13;
    else if(lPercentiles > 90) centralityBin=12;
    else if(lPercentiles > 81) centralityBin=11;
    else if(lPercentiles > 72) centralityBin=10;
    else if(lPercentiles > 63) centralityBin=9;
    else if(lPercentiles > 54) centralityBin=8;
    else if(lPercentiles > 45) centralityBin=7;
    else if(lPercentiles > 30) centralityBin=6;
    else centralityBin = 5;
  }
  if (((centralityBin+1) >fnMultBins) || ((zBin+1) > fzVertexBins)){ 
    //c cout<<" ##################  WARNING: I'm going to break bacause of dimensional issues ########################"<<endl;
  }


  fEventColl[zBin][centralityBin]->FifoClear();
  fEvt = fEventColl[zBin][centralityBin]->fEvt;


  //***********************************************************************************
  //-----------------------------------LOOP OVER THE TRACKS
  //***********************************************************************************
  AliVTrack *vtrackg = 0x0;
  AliVTrack *vtrack = 0x0;
  AliVTrack *vtrackgPos = 0x0;
  AliVTrack *vtrackgNeg = 0x0;

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
  Int_t VPdgTrig[50]={0};
  Int_t VParticleTrigLabel[50]={0};
  TClonesArray* AODMCTrackArray =0x0;  
  Float_t ptTriggerMinimoMC=10000;
  Double_t ptTriggerMassimoMC=0;
  Double_t etaTriggerMassimoMC=0;
  Double_t phiTriggerMassimoMC=0;

  TriggerPdgCode=0;
  
  //begin loop for trigger particles (MC truth analysis)
  if (fMCEvent){
    //    cout <<"Looping over the tracks: " << fMCEvent->GetNumberOfTracks()<< endl;
    for ( Int_t i = 0; i< fMCEvent->GetNumberOfTracks(); i++ )    {
      TParticle* trParticle = static_cast<TParticle*>( fMCEvent->Particle(i) );
      if (!trParticle) continue;       
      if (!trParticle->GetPDG()) continue;       
      fHistTrack->Fill(1);
      if((trParticle->GetPDG()->Charge())==0)continue;
      fHistTrack->Fill(2);
      if(TMath::Abs(trParticle->Eta())>fEtaTrigger)continue; //I need to select particles within this eta range!
      fHistTrack->Fill(3);
      if (!fMCEvent->IsPhysicalPrimary(i)) continue;
      fHistTrack->Fill(4);
      if (!(TMath::Abs(trParticle ->GetPdgCode()) == 211 || TMath::Abs(trParticle ->GetPdgCode()) == 321|| TMath::Abs(trParticle ->GetPdgCode()) ==2212 || TMath::Abs(trParticle ->GetPdgCode()) == 11|| TMath::Abs(trParticle ->GetPdgCode()) == 13)) continue;
      fHistTrack->Fill(5);
      NumberFirstParticleAllPtMC++; 
      if(trParticle->Pt()<= fminPtj || trParticle->Pt()>=fmaxPtj) continue;
      fHistTrack->Fill(6);
      fHistGeneratedTriggerPtPhi->Fill(trParticle->Pt(), trParticle->Phi(), lPercentiles);
      fHistGeneratedTriggerPtEta->Fill(trParticle->Pt(), trParticle->Eta(), lPercentiles);

      if(trParticle->Pt()< ptTriggerMinimoMC) ptTriggerMinimoMC=trParticle->Pt();
      if(trParticle->Pt()> ptTriggerMassimoMC){
	ptTriggerMassimoMC=trParticle->Pt();
	etaTriggerMassimoMC=trParticle->Eta();
	phiTriggerMassimoMC=trParticle->Phi();
	TriggerPdgCode =trParticle ->GetPdgCode();
      }
      NumberFirstParticleMC++;
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fCharge       = trParticle->GetPDG()->Charge()/3.;
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fPt           = trParticle->Pt();
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fEta          = trParticle->Eta();
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fPhi          = trParticle->Phi();
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fTheta        = trParticle->Theta();
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fDCAz         = 0;
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fDCAxy        = 0;
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fMultiplicity = lPercentiles;
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fZvertex      = lBestPrimaryVtxPos[2];
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].isP           = 1;
      fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fPDGcode      = trParticle->GetPdgCode();
    }
    //      cout << "I found " << 	NumberFirstParticleMC << " trigger particles (pT > 3 GeV/c, not only the highest pt one)" << endl;
    //      cout << " pt of trigger particle " << 	  ptTriggerMassimoMC << endl;
  }

  if (ptTriggerMassimoMC!=0)  fHistPtMaxvsMultBefAll->Fill(ptTriggerMassimoMC, lPercentiles);  
  if (fisInclusiveINELgt0) fHistPtMaxvsMultBefAll->Fill(0.0001, lPercentiles); //fake pt
  fHistPtTMinBefAllMC->Fill(ptTriggerMinimoMC);

  fHistTriggerMCTruth->Fill(NumberFirstParticleMC);
  fHistMultvsTriggerMCTruthBefAll->Fill(NumberFirstParticleMC, lPercentiles);

  fHistTrack->AddBinContent(7, NumberFirstParticleMC);
  if (NumberFirstParticleMC>0) fHistEventMult->Fill(6);   
  if (NumberFirstParticleMC>1)fHistEventMult->Fill(7); 
  
  NumberFirstParticleAll=NumberFirstParticleMC; ptTriggerMassimoAll=ptTriggerMassimoMC; NumberFirstParticleAllPt = NumberFirstParticleAllPtMC;
  if(NumberFirstParticleAll==0) fHistNumberChargedNoTrigger->Fill(lPercentiles, NumberCharged,0 );
  if(NumberFirstParticleAll!=0) fHistNumberChargedTrigger->Fill(lPercentiles, NumberCharged,ptTriggerMassimoAll );
  fHistNumberChargedAllEvents->Fill(lPercentiles, NumberCharged,ptTriggerMassimoAll );

  if(NumberFirstParticleAllPt==0){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3, fBkgTree);
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);
    // cout  << "event does not have Trigger particles " << endl;          
    return;
  }

  if(NumberFirstParticleAll==0){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3, fBkgTree);
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    PostData(6, fOutputList4);
    //    cout  << "event does not have Trigger particles with pT> 3 GeV/c " << endl;          
    return;
  }

  fHistMultForwardvsMidRap->Fill(lNchEta5, lPercentiles);
  fHistMultForwardvsMidRapvsPt->Fill(lNchEta5, lPercentiles, ptTriggerMassimoMC);
  fHist_multiplicity_EvwTrigger->Fill(lPercentiles);
  fHistTrack->AddBinContent(8, NumberFirstParticleMC);
  fHistTriggerCompositionMCTruth->Fill(TriggerPdgCode,1,ptTriggerMassimoMC);

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
  Float_t kctau[2]={0, 0};
  Int_t ParticleType=0;
  Int_t PDGCodeAssoc[3] = {310, 3122 , 3312};

  if(fV0=="K0s") ParticleType =0;
  if(fV0=="Lambda") ParticleType=1;
  if(fV0=="Xi") ParticleType=2;
  
  //  cout << " beginning loop for associated particles (MCtruth) " << endl;
  //begin MC truth loop for K0s particles as associated 
  if (fMCEvent){
    //	cout << " MC truth for associated particles " << endl;
    //  cout << "number of trakcs" << fMCEvent->GetNumberOfTracks() << endl;
    for ( Int_t i = 0; i< fMCEvent->GetNumberOfTracks(); i++ )    {
      Bool_t skipV0_MC=kFALSE;
      TParticle* particleV0 = static_cast<TParticle*>( fMCEvent->Particle(i) );
      if (!particleV0) continue;
      if (!particleV0->GetPDG()) continue;
      if(TMath::Abs(particleV0->GetPdgCode())!=PDGCodeAssoc[ParticleType]) continue;
      //selection on eta K0s applied online to reduce tree size
      if(TMath::Abs(particleV0->Eta())> fEtaV0Assoc)continue;
      if (!fMCEvent->IsPhysicalPrimary(i))           continue;

      Int_t PAntiP =1;
      if (ParticleType==2){
	PAntiP = particleV0->GetPDG()->Charge()/3.;
      }

      fHistGeneratedV0PtTMaxPhi->Fill(PAntiP*ptTriggerMassimoMC,particleV0->Phi(), lPercentiles );
      fHistGeneratedV0PtTMaxEta->Fill(PAntiP*ptTriggerMassimoMC,particleV0->Eta(), lPercentiles );
      fHistGeneratedV0PtPtTMax->Fill(particleV0->Pt(),PAntiP*ptTriggerMassimoMC, lPercentiles );
      fHistGeneratedV0PtEta->Fill(particleV0->Pt(), particleV0->Eta(), lPercentiles );
      fHistGeneratedV0Pt->Fill(particleV0->Pt(), lPercentiles);

      if(!(particleV0->Pt()> fminPtV0 && particleV0->Pt()<fmaxPtV0) )continue;
	
      if ((particleV0->Pt())>=ptTriggerMassimoMC) {
	skipV0_MC=kTRUE;
	NumberSecondParticleMCNoAssoc++;
      }

      NumberSecondParticleMC++;

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
      fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sCharge       = particleV0->GetPDG()->Charge()/3.;
      fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDCAz         = 0;
      fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDCAxy        = 0;
      fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sAssocOrNot   = skipV0_MC;
      fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sIsCommonParton        = 0;
      fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPdgCommonParton       = 0;
    }
    //	  cout << "I found " << 	NumberSecondParticleMC - 	    NumberSecondParticleMCNoAssoc<< " associated particles (pT < pt,Trig)" << endl;
    //    cout << "I found " << 	NumberSecondParticleMC << " associated particles (all pt)" << endl;
  }
  fHistEventV0->AddBinContent(20, NumberSecondParticle);    
  fHistEventV0->AddBinContent(21, NumberSecondParticleRecoTrue);    
  fHistEventV0->AddBinContent(22, NumberSecondParticleMC);  
  fHistEventV0->AddBinContent(23, NumberSecondParticleNoAssoc);    
  fHistEventV0->AddBinContent(24, NumberSecondParticleMCNoAssoc);      
  fHistMultvsV0All->Fill(NumberSecondParticle,lPercentiles);
  fHistMultvsV0AllTruth->Fill(NumberSecondParticleRecoTrue,lPercentiles);
  fHistMultvsV0MCAll->Fill(NumberSecondParticleMC,lPercentiles);
  fHistSecondParticleAll->Fill(NumberSecondParticle, NumberSecondParticleMC);
  fHistSecondParticleTruthAll->Fill(NumberSecondParticleRecoTrue, NumberSecondParticleMC);
  fHistTriggerNotLeading->Fill(NumberSecondParticleNoAssoc,lPercentiles, ptTriggerMassimoDati);
  fHistTriggerNotLeadingMC->Fill(NumberSecondParticleMCNoAssoc,lPercentiles, ptTriggerMassimoMC);

  NumberSecondParticleAll=NumberSecondParticleMC;
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

  fHistEventV0->AddBinContent(25, NumberSecondParticle);    
  fHistEventV0->AddBinContent(26, NumberSecondParticleRecoTrue);    
  fHistEventV0->AddBinContent(27, NumberSecondParticleMC);    
  
  fHistEventMult->Fill(8);  
  fEvt->fNumberCandidateFirst = NumberFirstParticleMC;
  fEvt->fNumberCandidateSecond = NumberSecondParticleMC;

  //Fill histos for selected events (events with which I perform the same-event and the mixed-event angular correlation, that is having NT> 0 and NV0>0
  fHistTriggerwV0->Fill(NumberFirstParticle);
  fHistTriggerwV0MCTruth->Fill(NumberFirstParticleMC);
  fHist_multiplicity->Fill(lPercentiles);
  fHistZvertex->Fill(lBestPrimaryVtxPos[2]);
  fHistTrack->AddBinContent(9, NumberFirstParticleMC);
  fHistMultvsTrigger->Fill(NumberFirstParticle, lPercentiles);
  fHistMultvsTriggerMCTruth->Fill(NumberFirstParticleMC, lPercentiles);
  fHistMultiplicityVsVertexZ->Fill(lBestPrimaryVtxPos[2], lPercentiles);
 
  for(Int_t i=0; i< 100; i++){
    if (fisHM){
      if(( lPercentiles <((float)i+1)/1000 ) && (lPercentiles >= (float)i/1000) ) fHistTriggervsMult->AddBinContent(i+1,NumberFirstParticle);  
    }
    else {
      if(( lPercentiles <i+1) && (lPercentiles >= i) ) fHistTriggervsMult->AddBinContent(i+1,NumberFirstParticle);  
    }
  }
  for(Int_t i=0; i< 100; i++){
    if (fisHM){
      if(( lPercentiles <((float)i+1)/1000 ) && (lPercentiles >= (float)i/1000) ) fHistTriggervsMultMC->AddBinContent(i+1,NumberFirstParticleMC);  
    }
    else{
      if(( lPercentiles <i+1) && (lPercentiles >= i) ) fHistTriggervsMultMC->AddBinContent(i+1,NumberFirstParticleMC);  
    }
  }


  for(Int_t l=0; l< NumberFirstParticleAll; l++){
    fHistPt->Fill(fEvt->fReconstructedFirst[l].fPt);
    fHistPtvsMult->Fill(fEvt->fReconstructedFirst[l].fPt, lPercentiles);
    fHist_eta_phi->Fill(fEvt->fReconstructedFirst[l].fPhi, fEvt->fReconstructedFirst[l].fEta);
 
  }

  fHistPtMaxvsMult->Fill(ptTriggerMassimoMC, lPercentiles);
  fHist_eta_phi_PtMax->Fill(phiTriggerMassimoMC, etaTriggerMassimoMC);
  
  //--------------------------------------------------------------
  Double_t ptTriggerMassimo=0;
  ptTriggerMassimo=ptTriggerMassimoMC;

  /*
    cout << "\npttrigger massimo dati (per hybrid, data, MC reco) " << ptTriggerMassimoDati << endl;
    cout << "pttrigger massimo mc " << ptTriggerMassimoMC << endl;
    cout << " number candidates first " << fEvt->fNumberCandidateFirst << endl;
    cout << " number first particle " <<   NumberFirstParticleAll << " mc " << NumberFirstParticleMC << " data " << NumberFirstParticle << endl;
    cout << " number candidates second " << fEvt->fNumberCandidateSecond << endl;
    cout << " number first particle " <<   NumberSecondParticleAll << " mc " << NumberSecondParticleMC << " data " << NumberSecondParticle << endl;
  */

  DoPairsh1h2(lPercentiles, lNchVZEROA, lNchVZEROC, fieldsign, lBestPrimaryVtxPos[2], ptTriggerMassimo);  

  PostData(1, fOutputList);     
  PostData(2, fSignalTree);
  PostData(3, fBkgTree);
  PostData(4, fOutputList2);  
  PostData(5, fOutputList3);     
  PostData(6, fOutputList4);     

  fEventColl[zBin][centralityBin]->FifoShift();       
}


//----------------------------------------------------------------------------------------------------

void AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::DoPairsh1h2 ( const Float_t lPercentiles, const Float_t lNchVZEROA, const Float_t lNchVZEROC, int fieldsign, Double_t lBestPrimaryVtxPos, Double_t ptTriggerMassimo)  {

  double dphi  = -999.;
  double dphis = -999.;
  double deta  = -999.;
  double dtheta = -999.;
  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
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
	  fTreeVariablePDGCodeTrigger         = fEvt->fReconstructedFirst[i].fPDGcode;
	  fTreeVariableLabelTrigger           = fEvt->fReconstructedFirst[i].fLabel ;
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
	  fTreeVariableMultiplicityV0A	      = lNchVZEROA;
	  fTreeVariableMultiplicityV0C	      = lNchVZEROC;
	  fTreeVariableZvertex                = lBestPrimaryVtxPos;
	  fTreeVariablePDGCodeAssoc           = fEvt->fReconstructedSecond[j].sPDGcode;
	  fTreeVariableisPrimaryTrigger       = fEvt->fReconstructedFirst[i].isP;
	  fTreeVariableisPrimaryV0            = fEvt->fReconstructedSecond[j].isP;
	  fTreeVariableSkipAssoc              = fEvt->fReconstructedSecond[j].sAssocOrNot;
	  fTreeVariableIsCommonParton         = fEvt->fReconstructedSecond[j].sIsCommonParton;
          fTreeVariablePdgCommonParton        = fEvt->fReconstructedSecond[j].sPdgCommonParton;
	  fTreeVariableLabelPos               = fEvt->fReconstructedSecond[j].sLabelPos;
          fTreeVariableLabelNeg               = fEvt->fReconstructedSecond[j].sLabelNeg;
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
	  fTreeVariableLabelTrigger           = fEvt->fReconstructedFirst[i].fLabel ;
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
	  fTreeVariableisPrimaryV0            = (fEvt+eventNumber)->fReconstructedSecond[j].isP;
	  fTreeVariableSkipAssoc              = (fEvt+eventNumber)->fReconstructedSecond[j].sAssocOrNot ;
	  //	  fTreeVariableEtaDaughterPos         =  (fEvt+eventNumber)->fReconstructedSecond[j].sEtaPos ;
	  //	  fTreeVariableEtaDaughterNeg         =  (fEvt+eventNumber)->fReconstructedSecond[j].sEtaNeg ;
	  fTreeVariableIsCommonParton         = (fEvt+eventNumber)->fReconstructedSecond[j].sIsCommonParton;
          fTreeVariablePdgCommonParton        = (fEvt+eventNumber)->fReconstructedSecond[j].sPdgCommonParton;
	  fTreeVariableLabelPos               = (fEvt+eventNumber)->fReconstructedSecond[j].sLabelPos;
          fTreeVariableLabelNeg               = (fEvt+eventNumber)->fReconstructedSecond[j].sLabelNeg;

	  fTreeVariableDeltaEta	       	      =deta;  
	  fTreeVariableDeltaPhi		      =dphi;
	  fTreeVariableDeltaTheta             =dtheta;      
	  fTreeVariableMultiplicity	      = lPercentiles;
	  fTreeVariableMultiplicityV0A	      = lNchVZEROA;
	  fTreeVariableMultiplicityV0C	      = lNchVZEROC;
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

//-----------------------------------------------------------------------------------------------
double AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::CalculateDPhiStar(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2,Double_t rad) { //AliFemtoUser/AliFemtoPairCutDetaDphi.h

  const Double_t unit_factor = 0.299792458 / 2.0;
  const Double_t b_field = 0.5006670488586 * magSign;

  Double_t  shift1 = TMath::ASin(unit_factor * chg1 * b_field * rad / ptv1);
  Double_t  shift2 = TMath::ASin(unit_factor * chg2 * b_field * rad / ptv2);

  double dps = (phi1 + shift1) - (phi2 + shift2);
  
  //  dps = TVector2::Phi_mpi_pi(dps); //to be checked

  return dps; //deltaphi;
  
}

//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::CalculateDeltaEta( Double_t eta1, Double_t eta2 )  {   //AliFemtoUser/AliFemtoPairCutDetaDphi.h
  const double deta = eta2 - eta1;
  return deta;
}
//_______________________________________________________________
Double_t AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::CalculateDeltaTheta( Double_t theta1, Double_t theta2 )  {  
  const double dtheta = theta2 - theta1;
  return dtheta;
}
//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::CalculateDeltaPhi( Double_t phi1, Double_t phi2 )  {   //AliFemtoUser/AliFemtoPairCutDetaDphi.h
  const double dphi = phi2 - phi1;
  return dphi;
}

//_____________________________________________________________________________
void AliAnalysisTaskCorrelationhK0sXi_PureMCOnly::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________

