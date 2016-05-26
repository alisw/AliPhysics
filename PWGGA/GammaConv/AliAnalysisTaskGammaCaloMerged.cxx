/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                               *
 * Author: Friederike Bock                                                         *
 * Version 1.0                                                                           *
 *                                                                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.                           *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on conversion photons + calo photons
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliAnalysisTaskGammaCaloMerged.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskGammaCaloMerged)

//________________________________________________________________________
AliAnalysisTaskGammaCaloMerged::AliAnalysisTaskGammaCaloMerged(): AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fHeaderNameList(NULL),
  fOutputContainer(NULL),
  fNClusterCandidates(0),
  fNClusterMergedCandidates(0),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fClusterCutArray(NULL),
  fClusterMergedCutArray(NULL),
  fMesonCutArray(NULL),
  fMesonCuts(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherPtY(NULL),
  fHistoMotherPtAlpha(NULL),
  fHistoMotherPtOpenAngle(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusMergedPtvsM02(NULL),
  fHistoClusMergedPtvsM02Accepted(NULL),
  fHistoClusNCellsPt(NULL),
  fHistoClusMergedNCellsPt(NULL),
  fHistoClusMergedNParticlePt(NULL),
  fHistoClusMergedNCellsAroundPt(NULL),
  fHistoClusMergedNCellsAroundAndInPt(NULL),
  fHistoClusMergedEAroundE(NULL),
  fHistoMCHeaders(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCPi0DalitzPt(NULL),
  fHistoMCPi0DalitzWOWeightPt(NULL),
  fHistoMCPi0DalitzWOEvtWeightPt(NULL),
  fHistoMCEtaDalitzPt(NULL),
  fHistoMCEtaDalitzWOWeightPt(NULL),
  fHistoMCEtaDalitzWOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCPi0DalitzInAccPt(NULL),
  fHistoMCEtaDalitzInAccPt(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCEtaPtJetPt(NULL),
  fHistoTrueClusMergedPtvsM02(NULL),  
  fHistoTrueClusPi0PtvsM02(NULL),
  fHistoTrueClusPrimPi0PtvsM02(NULL),
  fHistoTrueClusSecPi0PtvsM02(NULL),
  fHistoTrueClusSecPi0FromK0sPtvsM02(NULL),
  fHistoTrueClusSecPi0FromLambdaPtvsM02(NULL),
  fHistoTrueClusEtaPtvsM02(NULL),
  fHistoTrueClusMergedPartConvPtvsM02(NULL),
  fHistoTrueClusMergedPartConvELeadPtvsM02(NULL),
  fHistoTrueClusPartConvPi0PtvsM02(NULL),
  fHistoTrueClusPartConvPrimPi0PtvsM02(NULL),
  fHistoTrueClusPartConvSecPi0PtvsM02(NULL),
  fHistoTrueClusPartConvSecPi0FromK0sPtvsM02(NULL),
  fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02(NULL),
  fHistoTrueClusPartConvEtaPtvsM02(NULL),
  fHistoTrueClusPartConvGammaPtvsM02(NULL),
  fHistoTrueClusBGPtvsM02(NULL),
  fHistoTrueClusGammaPtvsM02(NULL),
  fHistoTrueClusGammaFromPi0PtvsM02(NULL),
  fHistoTrueClusGammaFromEtaPtvsM02(NULL),
  fHistoTrueClusElectronPtvsM02(NULL),
  fHistoTrueClusElectronFromPi0PtvsM02(NULL),
  fHistoTrueClusElectronFromEtaPtvsM02(NULL),
  fHistoTrueClusElectronFromGammaPtvsM02(NULL),
  fHistoTrueClusMergedInvMassvsPt(NULL),  
  fHistoTrueClusPi0InvMassvsPt(NULL),
  fHistoTrueClusPrimPi0InvMassvsPt(NULL),
  fHistoTrueClusSecPi0InvMassvsPt(NULL),
  fHistoTrueClusSecPi0FromK0sInvMassvsPt(NULL),
  fHistoTrueClusSecPi0FromLambdaInvMassvsPt(NULL),
  fHistoTrueClusEtaInvMassvsPt(NULL),
  fHistoTrueClusMergedPartConvInvMassvsPt(NULL),
  fHistoTrueClusMergedPartConvELeadInvMassvsPt(NULL),
  fHistoTrueClusPartConvPi0InvMassvsPt(NULL),
  fHistoTrueClusPartConvPrimPi0InvMassvsPt(NULL),
  fHistoTrueClusPartConvSecPi0InvMassvsPt(NULL),
  fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt(NULL),
  fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt(NULL),
  fHistoTrueClusPartConvEtaInvMassvsPt(NULL),
  fHistoTrueClusBGInvMassvsPt(NULL),
  fHistoTrueClusGammaInvMassvsPt(NULL),  
  fHistoTrueClusElectronInvMassvsPt(NULL),  
  fHistoTrueClusBGPtvsSource(NULL),
  fHistoTrueClusGammaPtvsSource(NULL),
  fHistoTrueClusElectronPtvsSource(NULL),
  fHistoTrueMergedMissedPDG(NULL),
  fHistoTrueMergedPartConvMissedPDG(NULL),
  fHistoTrueMergedPartConvNonLeadingPtvsM02(NULL),
  fHistoTruePi0PtY(NULL),
  fHistoTrueEtaPtY(NULL),
  fHistoTruePi0PtAlpha(NULL),
  fHistoTrueEtaPtAlpha(NULL),
  fHistoTruePi0PtOpenAngle(NULL),
  fHistoTrueEtaPtOpenAngle(NULL),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoNClusterCandidates(NULL),
  fHistoNClusterMergedCandidates(NULL),
  fHistoNGoodESDTracksVsNClusterCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fDoMesonQA(0),
  fDoClusterQA(0),
  fIsFromMBHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fSetPlotHistsExtQA(kFALSE),
  fWeightJetJetMC(1),
  fSelectedMesonID(0),
  fEnableDetailedPrintOut(kFALSE)
{
  
}

//________________________________________________________________________
AliAnalysisTaskGammaCaloMerged::AliAnalysisTaskGammaCaloMerged(const char *name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fHeaderNameList(NULL),
  fOutputContainer(NULL),
  fNClusterCandidates(0),
  fNClusterMergedCandidates(0),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fClusterCutArray(NULL),
  fClusterMergedCutArray(NULL),
  fMesonCutArray(NULL),
  fMesonCuts(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherPtY(NULL),
  fHistoMotherPtAlpha(NULL),
  fHistoMotherPtOpenAngle(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusMergedPtvsM02(NULL),
  fHistoClusMergedPtvsM02Accepted(NULL),
  fHistoClusNCellsPt(NULL),
  fHistoClusMergedNCellsPt(NULL),
  fHistoClusMergedNParticlePt(NULL),
  fHistoClusMergedNCellsAroundPt(NULL),
  fHistoClusMergedNCellsAroundAndInPt(NULL),
  fHistoClusMergedEAroundE(NULL),
  fHistoMCHeaders(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCPi0DalitzPt(NULL),
  fHistoMCPi0DalitzWOWeightPt(NULL),
  fHistoMCPi0DalitzWOEvtWeightPt(NULL),
  fHistoMCEtaDalitzPt(NULL),
  fHistoMCEtaDalitzWOWeightPt(NULL),
  fHistoMCEtaDalitzWOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCPi0DalitzInAccPt(NULL),
  fHistoMCEtaDalitzInAccPt(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCEtaPtJetPt(NULL),
  fHistoTrueClusMergedPtvsM02(NULL),  
  fHistoTrueClusPi0PtvsM02(NULL),
  fHistoTrueClusPrimPi0PtvsM02(NULL),
  fHistoTrueClusSecPi0PtvsM02(NULL),
  fHistoTrueClusSecPi0FromK0sPtvsM02(NULL),
  fHistoTrueClusSecPi0FromLambdaPtvsM02(NULL),
  fHistoTrueClusEtaPtvsM02(NULL),
  fHistoTrueClusMergedPartConvPtvsM02(NULL),
  fHistoTrueClusMergedPartConvELeadPtvsM02(NULL),
  fHistoTrueClusPartConvPi0PtvsM02(NULL),
  fHistoTrueClusPartConvPrimPi0PtvsM02(NULL),
  fHistoTrueClusPartConvSecPi0PtvsM02(NULL),
  fHistoTrueClusPartConvSecPi0FromK0sPtvsM02(NULL),
  fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02(NULL),
  fHistoTrueClusPartConvEtaPtvsM02(NULL),
  fHistoTrueClusPartConvGammaPtvsM02(NULL),
  fHistoTrueClusBGPtvsM02(NULL),
  fHistoTrueClusGammaPtvsM02(NULL),
  fHistoTrueClusGammaFromPi0PtvsM02(NULL),
  fHistoTrueClusGammaFromEtaPtvsM02(NULL),
  fHistoTrueClusElectronPtvsM02(NULL),
  fHistoTrueClusElectronFromPi0PtvsM02(NULL),
  fHistoTrueClusElectronFromEtaPtvsM02(NULL),
  fHistoTrueClusElectronFromGammaPtvsM02(NULL),
  fHistoTrueClusMergedInvMassvsPt(NULL),  
  fHistoTrueClusPi0InvMassvsPt(NULL),
  fHistoTrueClusPrimPi0InvMassvsPt(NULL),
  fHistoTrueClusSecPi0InvMassvsPt(NULL),
  fHistoTrueClusSecPi0FromK0sInvMassvsPt(NULL),
  fHistoTrueClusSecPi0FromLambdaInvMassvsPt(NULL),
  fHistoTrueClusEtaInvMassvsPt(NULL),
  fHistoTrueClusMergedPartConvInvMassvsPt(NULL),
  fHistoTrueClusMergedPartConvELeadInvMassvsPt(NULL),
  fHistoTrueClusPartConvPi0InvMassvsPt(NULL),
  fHistoTrueClusPartConvPrimPi0InvMassvsPt(NULL),
  fHistoTrueClusPartConvSecPi0InvMassvsPt(NULL),
  fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt(NULL),
  fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt(NULL),
  fHistoTrueClusPartConvEtaInvMassvsPt(NULL),
  fHistoTrueClusBGInvMassvsPt(NULL),
  fHistoTrueClusGammaInvMassvsPt(NULL),
  fHistoTrueClusElectronInvMassvsPt(NULL),
  fHistoTrueClusBGPtvsSource(NULL),
  fHistoTrueClusGammaPtvsSource(NULL),
  fHistoTrueClusElectronPtvsSource(NULL),
  fHistoTrueMergedMissedPDG(NULL),
  fHistoTrueMergedPartConvMissedPDG(NULL),
  fHistoTrueMergedPartConvNonLeadingPtvsM02(NULL),
  fHistoTruePi0PtY(NULL),
  fHistoTrueEtaPtY(NULL),
  fHistoTruePi0PtAlpha(NULL),
  fHistoTrueEtaPtAlpha(NULL),
  fHistoTruePi0PtOpenAngle(NULL),
  fHistoTrueEtaPtOpenAngle(NULL),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoNClusterCandidates(NULL),
  fHistoNClusterMergedCandidates(NULL),
  fHistoNGoodESDTracksVsNClusterCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fDoMesonQA(0),
  fDoClusterQA(0),
  fIsFromMBHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fSetPlotHistsExtQA(kFALSE),
  fWeightJetJetMC(1),
  fSelectedMesonID(0),
  fEnableDetailedPrintOut(kFALSE)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaCaloMerged::~AliAnalysisTaskGammaCaloMerged()
{

}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::UserCreateOutputObjects(){
  
  Int_t invMassBins                           = 800;
  Float_t startMass                           = 0;
  Float_t endMass                             = 0.8;
  
  if (GetSelectedMesonID() == 1) {
    invMassBins                               = 400;
    startMass                                 = 0;
    endMass                                   = 0.4;
  } else if (GetSelectedMesonID() == 2) {
    invMassBins                               = 800;
    startMass                                 = 0.;
    endMass                                   = 0.8;    
  }
    
    
  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer                          = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer                          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
    
  fCutFolder                                  = new TList*[fnCuts];
  fESDList                                    = new TList*[fnCuts];

  fHistoNEvents                               = new TH1F*[fnCuts];
  if(fIsMC == 2){
    fHistoNEventsWOWeight                     = new TH1F*[fnCuts];
    fProfileJetJetXSection                    = new TProfile*[fnCuts];
    fHistoJetJetNTrials                       = new TH1F*[fnCuts];
  }
  
  fHistoNGoodESDTracks                        = new TH1F*[fnCuts];
  fHistoVertexZ                               = new TH1F*[fnCuts];
  fHistoNClusterCandidates                    = new TH1F*[fnCuts];
  fHistoNClusterMergedCandidates              = new TH1F*[fnCuts];
  fHistoNGoodESDTracksVsNClusterCandidates    = new TH2F*[fnCuts];
  fHistoSPDClusterTrackletBackground          = new TH2F*[fnCuts];
  fHistoNV0Tracks                             = new TH1F*[fnCuts];
  fProfileEtaShift                            = new TProfile*[fnCuts];
    
  fHistoMotherInvMassPt                       = new TH2F*[fnCuts];
  if (fDoMesonQA > 0 ){
    fHistoMotherPtY                           =  new TH2F*[fnCuts];
    fHistoMotherPtAlpha                       =  new TH2F*[fnCuts];
    fHistoMotherPtOpenAngle                   = new TH2F*[fnCuts];
  }
    
  fHistoClusGammaPt                           = new TH1F*[fnCuts];
  fHistoClusOverlapHeadersGammaPt             = new TH1F*[fnCuts];
  fHistoClusMergedPtvsM02                     = new TH2F*[fnCuts];
  fHistoClusMergedPtvsM02Accepted             = new TH2F*[fnCuts];
  if (fDoClusterQA > 0){
    fHistoClusNCellsPt                        = new TH2F*[fnCuts];
    fHistoClusMergedNCellsPt                  = new TH2F*[fnCuts];
    fHistoClusMergedNCellsAroundPt            = new TH2F*[fnCuts];
    fHistoClusMergedNCellsAroundAndInPt       = new TH2F*[fnCuts];
    fHistoClusMergedEAroundE                  = new TH2F*[fnCuts];

    if (fIsMC > 0){
      fHistoClusMergedNParticlePt             = new TH2F*[fnCuts];
    }
    
  }  
  
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent                        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringCalo                         = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringCaloMerged                   = ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson                        = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
    
    Int_t nLMCut   = ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetMinNLMCut();
    
    fCutFolder[iCut]                              = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s_%s",cutstringEvent.Data(), cutstringCalo.Data(), cutstringCaloMerged.Data(), cutstringMeson.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut]                                = new TList();
    fESDList[iCut]->SetName(Form("%s_%s_%s_%s ESD histograms",cutstringEvent.Data(), cutstringCalo.Data(), cutstringCaloMerged.Data(), cutstringMeson.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);
    
    fHistoNEvents[iCut]                           = new TH1F("NEvents","NEvents",12,-0.5,11.5);
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){ 
      TString TriggerNames = "Not Trigger: ";
      TriggerNames = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    } else {
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problems");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);
    
    if (fIsMC == 2){
      fHistoNEventsWOWeight[iCut]                 = new TH1F("NEventsWOWeight","NEventsWOWeight",12,-0.5,11.5);
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
      if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
        TString TriggerNames = "Not Trigger: ";
        TriggerNames = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
        fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
      } else {
        fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
      }
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
      fESDList[iCut]->Add(fHistoNEventsWOWeight[iCut]);

      fProfileJetJetXSection[iCut]                = new TProfile("XSection","XSection",1,-0.5,0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fHistoJetJetNTrials[iCut]                   = new TH1F("NTrials","#sum{NTrials}",1,0,1);
      fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
    }

    if(fIsHeavyIon == 1) 
      fHistoNGoodESDTracks[iCut]                  = new TH1F("GoodESDTracks","GoodESDTracks",4000,-0.5,3999.5);
    else if(fIsHeavyIon == 2) 
      fHistoNGoodESDTracks[iCut]                  = new TH1F("GoodESDTracks","GoodESDTracks",400,-0.5,399.5);
    else 
      fHistoNGoodESDTracks[iCut]                  = new TH1F("GoodESDTracks","GoodESDTracks",200,-0.5,199.5);
    fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

    fHistoVertexZ[iCut]                           = new TH1F("VertexZ","VertexZ",1000,-50,50);
    fESDList[iCut]->Add(fHistoVertexZ[iCut]);

    if(fIsHeavyIon == 1) 
      fHistoNClusterCandidates[iCut]              = new TH1F("GammaCandidates","GammaCandidates",600,-0.5,599.5);
    else if(fIsHeavyIon == 2) 
      fHistoNClusterCandidates[iCut]              = new TH1F("GammaCandidates","GammaCandidates",400,-0.5,399.5);
    else 
      fHistoNClusterCandidates[iCut]              = new TH1F("GammaCandidates","GammaCandidates",100,-0.5,99.5);
    fESDList[iCut]->Add(fHistoNClusterCandidates[iCut]);

    if(fIsHeavyIon == 1) 
      fHistoNClusterMergedCandidates[iCut]        = new TH1F("MergedCandidates","MergedCandidates",300,-0.5,299.5);
    else if(fIsHeavyIon == 2) 
      fHistoNClusterMergedCandidates[iCut]        = new TH1F("MergedCandidates","MergedCandidates",200,-0.5,199.5);
    else 
      fHistoNClusterMergedCandidates[iCut]        = new TH1F("MergedCandidates","MergedCandidates",50,-0.5,49.5);
    fESDList[iCut]->Add(fHistoNClusterMergedCandidates[iCut]);

    
    if(fIsHeavyIon == 1) 
      fHistoNGoodESDTracksVsNClusterCandidates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",4000,-0.5,3999.5,600,-0.5,599.5);
    else if(fIsHeavyIon == 2) 
      fHistoNGoodESDTracksVsNClusterCandidates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",400,-0.5,399.5,400,-0.5,399.5);
    else 
      fHistoNGoodESDTracksVsNClusterCandidates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",200,-0.5,199.5,100,-0.5,99.5);
    fESDList[iCut]->Add(fHistoNGoodESDTracksVsNClusterCandidates[iCut]);

    fHistoSPDClusterTrackletBackground[iCut]      = new TH2F("SPD tracklets vs SPD clusters","SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
    fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);
    
    if(fIsHeavyIon == 1) 
      fHistoNV0Tracks[iCut]                       = new TH1F("V0 Multiplicity","V0 Multiplicity",30000,-0.5,29999.5);
    else if(fIsHeavyIon == 2) 
      fHistoNV0Tracks[iCut]                       = new TH1F("V0 Multiplicity","V0 Multiplicity",2500,-0.5,2499.5);
    else 
      fHistoNV0Tracks[iCut]                       = new TH1F("V0 Multiplicity","V0 Multiplicity",1500,-0.5,1499.5);
    fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);
    fProfileEtaShift[iCut]                        = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
    fESDList[iCut]->Add(fProfileEtaShift[iCut]);

    if (fIsMC == 2){
      fHistoNEvents[iCut]->Sumw2();
      fHistoNGoodESDTracks[iCut]->Sumw2();
      fHistoVertexZ[iCut]->Sumw2();
      fHistoNClusterCandidates[iCut]->Sumw2();
      fHistoNGoodESDTracksVsNClusterCandidates[iCut]->Sumw2();
      fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
      fHistoNV0Tracks[iCut]->Sumw2();
      fProfileEtaShift[iCut]->Sumw2();
    }

    fHistoClusGammaPt[iCut]                       = new TH1F("ClusGamma_Pt","ClusGamma_Pt",500, 0, 50);
    fESDList[iCut]->Add(fHistoClusGammaPt[iCut]);
    fHistoClusOverlapHeadersGammaPt[iCut]         = new TH1F("ClusGammaOverlapHeaders_Pt","ClusGammaOverlapHeaders_Pt",500, 0, 50);
    fESDList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);
    fHistoClusMergedPtvsM02[iCut]                 = new TH2F("ClusMerged_Pt_M02","ClusMerged_Pt_M02",500,0,50,500,0,5);
    fESDList[iCut]->Add(fHistoClusMergedPtvsM02[iCut]);
    fHistoClusMergedPtvsM02Accepted[iCut]         = new TH2F("ClusMerged_Pt_M02_AcceptedMeson","ClusMerged_Pt_M02_AcceptedMeson",500,0,50,500,0,5);
    fESDList[iCut]->Add(fHistoClusMergedPtvsM02Accepted[iCut]);

    if (fIsMC == 2){
      fHistoClusGammaPt[iCut]->Sumw2();
      fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
      fHistoClusMergedPtvsM02[iCut]->Sumw2();
      fHistoClusMergedPtvsM02Accepted[iCut]->Sumw2();
    }
    
    if (fDoClusterQA > 0){
      fHistoClusNCellsPt[iCut]                      = new TH2F("Clus_NCells_Pt","Clus_NCells_Pt",100,-0.5,99.5,500,0,50);
      fESDList[iCut]->Add(fHistoClusNCellsPt[iCut]);
      fHistoClusMergedNCellsPt[iCut]                = new TH2F("ClusMerged_NCells_Pt","ClusMerged_NCells_Pt",100,-0.5,99.5,500,0,50);
      fESDList[iCut]->Add(fHistoClusMergedNCellsPt[iCut]);
      fHistoClusMergedNCellsAroundPt[iCut]          = new TH2F("ClusMerged_NCellsAroundClus_Pt","ClusMerged_NCellsAroundClus_Pt",100,-0.5,99.5,500,0,50);
      fESDList[iCut]->Add(fHistoClusMergedNCellsAroundPt[iCut]);
      fHistoClusMergedNCellsAroundAndInPt[iCut]     = new TH2F("ClusMerged_NCellsAroundAndInClus_Pt","ClusMerged_NCellsAroundAndInClus_Pt",100,-0.5,99.5,500,0,50);
      fESDList[iCut]->Add(fHistoClusMergedNCellsAroundAndInPt[iCut]);
      fHistoClusMergedEAroundE[iCut]                = new TH2F("ClusMerged_EAroundClus_E","ClusMerged_EAroundClus_E",500,0,100,500,0,50);
      fESDList[iCut]->Add(fHistoClusMergedEAroundE[iCut]);
  
      if (fIsMC > 0){
        fHistoClusMergedNParticlePt[iCut]           = new TH2F("ClusMerged_NPart_Pt","ClusMerged_NPart_Pt",100,-0.5,99.5,500,0,50);
        fESDList[iCut]->Add(fHistoClusMergedNParticlePt[iCut]);
        if (fIsMC == 2){
          fHistoClusNCellsPt[iCut]->Sumw2();
          fHistoClusMergedNCellsPt[iCut]->Sumw2();
          fHistoClusMergedNCellsAroundPt[iCut]->Sumw2();
          fHistoClusMergedNCellsAroundAndInPt[iCut]->Sumw2();
          fHistoClusMergedEAroundE[iCut]->Sumw2();
          fHistoClusMergedNParticlePt[iCut]->Sumw2();
          fHistoClusMergedEAroundE[iCut]->Sumw2();
        }  
      }
    }  
  
    fHistoMotherInvMassPt[iCut]                   = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
    fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);

    if (fIsMC == 2){
      fHistoMotherInvMassPt[iCut]->Sumw2();
    }
    
    if (fDoMesonQA > 0 ){
      fHistoMotherPtY[iCut]                       = new TH2F("ESD_Mother_Pt_Y","ESD_Mother_Pt_Y",350,0.03,35.,150,-1.5,1.5);
      SetLogBinningXTH2(fHistoMotherPtY[iCut]);
      fESDList[iCut]->Add(fHistoMotherPtY[iCut]);
      fHistoMotherPtAlpha[iCut]                   = new TH2F("ESD_Mother_Pt_Alpha","ESD_Mother_Pt_Alpha",350,0.03,35.,100,0,1);
      SetLogBinningXTH2(fHistoMotherPtAlpha[iCut]);
      fESDList[iCut]->Add(fHistoMotherPtAlpha[iCut]);
      fHistoMotherPtOpenAngle[iCut]               = new TH2F("ESD_Mother_Pt_OpenAngle","ESD_Mother_Pt_OpenAngle",350,0.03,35.,100,0, 0.5);
      SetLogBinningXTH2(fHistoMotherPtOpenAngle[iCut]);
      fESDList[iCut]->Add(fHistoMotherPtOpenAngle[iCut]);
      if (fIsMC == 2){
        fHistoMotherPtY[iCut]->Sumw2();
        fHistoMotherPtAlpha[iCut]->Sumw2();
        fHistoMotherPtOpenAngle[iCut]->Sumw2();
      }
    }
  }
  
  if (fIsMC > 0){
    fMCList = new TList*[fnCuts];
    fTrueList = new TList*[fnCuts];

    if (GetSelectedMesonID() != 2){
      fHistoMCPi0Pt                               = new TH1F*[fnCuts];
      fHistoMCPi0WOWeightPt                       = new TH1F*[fnCuts];
      fHistoMCPi0InAccPt                          = new TH1F*[fnCuts];
      fHistoMCPi0DalitzPt                         = new TH1F*[fnCuts];
      fHistoMCPi0DalitzWOWeightPt                 = new TH1F*[fnCuts];
      fHistoMCPi0DalitzInAccPt                    = new TH1F*[fnCuts];
    }
    if (GetSelectedMesonID() != 1){
      fHistoMCEtaPt                               = new TH1F*[fnCuts];
      fHistoMCEtaWOWeightPt                       = new TH1F*[fnCuts];
      fHistoMCEtaInAccPt                          = new TH1F*[fnCuts];
      fHistoMCEtaDalitzPt                         = new TH1F*[fnCuts];
      fHistoMCEtaDalitzWOWeightPt                 = new TH1F*[fnCuts];
      fHistoMCEtaDalitzInAccPt                    = new TH1F*[fnCuts];
    }
    if (fIsMC == 2){
      if (GetSelectedMesonID() != 2)  
        fHistoMCPi0WOEvtWeightPt                  = new TH1F*[fnCuts];
        fHistoMCPi0DalitzWOEvtWeightPt            = new TH1F*[fnCuts];
      if (GetSelectedMesonID() != 1)
        fHistoMCEtaWOEvtWeightPt                  = new TH1F*[fnCuts];
        fHistoMCEtaDalitzWOEvtWeightPt            = new TH1F*[fnCuts];
      if (fDoMesonQA > 0){
        if (GetSelectedMesonID() != 2)
          fHistoMCPi0PtJetPt                      = new TH2F*[fnCuts];
        if (GetSelectedMesonID() != 1)
          fHistoMCEtaPtJetPt                      = new TH2F*[fnCuts];
      }
    }

    fHistoTrueClusMergedPtvsM02                   = new TH2F*[fnCuts];
    fHistoTrueClusPi0PtvsM02                      = new TH2F*[fnCuts];
    if (GetSelectedMesonID() < 2){
      fHistoTrueClusPrimPi0PtvsM02                = new TH2F*[fnCuts];
      fHistoTrueClusSecPi0PtvsM02                 = new TH2F*[fnCuts];
      fHistoTrueClusSecPi0FromK0sPtvsM02          = new TH2F*[fnCuts];
      fHistoTrueClusSecPi0FromLambdaPtvsM02       = new TH2F*[fnCuts];
    }
    fHistoTrueClusEtaPtvsM02                      = new TH2F*[fnCuts];
    fHistoTrueClusMergedPartConvPtvsM02           = new TH2F*[fnCuts];
    fHistoTrueClusMergedPartConvELeadPtvsM02      = new TH2F*[fnCuts];
    fHistoTrueClusPartConvPi0PtvsM02              = new TH2F*[fnCuts];
    fHistoTrueClusPartConvEtaPtvsM02              = new TH2F*[fnCuts];
    fHistoTrueClusPartConvGammaPtvsM02            = new TH2F*[fnCuts];
    if (GetSelectedMesonID() < 2){
      fHistoTrueClusPartConvPrimPi0PtvsM02          = new TH2F*[fnCuts];
      fHistoTrueClusPartConvSecPi0PtvsM02           = new TH2F*[fnCuts];
      fHistoTrueClusPartConvSecPi0FromK0sPtvsM02    = new TH2F*[fnCuts];
      fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02 = new TH2F*[fnCuts];
    }
    fHistoTrueClusBGPtvsM02                       = new TH2F*[fnCuts];
    fHistoTrueClusGammaPtvsM02                    = new TH2F*[fnCuts];
    fHistoTrueClusGammaFromPi0PtvsM02             = new TH2F*[fnCuts];
    fHistoTrueClusGammaFromEtaPtvsM02             = new TH2F*[fnCuts];
    fHistoTrueClusElectronPtvsM02                 = new TH2F*[fnCuts];
    fHistoTrueClusElectronFromPi0PtvsM02          = new TH2F*[fnCuts];
    fHistoTrueClusElectronFromEtaPtvsM02          = new TH2F*[fnCuts];
    fHistoTrueClusElectronFromGammaPtvsM02        = new TH2F*[fnCuts];
    
    fHistoTrueClusBGPtvsSource                    = new TH2F*[fnCuts];
    fHistoTrueClusGammaPtvsSource                 = new TH2F*[fnCuts];
    fHistoTrueClusElectronPtvsSource              = new TH2F*[fnCuts];
    fHistoTrueMergedMissedPDG                     = new TH1F*[fnCuts];
    fHistoTrueMergedPartConvMissedPDG             = new TH1F*[fnCuts];
    fHistoTrueMergedPartConvNonLeadingPtvsM02     = new TH2F*[fnCuts];

    if (fDoMesonQA > 1){
      fHistoTrueClusMergedInvMassvsPt               = new TH2F*[fnCuts];
      fHistoTrueClusPi0InvMassvsPt                  = new TH2F*[fnCuts];
      fHistoTrueClusEtaInvMassvsPt                  = new TH2F*[fnCuts];
      fHistoTrueClusMergedPartConvInvMassvsPt       = new TH2F*[fnCuts];
      fHistoTrueClusPartConvPi0InvMassvsPt          = new TH2F*[fnCuts];
      fHistoTrueClusPartConvEtaInvMassvsPt          = new TH2F*[fnCuts];
      fHistoTrueClusBGInvMassvsPt                   = new TH2F*[fnCuts];
      fHistoTrueClusGammaInvMassvsPt                = new TH2F*[fnCuts];
      fHistoTrueClusElectronInvMassvsPt             = new TH2F*[fnCuts];
      if (GetSelectedMesonID() < 2){
        fHistoTrueClusPrimPi0InvMassvsPt                  = new TH2F*[fnCuts];
        fHistoTrueClusSecPi0InvMassvsPt                   = new TH2F*[fnCuts];
        fHistoTrueClusSecPi0FromK0sInvMassvsPt            = new TH2F*[fnCuts];
        fHistoTrueClusSecPi0FromLambdaInvMassvsPt         = new TH2F*[fnCuts];                                    
        fHistoTrueClusPartConvPrimPi0InvMassvsPt          = new TH2F*[fnCuts];
        fHistoTrueClusPartConvSecPi0InvMassvsPt           = new TH2F*[fnCuts];
        fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt    = new TH2F*[fnCuts];
        fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt = new TH2F*[fnCuts];
      }
      fHistoTrueClusMergedPartConvELeadInvMassvsPt= new TH2F*[fnCuts];
    }
    
    if (fDoMesonQA > 0){
      if (GetSelectedMesonID() != 2){
        fHistoTruePi0PtY                          = new TH2F*[fnCuts];
        fHistoTruePi0PtAlpha                      = new TH2F*[fnCuts];
        fHistoTruePi0PtOpenAngle                  = new TH2F*[fnCuts];
      }
      if (GetSelectedMesonID() != 1){
        fHistoTrueEtaPtY                          = new TH2F*[fnCuts];
        fHistoTrueEtaPtAlpha                      = new TH2F*[fnCuts];
        fHistoTrueEtaPtOpenAngle                  = new TH2F*[fnCuts];                  
      }
    }
    
    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent                        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo                         = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringCaloMerged                   = ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson                        = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
      
      
      fMCList[iCut]                                 = new TList();
      fMCList[iCut]->SetName(Form("%s_%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringCaloMerged.Data(), cutstringMeson.Data()));
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);

      if (GetSelectedMesonID() != 2){
        fHistoMCPi0Pt[iCut]                         = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",500, 0, 50);
        fHistoMCPi0Pt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
        fHistoMCPi0WOWeightPt[iCut]                 = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",500, 0, 50);
        fHistoMCPi0WOWeightPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);
        fHistoMCPi0InAccPt[iCut]                    = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",500, 0, 50);
        fHistoMCPi0InAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
        fHistoMCPi0DalitzPt[iCut]                   = new TH1F("MC_Pi0Dalitz_Pt","MC_Pi0Dalitz_Pt",500, 0, 50);
        fHistoMCPi0DalitzPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0DalitzPt[iCut]);
        fHistoMCPi0DalitzWOWeightPt[iCut]           = new TH1F("MC_Pi0Dalitz_WOWeights_Pt","MC_Pi0Dalitz_WOWeights_Pt",500, 0, 50);
        fHistoMCPi0DalitzWOWeightPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0DalitzWOWeightPt[iCut]);
        fHistoMCPi0DalitzInAccPt[iCut]              = new TH1F("MC_Pi0DalitzInAcc_Pt","MC_Pi0DalitzInAcc_Pt",500, 0, 50);
        fHistoMCPi0DalitzInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0DalitzInAccPt[iCut]);
      }
      if (GetSelectedMesonID() != 1){
        fHistoMCEtaPt[iCut]                         = new TH1F("MC_Eta_Pt","MC_Eta_Pt",500, 0, 50);
        fHistoMCEtaPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaPt[iCut]);
        fHistoMCEtaWOWeightPt[iCut]                 = new TH1F("MC_Eta_WOWeights_Pt","MC_Eta_WOWeights_Pt",500, 0, 50);
        fHistoMCEtaWOWeightPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaWOWeightPt[iCut]);
        fHistoMCEtaInAccPt[iCut]                    = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",500, 0, 50);
        fHistoMCEtaInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);
        fHistoMCEtaDalitzPt[iCut]                   = new TH1F("MC_EtaDalitz_Pt","MC_EtaDalitz_Pt",500, 0, 50);
        fHistoMCEtaDalitzPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaDalitzPt[iCut]);
        fHistoMCEtaDalitzWOWeightPt[iCut]           = new TH1F("MC_EtaDalitz_WOWeights_Pt","MC_EtaDalitz_WOWeights_Pt",500, 0, 50);
        fHistoMCEtaDalitzWOWeightPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaDalitzWOWeightPt[iCut]);
        fHistoMCEtaDalitzInAccPt[iCut]              = new TH1F("MC_EtaDalitzInAcc_Pt","MC_EtaDalitzInAcc_Pt",500, 0, 50);
        fHistoMCEtaDalitzInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaDalitzInAccPt[iCut]);

        
      }
      if (fIsMC == 2){
        if (GetSelectedMesonID() != 2){
          fHistoMCPi0WOWeightPt[iCut]->Sumw2();
          fHistoMCPi0WOEvtWeightPt[iCut]            = new TH1F("MC_Pi0_WOEventWeights_Pt","MC_Pi0_WOEventWeights_Pt",500, 0, 50);
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);
          fHistoMCPi0DalitzWOWeightPt[iCut]->Sumw2();
          fHistoMCPi0DalitzWOEvtWeightPt[iCut]      = new TH1F("MC_Pi0Dalitz_WOEventWeights_Pt","MC_Pi0Dalitz_WOEventWeights_Pt",500, 0, 50);
          fMCList[iCut]->Add(fHistoMCPi0DalitzWOEvtWeightPt[iCut]);

          
        }
        if (GetSelectedMesonID() != 1){
          fHistoMCEtaWOWeightPt[iCut]->Sumw2();
          fHistoMCEtaWOEvtWeightPt[iCut]            = new TH1F("MC_Eta_WOEventWeights_Pt","MC_Eta_WOEventWeights_Pt",500, 0, 50);
          fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightPt[iCut]);
          fHistoMCEtaDalitzWOWeightPt[iCut]->Sumw2();
          fHistoMCEtaDalitzWOEvtWeightPt[iCut]      = new TH1F("MC_EtaDalitz_WOEventWeights_Pt","MC_EtaDalitz_WOEventWeights_Pt",500, 0, 50);
          fMCList[iCut]->Add(fHistoMCEtaDalitzWOEvtWeightPt[iCut]);
        }
        if (fDoMesonQA > 0){
          if (GetSelectedMesonID() != 2){
            fHistoMCPi0PtJetPt[iCut]                = new TH2F("MC_Pi0_Pt_JetPt","MC_Pi0_Pt_JetPt",350,0.03,35.,200,-0.5,199.5);
            fHistoMCPi0PtJetPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoMCPi0PtJetPt[iCut]);
            fMCList[iCut]->Add(fHistoMCPi0PtJetPt[iCut]);
          }
          if (GetSelectedMesonID() != 1){
            fHistoMCEtaPtJetPt[iCut]                = new TH2F("MC_Eta_Pt_JetPt","MC_Eta_Pt_JetPt",350,0.03,35.,200,-0.5,199.5);
            fHistoMCEtaPtJetPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoMCEtaPtJetPt[iCut]);
            fMCList[iCut]->Add(fHistoMCEtaPtJetPt[iCut]);
          }
        }
      }
        
      fTrueList[iCut]                               = new TList();
      fTrueList[iCut]->SetName(Form("%s_%s_%s_%s True histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringCaloMerged.Data(), cutstringMeson.Data()));
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      fHistoTrueClusMergedPtvsM02[iCut]             = new TH2F("ESD_TrueClusMerged_Pt_M02","ESD_TrueClusMerged_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusMergedPtvsM02[iCut]); 
      fHistoTrueClusPi0PtvsM02 [iCut]               = new TH2F("ESD_TrueClusFromPi0_Pt_M02","ESD_TrueClusFromPi0_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusPi0PtvsM02[iCut]);
      fHistoTrueClusEtaPtvsM02[iCut]                = new TH2F("ESD_TrueClusFromEta_Pt_M02","ESD_TrueClusFromEta_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusEtaPtvsM02[iCut]);                  
      fHistoTrueClusMergedPartConvPtvsM02[iCut]     = new TH2F("ESD_TrueClusMergedPartConv_Pt_M02","ESD_TrueClusMergedPartConv_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvPtvsM02[iCut]);
      fHistoTrueClusMergedPartConvELeadPtvsM02[iCut]= new TH2F("ESD_TrueClusMergedPartConvLeadE_Pt_M02","ESD_TrueClusMergedPartConvLeadE_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvELeadPtvsM02[iCut]);
      fHistoTrueClusPartConvPi0PtvsM02[iCut]        = new TH2F("ESD_TrueClusPartConvFromPi0_Pt_M02","ESD_TrueClusPartConvFromPi0_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusPartConvPi0PtvsM02[iCut]);
      fHistoTrueClusPartConvEtaPtvsM02[iCut]        = new TH2F("ESD_TrueClusPartConvFromEta_Pt_M02","ESD_TrueClusPartConvFromEta_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusPartConvEtaPtvsM02[iCut]);
      fHistoTrueClusPartConvGammaPtvsM02[iCut]      = new TH2F("ESD_TrueClusPartConvFromGamma_Pt_M02","ESD_TrueClusPartConvFromEta_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusPartConvGammaPtvsM02[iCut]);
      fHistoTrueMergedPartConvNonLeadingPtvsM02[iCut]= new TH2F("ESD_TrueClusPartConvNonLeading_Pt_M02","ESD_TrueClusPartConvNonLeading_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueMergedPartConvNonLeadingPtvsM02[iCut]);

      fHistoTrueClusBGPtvsM02[iCut]                 = new TH2F("ESD_TrueClusBG_Pt_M02","ESD_TrueClusBG_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusBGPtvsM02[iCut]);
      fHistoTrueClusGammaPtvsM02[iCut]              = new TH2F("ESD_TrueClusGamma_Pt_M02","ESD_TrueClusGamma_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusGammaPtvsM02[iCut]);                  
      fHistoTrueClusGammaFromPi0PtvsM02[iCut]       = new TH2F("ESD_TrueClusGamma_FromPi0_Pt_M02","ESD_TrueClusGamma_FromPi0_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusGammaFromPi0PtvsM02[iCut]);                  
      fHistoTrueClusGammaFromEtaPtvsM02[iCut]       = new TH2F("ESD_TrueClusGamma_FromEta_Pt_M02","ESD_TrueClusGamma_FromEta_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusGammaFromEtaPtvsM02[iCut]);                  
      fHistoTrueClusElectronPtvsM02[iCut]           = new TH2F("ESD_TrueClusElectron_Pt_M02","ESD_TrueClusElectron_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusElectronPtvsM02[iCut]);                  
      fHistoTrueClusElectronFromPi0PtvsM02[iCut]    = new TH2F("ESD_TrueClusElectron_FromPi0_Pt_M02","ESD_TrueClusElectron_FromPi0_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusElectronFromPi0PtvsM02[iCut]);                  
      fHistoTrueClusElectronFromEtaPtvsM02[iCut]    = new TH2F("ESD_TrueClusElectron_FromEta_Pt_M02","ESD_TrueClusElectron_FromEta_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusElectronFromEtaPtvsM02[iCut]);                  
      fHistoTrueClusElectronFromGammaPtvsM02[iCut]  = new TH2F("ESD_TrueClusElectron_FromGamma_Pt_M02","ESD_TrueClusElectron_FromGamma_Pt_M02",500,0,50,500,0,5);
      fTrueList[iCut]->Add(fHistoTrueClusElectronFromGammaPtvsM02[iCut]);                  
      
      if (GetSelectedMesonID() < 2){
        fHistoTrueClusPrimPi0PtvsM02[iCut]                  = new TH2F("ESD_TrueClusFromPrimPi0_Pt_M02","ESD_TrueClusFromPrimPi0_Pt_M02",500,0,50,500,0,5);
        fTrueList[iCut]->Add(fHistoTrueClusPrimPi0PtvsM02[iCut]);
        fHistoTrueClusSecPi0PtvsM02[iCut]                   = new TH2F("ESD_TrueClusFromSecPi0_Pt_M02","ESD_TrueClusFromSecPi0_Pt_M02",500,0,50,500,0,5);
        fTrueList[iCut]->Add(fHistoTrueClusSecPi0PtvsM02[iCut]);
        fHistoTrueClusSecPi0FromK0sPtvsM02[iCut]            = new TH2F("ESD_TrueClusFromSecPi0FromK0s_Pt_M02","ESD_TrueClusFromSecPi0FromK0s_Pt_M02",500,0,50,500,0,5);
        fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromK0sPtvsM02[iCut]);
        fHistoTrueClusSecPi0FromLambdaPtvsM02[iCut]         = new TH2F("ESD_TrueClusFromSecPi0FromLambda_Pt_M02","ESD_TrueClusFromSecPi0FromLambda_Pt_M02",500,0,50,500,0,5);
        fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromLambdaPtvsM02[iCut]);

        fHistoTrueClusPartConvPrimPi0PtvsM02[iCut]          = new TH2F("ESD_TrueClusPartConvFromPrimPi0_Pt_M02","ESD_TrueClusPartConvFromPrimPi0_Pt_M02",500,0,50,500,0,5);
        fTrueList[iCut]->Add(fHistoTrueClusPartConvPrimPi0PtvsM02[iCut]);
        fHistoTrueClusPartConvSecPi0PtvsM02[iCut]           = new TH2F("ESD_TrueClusPartConvFromSecPi0_Pt_M02","ESD_TrueClusPartConvFromSecPi0_Pt_M02",500,0,50,500,0,5);
        fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0PtvsM02[iCut]);
        fHistoTrueClusPartConvSecPi0FromK0sPtvsM02[iCut]    = new TH2F("ESD_TrueClusPartConvFromSecPi0FromK0s_Pt_M02","ESD_TrueClusPartConvFromSecPi0FromK0s_Pt_M02",500,0,50,500,0,5);
        fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0FromK0sPtvsM02[iCut]);
        fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02[iCut] = new TH2F("ESD_TrueClusPartConvFromSecPi0FromLamdba_Pt_M02","ESD_TrueClusPartConvFromSecPi0FromLamdba_Pt_M02",500,0,50,500,0,5);
        fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02[iCut]);
      }
      
        
      fHistoTrueClusBGPtvsSource[iCut]                      = new TH2F("ESD_TrueClusBG_Pt_Source","ESD_TrueClusBG_Pt_Source",500, 0, 50,10, 0, 10);
      fTrueList[iCut]->Add(fHistoTrueClusBGPtvsSource[iCut]);                  
      fHistoTrueClusGammaPtvsSource[iCut]                   = new TH2F("ESD_TrueClusGamma_Pt_Source","ESD_TrueClusGamma_Pt_Source",500, 0, 50,8, 0, 8);
      fTrueList[iCut]->Add(fHistoTrueClusGammaPtvsSource[iCut]);                  
      fHistoTrueClusElectronPtvsSource[iCut]                = new TH2F("ESD_TrueClusElectron_Pt_Source","ESD_TrueClusElectron_Pt_Source",500, 0, 50,9, 0, 9);
      fTrueList[iCut]->Add(fHistoTrueClusElectronPtvsSource[iCut]);                  
      fHistoTrueMergedMissedPDG[iCut]                       = new TH1F("ESD_TrueMergedMissed_PDG","ESD_TrueMergedMissed_PDG",10000, -1.5, 9998.5);
      fTrueList[iCut]->Add(fHistoTrueMergedMissedPDG[iCut]);                  
      fHistoTrueMergedPartConvMissedPDG[iCut]               = new TH1F("ESD_TrueMergedPartConvMissed_PDG","ESD_TrueMergedPartConvMissed_PDG",10000, -1.5, 9998.5);
      fTrueList[iCut]->Add(fHistoTrueMergedPartConvMissedPDG[iCut]);                  

      if (fDoMesonQA > 1){
        fHistoTrueClusMergedInvMassvsPt[iCut]                 = new TH2F("ESD_TrueClusMerged_InvMass_Pt","ESD_TrueClusMerged_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
        fTrueList[iCut]->Add(fHistoTrueClusMergedInvMassvsPt[iCut]); 
        fHistoTrueClusPi0InvMassvsPt [iCut]                   = new TH2F("ESD_TrueClusFromPi0_InvMass_Pt","ESD_TrueClusFromPi0_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
        fTrueList[iCut]->Add(fHistoTrueClusPi0InvMassvsPt[iCut]);
        fHistoTrueClusEtaInvMassvsPt[iCut]                    = new TH2F("ESD_TrueClusFromEta_InvMass_Pt","ESD_TrueClusFromEta_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
        fTrueList[iCut]->Add(fHistoTrueClusEtaInvMassvsPt[iCut]);                  
        fHistoTrueClusMergedPartConvInvMassvsPt[iCut]         = new TH2F("ESD_TrueClusMergedPartConv_InvMass_Pt","ESD_TrueClusMergedPartConv_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
        fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvInvMassvsPt[iCut]);
        fHistoTrueClusPartConvPi0InvMassvsPt[iCut]            = new TH2F("ESD_TrueClusPartConvFromPi0_InvMass_Pt","ESD_TrueClusPartConvFromPi0_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
        fTrueList[iCut]->Add(fHistoTrueClusPartConvPi0InvMassvsPt[iCut]);
        fHistoTrueClusPartConvEtaInvMassvsPt[iCut]            = new TH2F("ESD_TrueClusPartConvFromEta_InvMass_Pt","ESD_TrueClusPartConvFromEta_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
        fTrueList[iCut]->Add(fHistoTrueClusPartConvEtaInvMassvsPt[iCut]);
        fHistoTrueClusBGInvMassvsPt[iCut]                     = new TH2F("ESD_TrueClusBG_InvMass_Pt","ESD_TrueClusBG_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
        fTrueList[iCut]->Add(fHistoTrueClusBGInvMassvsPt[iCut]);
        fHistoTrueClusGammaInvMassvsPt[iCut]                  = new TH2F("ESD_TrueClusGamma_InvMass_Pt","ESD_TrueClusGamma_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
        fTrueList[iCut]->Add(fHistoTrueClusGammaInvMassvsPt[iCut]);                  
        fHistoTrueClusElectronInvMassvsPt[iCut]               = new TH2F("ESD_TrueClusElectron_InvMass_Pt","ESD_TrueClusElectron_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
        fTrueList[iCut]->Add(fHistoTrueClusElectronInvMassvsPt[iCut]);                  
        fHistoTrueClusMergedPartConvELeadInvMassvsPt[iCut]          = new TH2F("ESD_TrueClusMergedPartConvLeadE_InvMass_Pt","ESD_TrueClusMergedPartConvLeadE_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
        fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvELeadInvMassvsPt[iCut]);
        if (GetSelectedMesonID() < 2){
          fHistoTrueClusPrimPi0InvMassvsPt[iCut]                    = new TH2F("ESD_TrueClusFromPrimPi0_InvMass_Pt","ESD_TrueClusFromPrimPi0_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
          fTrueList[iCut]->Add(fHistoTrueClusPrimPi0InvMassvsPt[iCut]);
          fHistoTrueClusSecPi0InvMassvsPt[iCut]                     = new TH2F("ESD_TrueClusFromSecPi0_InvMass_Pt","ESD_TrueClusFromSecPi0_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
          fTrueList[iCut]->Add(fHistoTrueClusSecPi0InvMassvsPt[iCut]);
          fHistoTrueClusSecPi0FromK0sInvMassvsPt[iCut]              = new TH2F("ESD_TrueClusFromSecPi0FromK0s_InvMass_Pt","ESD_TrueClusFromSecPi0FromK0s_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
          fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromK0sInvMassvsPt[iCut]);
          fHistoTrueClusSecPi0FromLambdaInvMassvsPt[iCut]           = new TH2F("ESD_TrueClusFromSecPi0FromLambda_InvMass_Pt","ESD_TrueClusFromSecPi0FromLambda_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
          fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromLambdaInvMassvsPt[iCut]);
          fHistoTrueClusPartConvPrimPi0InvMassvsPt[iCut]            = new TH2F("ESD_TrueClusPartConvFromPrimPi0_InvMass_Pt","ESD_TrueClusPartConvFromPrimPi0_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
          fTrueList[iCut]->Add(fHistoTrueClusPartConvPrimPi0InvMassvsPt[iCut]);
          fHistoTrueClusPartConvSecPi0InvMassvsPt[iCut]             = new TH2F("ESD_TrueClusPartConvFromSecPi0_InvMass_Pt","ESD_TrueClusPartConvFromSecPi0_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
          fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0InvMassvsPt[iCut]);
          fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt[iCut]      = new TH2F("ESD_TrueClusPartConvFromSecPi0FromK0s_InvMass_Pt","ESD_TrueClusPartConvFromSecPi0FromK0s_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
          fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt[iCut]);
          fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt[iCut]   = new TH2F("ESD_TrueClusPartConvFromSecPi0FromLamdba_InvMass_Pt","ESD_TrueClusPartConvFromSecPi0FromLamdba_InvMass_Pt",invMassBins, startMass, endMass,500, 0, 50);
          fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt[iCut]);
        }
      }
      
      if (fDoMesonQA > 0){
        if (GetSelectedMesonID() != 2){
          fHistoTruePi0PtY[iCut]                            = new TH2F("ESD_TruePi0_Pt_Y","ESD_TruePi0_Pt_Y",350,0.03,35.,150,-1.5,1.5);
          SetLogBinningXTH2(fHistoTruePi0PtY[iCut]);
          fTrueList[iCut]->Add(fHistoTruePi0PtY[iCut]);
          fHistoTruePi0PtAlpha[iCut]                        = new TH2F("ESD_TruePi0_Pt_Alpha","ESD_TruePi0_Pt_Alpha",350,0.03,35.,100,0,1);
          SetLogBinningXTH2(fHistoTruePi0PtAlpha[iCut]);
          fTrueList[iCut]->Add(fHistoTruePi0PtAlpha[iCut]);
          fHistoTruePi0PtOpenAngle[iCut]                    = new TH2F("ESD_TruePi0_Pt_OpenAngle","ESD_TruePi0_Pt_OpenAngle",350,0.03,35.,100,0,0.5);
          SetLogBinningXTH2(fHistoTruePi0PtOpenAngle[iCut]);
          fTrueList[iCut]->Add(fHistoTruePi0PtOpenAngle[iCut]);
        }
        if (GetSelectedMesonID() != 1){
          fHistoTrueEtaPtY[iCut]                            = new TH2F("ESD_TrueEta_Pt_Y","ESD_TrueEta_Pt_Y",350,0.03,35.,150,-1.5,1.5);
          SetLogBinningXTH2(fHistoTrueEtaPtY[iCut]);
          fTrueList[iCut]->Add(fHistoTrueEtaPtY[iCut]);                  
          fHistoTrueEtaPtAlpha[iCut]                        = new TH2F("ESD_TrueEta_Pt_Alpha","ESD_TrueEta_Pt_Alpha",350,0.03,35.,100,0,1);
          SetLogBinningXTH2(fHistoTrueEtaPtAlpha[iCut]);
          fTrueList[iCut]->Add(fHistoTrueEtaPtAlpha[iCut]);                  
          fHistoTrueEtaPtOpenAngle[iCut]                    = new TH2F("ESD_TrueEta_Pt_OpenAngle","ESD_TrueEta_Pt_OpenAngle",350,0.03,35.,180,0,1.8);
          SetLogBinningXTH2(fHistoTrueEtaPtOpenAngle[iCut]);
          fTrueList[iCut]->Add(fHistoTrueEtaPtOpenAngle[iCut]);
        }
      }
      
      if (fIsMC == 2){
        fHistoTrueClusMergedPtvsM02[iCut]->Sumw2();
        fHistoTrueClusPi0PtvsM02[iCut]->Sumw2();
        fHistoTrueClusEtaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusMergedPartConvPtvsM02[iCut]->Sumw2();
        fHistoTrueClusMergedPartConvELeadPtvsM02[iCut]->Sumw2();
        fHistoTrueClusPartConvPi0PtvsM02[iCut]->Sumw2();                  
        fHistoTrueClusPartConvEtaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusPartConvGammaPtvsM02[iCut]->Sumw2();
        fHistoTrueMergedMissedPDG[iCut]->Sumw2();
        fHistoTrueMergedPartConvMissedPDG[iCut]->Sumw2();
        fHistoTrueMergedPartConvNonLeadingPtvsM02[iCut]->Sumw2();
        fHistoTrueClusBGPtvsM02[iCut]->Sumw2();
        fHistoTrueClusGammaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusGammaFromPi0PtvsM02[iCut]->Sumw2();
        fHistoTrueClusGammaFromEtaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusElectronPtvsM02[iCut]->Sumw2();
        fHistoTrueClusElectronFromPi0PtvsM02[iCut]->Sumw2();
        fHistoTrueClusElectronFromEtaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusElectronFromGammaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusBGPtvsSource[iCut]->Sumw2();
        fHistoTrueClusGammaPtvsSource[iCut]->Sumw2();
        fHistoTrueClusElectronPtvsSource[iCut]->Sumw2();
        if (GetSelectedMesonID() < 2){
          fHistoTrueClusPrimPi0PtvsM02[iCut]->Sumw2();
          fHistoTrueClusSecPi0PtvsM02[iCut]->Sumw2();
          fHistoTrueClusSecPi0FromK0sPtvsM02[iCut]->Sumw2();
          fHistoTrueClusSecPi0FromLambdaPtvsM02[iCut]->Sumw2();
          fHistoTrueClusPartConvPrimPi0PtvsM02[iCut]->Sumw2();
          fHistoTrueClusPartConvSecPi0PtvsM02[iCut]->Sumw2();
          fHistoTrueClusPartConvSecPi0FromK0sPtvsM02[iCut]->Sumw2();
          fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02[iCut]->Sumw2();                    
        }
  
        if (fDoMesonQA > 1){
          fHistoTrueClusMergedInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusPi0InvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusEtaInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusMergedPartConvInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusPartConvPi0InvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusPartConvEtaInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusBGInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusGammaInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusElectronInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusMergedPartConvELeadInvMassvsPt[iCut]->Sumw2();
          if (GetSelectedMesonID() < 2){
            fHistoTrueClusPrimPi0InvMassvsPt[iCut]->Sumw2();  
            fHistoTrueClusSecPi0InvMassvsPt[iCut]->Sumw2();
            fHistoTrueClusSecPi0FromK0sInvMassvsPt[iCut]->Sumw2();
            fHistoTrueClusSecPi0FromLambdaInvMassvsPt[iCut]->Sumw2();
            fHistoTrueClusPartConvPrimPi0InvMassvsPt[iCut]->Sumw2();
            fHistoTrueClusPartConvSecPi0InvMassvsPt[iCut]->Sumw2();
            fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt[iCut]->Sumw2();
            fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt[iCut]->Sumw2();                    
          }
        }
        
        if (fDoMesonQA > 0 ){
          if (GetSelectedMesonID() != 2){
            fHistoTruePi0PtY[iCut]->Sumw2();
            fHistoTruePi0PtAlpha[iCut]->Sumw2();
            fHistoTruePi0PtOpenAngle[iCut]->Sumw2();
          }
          if (GetSelectedMesonID() != 1){
            fHistoTrueEtaPtY[iCut]->Sumw2();  
            fHistoTrueEtaPtAlpha[iCut]->Sumw2();
            fHistoTrueEtaPtOpenAngle[iCut]->Sumw2();
          }
        }
      }
    }
    
  }
  

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
  
  if(fV0Reader)
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
  if(fV0Reader)
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

      
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))) continue;
    if(((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetCutHistograms());
    }

    if(fSetPlotHistsExtQA){
      if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetExtQAHistograms()){
        fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetExtQAHistograms());
      }
    }
    if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
    if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
    }
  }
  PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloMerged::Notify()
{
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }  
    
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift()){
      fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      continue; // No Eta Shift requested, continue
    }
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
      fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      continue;
    }
    else{
      printf(" Gamma Conversion Task %s :: Eta Shift Manually Set to %f \n\n",
          (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
      fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
    }
  }
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::UserExec(Option_t *)
{
  //
  // Called for each event
  //

  if(fIsMC> 0) fMCEvent = MCEvent();
  if(fMCEvent == NULL) fIsMC = 0;

  fInputEvent = InputEvent();

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fHistoNEvents[iCut]->Fill(eventQuality);
      if (fIsMC==2) 
        fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
    }
    return;
  }
  
  if(fIsMC>0 && fInputEvent->IsA()==AliESDEvent::Class()){
    fMCStack = fMCEvent->Stack();
    if(fMCStack == NULL) fIsMC = 0;
  }
  
  // ------------------- BeginEvent ----------------------------
    
  for(Int_t iCut = 0; iCut<fnCuts; iCut++){  
    fiCut = iCut;
    
    Bool_t isRunningEMCALrelAna = kFALSE;
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;
    
    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);
    
    if(fIsMC==2){
      Float_t xsection = -1.; Float_t ntrials = -1.;
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials);
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}",ntrials);
    }

    fWeightJetJetMC = 1;
    //     cout << fMCEvent << endl;
    Bool_t isMCJet = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC );
    if (!isMCJet){
      fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
      if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(10);
      continue;
    }

    Bool_t triggered = kTRUE;
    if(eventNotAccepted){
    // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      fHistoNEvents[iCut]->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
      if (eventNotAccepted==3 && fIsMC>0){
        triggered = kFALSE;
      } else {  
        continue;
      }
    }

    if(eventQuality != 0){// Event Not Accepted
      //cout << "event rejected due to: " <<eventQuality << endl;
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC);
      if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);

      continue;
    }
    if (triggered == kTRUE) {
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
      if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here

      fHistoNGoodESDTracks[iCut]->Fill(  fV0Reader->GetNumberOfPrimaryTracks(), 
                        fWeightJetJetMC);
      fHistoVertexZ[iCut]->Fill(  fInputEvent->GetPrimaryVertex()->GetZ(), 
                    fWeightJetJetMC);
      fHistoSPDClusterTrackletBackground[iCut]->Fill(  fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),
                              (fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)), 
                              fWeightJetJetMC);
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)
        fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A(),
                      fWeightJetJetMC);
      else 
        fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C(),
                      fWeightJetJetMC);
    }
    if(fIsMC> 0){
      // Process MC Particle
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                          ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                          fMCEvent);
        }
        else if(fInputEvent->IsA()==AliAODEvent::Class()){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                          ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                          fInputEvent);
        }

        if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader()){
          for(Int_t i = 0;i<(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader())->GetEntries();i++){
            TString nameBin= fHistoMCHeaders[iCut]->GetXaxis()->GetBinLabel(i+1);
            if (nameBin.CompareTo("")== 0){
              TString nameHeader = ((TObjString*)((TList*)((AliConvEventCuts*)fEventCutArray->At(iCut))
                                ->GetAcceptedHeader())->At(i))->GetString();
              fHistoMCHeaders[iCut]->GetXaxis()->SetBinLabel(i+1,nameHeader.Data());
            }
          }
        }
      }
    }
    if(fIsMC>0) ProcessMCParticles();
    
    if (triggered==kFALSE) continue;
    
    // it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
    ProcessClusters();                      // process calo clusters

    fHistoNClusterCandidates[iCut]->Fill(  fNClusterCandidates, 
                        fWeightJetJetMC);
    fHistoNClusterMergedCandidates[iCut]->Fill(  fNClusterMergedCandidates, 
                        fWeightJetJetMC);

    fHistoNGoodESDTracksVsNClusterCandidates[iCut]->Fill(  fV0Reader->GetNumberOfPrimaryTracks(),
                                fNClusterCandidates,
                                fWeightJetJetMC);
//   
//                   fVectorDoubleCountTruePi0s.clear();
//                   fVectorDoubleCountTrueEtas.clear();

  }
  
  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessClusters(){
  
  Int_t nclus = 0;
  nclus = fInputEvent->GetNumberOfCaloClusters();
  
//   cout << nclus << endl;
  
  if(nclus == 0)  return;
  
  fNClusterCandidates = 0;
  fNClusterMergedCandidates = 0;
  
  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent);

  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  // Loop over EMCal clusters
  for(Long_t i = 0; i < nclus; i++){
  
    // select clusters with open cuts for normalizations histos
    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));

    if (!clus){
      delete clus;
      continue;
    }
    
    // if open cluster cuts are not fullfilled I can abort
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fIsMC)){
      delete clus;
      continue;
    }
    
    // TLorentzvector with cluster  
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);
      
    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());
  
    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
    if(!PhotonCandidate){
      delete clus;
      delete tmpvec;
      if (PhotonCandidate) delete PhotonCandidate;
      continue;
    }
    
    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton();
    PhotonCandidate->SetCaloClusterRef(i);
    
    // get MC label
    if(fIsMC> 0){
      Int_t* mclabelsCluster  = clus->GetLabels();
//       cout << "cluster: " << i << endl;
      Int_t nValidClusters    = 0;
      
      if (clus->GetNLabels()>0){
        for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
//           TParticle *dummy    = NULL;
          if (mclabelsCluster[k]>0){
//             dummy             = fMCStack->Particle(mclabelsCluster[k]);
//             if (dummy->R() < 407.0){
              if (nValidClusters< 50)PhotonCandidate->SetCaloPhotonMCLabel(nValidClusters,mclabelsCluster[k]);
              nValidClusters++;
//             }  
          }
        }
      }
      PhotonCandidate->SetNCaloPhotonMCLabels(nValidClusters);
      
    }
    
    fIsFromMBHeader = kTRUE; 
    fIsOverlappingWithOtherHeader = kFALSE;

    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
      if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCStack, fInputEvent) == 0) fIsFromMBHeader = kFALSE;
      if (clus->GetNLabels()>1){
        Int_t* mclabelsCluster = clus->GetLabels();
        for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
          if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCStack, fInputEvent) == 0) fIsOverlappingWithOtherHeader = kTRUE;
        }
      }
    }
    
    // check whether photon is from correct header 
    if (fIsFromMBHeader && fIsOverlappingWithOtherHeader){
      fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC);
    }  
    // only cluster candidates from correct header will be processed fully  
    if (fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
      fHistoClusGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC);  
      if (fDoClusterQA > 0){
        fHistoClusNCellsPt[fiCut]->Fill(clus->GetNCells(), PhotonCandidate->Pt(), fWeightJetJetMC);
      }  
      fNClusterCandidates++;
    } else {
      delete clus;
      delete tmpvec;
      delete PhotonCandidate;
      continue;
    }
    
    // check whether photon fullfill merged cluster cuts as well
    if(!((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fIsMC)){
      delete clus;
      delete tmpvec;
      delete PhotonCandidate;
      continue;
    }else {
      fNClusterMergedCandidates++;   
    }

    AliAODCaloCluster* clusSub1 = new AliAODCaloCluster();
    AliAODCaloCluster* clusSub2 = new AliAODCaloCluster();
    
    
    // split clusters according to their shares in the cluster (NLM == 1) needs to be treated differently
    if (((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->GetMinNLMCut() == 1 && 
      ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->GetMaxNLMCut() == 1 ){
      
      Int_t absCellIdFirst    = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus, fInputEvent);
      Int_t absCellIdSecond   = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindSecondLargestCellInCluster(clus, fInputEvent);
      
      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->SplitEnergy(absCellIdFirst, absCellIdSecond,
                                      clus, fInputEvent, fIsMC, clusSub1, clusSub2);
      
      
    } else if (((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->GetMinNLMCut() > 1 ){
      const Int_t   nc = clus->GetNCells();    
      Int_t   absCellIdList[nc]; 
      Float_t maxEList[nc];     
      Int_t nlm = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetNumberOfLocalMaxima(clus, fInputEvent, absCellIdList, maxEList);
      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->SplitEnergy(absCellIdList[0], absCellIdList[1],
                                    clus, fInputEvent, fIsMC, clusSub1, clusSub2);
    }
    fHistoClusMergedPtvsM02[fiCut]->Fill( PhotonCandidate->Pt(),
                          clus->GetM02(),
                          fWeightJetJetMC);  
            
    // TLorentzvector with sub cluster 1
    TLorentzVector clusterVector1;
    clusSub1->GetMomentum(clusterVector1,vertex);
    TLorentzVector* tmpvec1 = new TLorentzVector();
    tmpvec1->SetPxPyPzE(clusterVector1.Px(),clusterVector1.Py(),clusterVector1.Pz(),clusterVector1.E());
    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate1=new AliAODConversionPhoton(tmpvec1);
    if(!PhotonCandidate1){
      delete clusSub1;
      delete tmpvec1;
      if (PhotonCandidate1) delete PhotonCandidate1;
      continue;
    }
    // Flag Photon as CaloPhoton
    PhotonCandidate1->SetIsCaloPhoton();
    // get MC label
//     if(fIsMC> 0){
//       Int_t* mclabelsCluster = clusSub1->GetLabels();
//       PhotonCandidate1->SetNCaloPhotonMCLabels(clusSub1->GetNLabels());
// 
//       if (clusSub1->GetNLabels()>0){
//         for (Int_t k =0; k< (Int_t)clusSub1->GetNLabels(); k++){
//           if (k< 50)PhotonCandidate1->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
//         }
//       }
//     }
    
    // TLorentzvector with sub cluster 2
    TLorentzVector clusterVector2;
    clusSub2->GetMomentum(clusterVector2,vertex);
    TLorentzVector* tmpvec2 = new TLorentzVector();
    tmpvec2->SetPxPyPzE(clusterVector2.Px(),clusterVector2.Py(),clusterVector2.Pz(),clusterVector2.E());
    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate2=new AliAODConversionPhoton(tmpvec2);
    if(!PhotonCandidate2){
      delete clusSub2;
      delete tmpvec2;
      if (PhotonCandidate2) delete PhotonCandidate2;
      continue;
    }
    // Flag Photon as CaloPhoton
    PhotonCandidate2->SetIsCaloPhoton();
    // get MC label
//     if(fIsMC> 0){
//       Int_t* mclabelsCluster = clusSub2->GetLabels();
//       PhotonCandidate2->SetNCaloPhotonMCLabels(clusSub2->GetNLabels());
// 
//       if (clusSub2->GetNLabels()>0){
//         for (Int_t k =0; k< (Int_t)clusSub2->GetNLabels(); k++){
//           if (k< 50)PhotonCandidate2->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
//         }
//       }
//     }
    
    AliAODConversionMother *pi0cand = new AliAODConversionMother(PhotonCandidate1,PhotonCandidate2);
            
    if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
      fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), fWeightJetJetMC);
      fHistoClusMergedPtvsM02Accepted[fiCut]->Fill( PhotonCandidate->Pt(), clus->GetM02(), fWeightJetJetMC);  
      if (fDoClusterQA > 0){
        fHistoClusMergedNCellsPt[fiCut]->Fill(clus->GetNCells(), PhotonCandidate->Pt(), fWeightJetJetMC);
        Double_t energyAround   = 0;
        Double_t nCellsAround   = 0;
        for (Int_t j = 0; j < nclus; j++){
          if (j == i) continue;
          AliVCluster* clusTemp   = NULL;
          if(fInputEvent->IsA()==AliESDEvent::Class()) clusTemp = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(j));
          else if(fInputEvent->IsA()==AliAODEvent::Class()) clusTemp = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(j));

          if (!clusTemp){
            delete clusTemp;
            continue;
          }
          
          Double_t R = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDistanceBetweenClusters(clus, clusTemp);
          
          if (R < 0.15){
            nCellsAround        = nCellsAround+clusTemp->GetNCells();
            energyAround        = energyAround+clusTemp->E();
          }
          delete clusTemp;
        }  
        fHistoClusMergedNCellsAroundPt[fiCut]->Fill(nCellsAround, PhotonCandidate->Pt(), fWeightJetJetMC);
        fHistoClusMergedNCellsAroundAndInPt[fiCut]->Fill(nCellsAround+clus->GetNCells(), PhotonCandidate->Pt(), fWeightJetJetMC);
        fHistoClusMergedEAroundE[fiCut]->Fill(energyAround, clus->E(), fWeightJetJetMC);
        if (fIsMC > 0){
          fHistoClusMergedNParticlePt[fiCut]->Fill(clus->GetNLabels(), PhotonCandidate->Pt(), fWeightJetJetMC);
        }  
      }  
      if (fDoMesonQA > 0 ){
        fHistoMotherPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
        fHistoMotherPtAlpha[fiCut]->Fill(pi0cand->Pt(),fabs(pi0cand->GetAlpha()), fWeightJetJetMC);
        fHistoMotherPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(), fWeightJetJetMC);
      }
    }
    
    
    if(fIsMC> 0 && PhotonCandidate && PhotonCandidate1 && PhotonCandidate2){
      ProcessTrueClusterCandidates(PhotonCandidate, clus->GetM02(), PhotonCandidate1, PhotonCandidate2);
    }

    delete pi0cand;                  
    delete clusSub1;
    delete tmpvec1;
    delete PhotonCandidate1;
    delete clusSub2;                  
    delete tmpvec2;    
    delete PhotonCandidate2;    
      
    delete clus;
    delete tmpvec;
    delete PhotonCandidate;
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessTrueClusterCandidates(AliAODConversionPhoton *TrueClusterCandidate, Float_t m02,
                                    AliAODConversionPhoton *TrueSubClusterCandidate1,
                                    AliAODConversionPhoton *TrueSubClusterCandidate2)
{

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  TParticle *Photon = NULL;
  if (!TrueClusterCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
  
  if (TrueClusterCandidate->GetNCaloPhotonMCLabels()>0) Photon = fMCStack->Particle(TrueClusterCandidate->GetCaloPhotonMCLabel(0));
    else return;
    
  if(Photon == NULL){
  //    cout << "no photon" << endl;
    return;
  }

  AliAODConversionMother *mesoncand = new AliAODConversionMother(TrueSubClusterCandidate1,TrueSubClusterCandidate2);
  Bool_t mesonIsSelected            = (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesoncand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()));
  if (!mesonIsSelected){
    delete mesoncand;
    return;
  }
  Int_t pdgCodeParticle             = Photon->GetPdgCode();
  TrueClusterCandidate->SetCaloPhotonMCFlags(fMCStack);
  if (TrueClusterCandidate->GetNCaloPhotonMCLabels()>1 && fEnableDetailedPrintOut){
    cout << endl << endl << "Cluster energy: " << TrueClusterCandidate->E() << endl;;
    TrueClusterCandidate->PrintCaloMCLabelsAndInfo(fMCStack);
    TrueClusterCandidate->PrintCaloMCFlags();
  }
  
  Bool_t filled                     = kFALSE;
  
  if(fIsFromMBHeader && !fIsOverlappingWithOtherHeader){  
    Int_t clusterClass    = 0;
    Bool_t isPrimary      = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, TrueClusterCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    
    // cluster classification:
    // 1    - nice merged cluster (2 gamma | contributions from 2 gamma) from pi0/eta
    // 2    - contribution from only 1 partner (1 gamma, 1 fully coverted gamma) from pi0/eta
    // 3    - contribution from part of 1 partner (1 electron) from pi0/eta
    if (TrueClusterCandidate->IsMerged() || TrueClusterCandidate->IsMergedPartConv()){
        clusterClass    = 1;
    } else if (TrueClusterCandidate->GetNCaloPhotonMotherMCLabels()> 0){
//       cout << TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0) << endl;
      if (TrueClusterCandidate->IsLargestComponentElectron() || TrueClusterCandidate->IsLargestComponentPhoton()){ 
        if (fMCStack->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0))->GetPdgCode() == 111 || fMCStack->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0))->GetPdgCode() == 221 ){
          if ( TrueClusterCandidate->IsConversion() && !TrueClusterCandidate->IsConversionFullyContained() ){
            clusterClass  = 3;
          } else {
            clusterClass  = 2;  
          }
        }  
      } else if (TrueClusterCandidate->IsSubLeadingEM()){
        if (TrueClusterCandidate->GetNCaloPhotonMotherMCLabels()> 1){
          if (fEnableDetailedPrintOut) cout << "Is Subleading EM: "<<  TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1) << endl;
          if ( TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1) > -1){
            if (fMCStack->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1))->GetPdgCode() == 111 || fMCStack->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1))->GetPdgCode() == 221 ){
              clusterClass  = 2;  
            }
          }
        }  
      }  
    }  
    if (fEnableDetailedPrintOut) cout << "cluster class: " << clusterClass << endl;
  
    if (TrueClusterCandidate->IsMerged() && !TrueClusterCandidate->IsMergedPartConv()){
      if (fEnableDetailedPrintOut) cout << "merged" << endl;
      filled = kTRUE;
      fHistoTrueClusMergedPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
      if (fDoMesonQA > 1)fHistoTrueClusMergedInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
      Int_t motherLab = Photon->GetMother(0);
      if (motherLab > -1){ 
        if (fMCStack->Particle(motherLab)->GetPdgCode() == 111){
          fHistoTrueClusPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);                  
          if (fDoMesonQA > 1)fHistoTrueClusPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
          if (fDoMesonQA > 0 && GetSelectedMesonID() != 2){
            fHistoTruePi0PtY[fiCut]->Fill(mesoncand->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
            fHistoTruePi0PtAlpha[fiCut]->Fill(mesoncand->Pt(),fabs(mesoncand->GetAlpha()), fWeightJetJetMC);
            fHistoTruePi0PtOpenAngle[fiCut]->Fill(mesoncand->Pt(),mesoncand->GetOpeningAngle(), fWeightJetJetMC);
          }
        
          if (GetSelectedMesonID() < 2) {
            if (isPrimary) {
              fHistoTrueClusPrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
              if (fDoMesonQA > 1)fHistoTrueClusPrimPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
            } else {
              fHistoTrueClusSecPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
              if (fDoMesonQA > 1)fHistoTrueClusSecPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
              Int_t grandMaLab = fMCStack->Particle(motherLab)->GetMother(0);
              if (grandMaLab > -1){
                if (abs(fMCStack->Particle(grandMaLab)->GetPdgCode()) == 310){
                  fHistoTrueClusSecPi0FromK0sPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
                  if (fDoMesonQA > 1)fHistoTrueClusSecPi0FromK0sInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
                  
                } else if (abs(fMCStack->Particle(grandMaLab)->GetPdgCode()) == 3122){
                  fHistoTrueClusSecPi0FromLambdaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
                  if (fDoMesonQA > 1)fHistoTrueClusSecPi0FromLambdaInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
                  
                } 
              }
            }
          }
        } else if (fMCStack->Particle(motherLab)->GetPdgCode() == 221){
          fHistoTrueClusEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);                  
          if (fDoMesonQA > 1)fHistoTrueClusEtaInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
          if ( fDoMesonQA > 0 && GetSelectedMesonID() != 1 ){
            fHistoTrueEtaPtY[fiCut]->Fill(mesoncand->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
            fHistoTrueEtaPtAlpha[fiCut]->Fill(mesoncand->Pt(),fabs(mesoncand->GetAlpha()), fWeightJetJetMC);
            fHistoTrueEtaPtOpenAngle[fiCut]->Fill(mesoncand->Pt(),mesoncand->GetOpeningAngle(), fWeightJetJetMC);                      
          }
        } else {
          fHistoTrueMergedMissedPDG[fiCut]->Fill(fMCStack->Particle(motherLab)->GetPdgCode(), fWeightJetJetMC);
//             cout << "offending particle pdg: " << fMCStack->Particle(motherLab)->GetPdgCode() << endl;
        } 
      } else {
        fHistoTrueMergedMissedPDG[fiCut]->Fill(-1, fWeightJetJetMC);
//         cout << "no mother " << endl;
      }  
    }
    if (TrueClusterCandidate->IsMergedPartConv() && !filled){
      if (fEnableDetailedPrintOut) cout << "merged part conv" << endl;
      filled = kTRUE;
      fHistoTrueClusMergedPartConvPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
      if (fDoMesonQA > 1)fHistoTrueClusMergedPartConvInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
      
      if (TrueClusterCandidate->IsLargestComponentElectron()){
        fHistoTrueClusMergedPartConvELeadPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
        if (fDoMesonQA > 1) fHistoTrueClusMergedPartConvELeadInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
      }
      Int_t motherLab1 = Photon->GetMother(0);
      Int_t motherLab = Photon->GetMother(0);
      if (TrueClusterCandidate->IsLargestComponentElectron()){
        if (motherLab > -1){
          if (fMCStack->Particle(motherLab)->GetPdgCode() == 22)
            motherLab = fMCStack->Particle(Photon->GetMother(0))->GetMother(0);
        }
      }
      if (motherLab > -1){ 
        if (fMCStack->Particle(motherLab)->GetPdgCode() == 111){
          fHistoTrueClusPartConvPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
          if (fDoMesonQA > 1) fHistoTrueClusPartConvPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
          if ( fDoMesonQA > 0 && GetSelectedMesonID() != 2){
            fHistoTruePi0PtY[fiCut]->Fill(mesoncand->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
            fHistoTruePi0PtAlpha[fiCut]->Fill(mesoncand->Pt(),fabs(mesoncand->GetAlpha()), fWeightJetJetMC);
            fHistoTruePi0PtOpenAngle[fiCut]->Fill(mesoncand->Pt(),mesoncand->GetOpeningAngle(), fWeightJetJetMC);                      
          }
          if (GetSelectedMesonID() < 2) {
            if (isPrimary) {
              fHistoTrueClusPartConvPrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
              if (fDoMesonQA > 1) fHistoTrueClusPartConvPrimPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
            } else {
              fHistoTrueClusPartConvSecPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
              if (fDoMesonQA > 1) fHistoTrueClusPartConvSecPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
              Int_t grandMaLab = fMCStack->Particle(motherLab)->GetMother(0);
              if (grandMaLab > -1){
                if (abs(fMCStack->Particle(grandMaLab)->GetPdgCode()) == 310){
                  fHistoTrueClusPartConvSecPi0FromK0sPtvsM02[fiCut]->Fill( TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
                  if (fDoMesonQA > 1) fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
                } else if (abs(fMCStack->Particle(grandMaLab)->GetPdgCode()) == 3122){
                  fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
                  if (fDoMesonQA > 1) fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
                } 
              }
            }
          }  
        } else if (fMCStack->Particle(motherLab)->GetPdgCode() == 221){
          fHistoTrueClusPartConvEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
          if (fDoMesonQA > 1) fHistoTrueClusPartConvEtaInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
          if ( fDoMesonQA > 0 && GetSelectedMesonID() != 1 ){
            fHistoTrueEtaPtY[fiCut]->Fill(mesoncand->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
            fHistoTrueEtaPtAlpha[fiCut]->Fill(mesoncand->Pt(),fabs(mesoncand->GetAlpha()), fWeightJetJetMC);
            fHistoTrueEtaPtOpenAngle[fiCut]->Fill(mesoncand->Pt(),mesoncand->GetOpeningAngle(), fWeightJetJetMC);                      
          }
        } else if (TrueClusterCandidate->IsLargestComponentElectron() || TrueClusterCandidate->IsLargestComponentPhoton()) {
          fHistoTrueMergedPartConvMissedPDG[fiCut]->Fill(fMCStack->Particle(motherLab)->GetPdgCode(), fWeightJetJetMC);
//             cout << "offending particle pdg: " << fMCStack->Particle(motherLab)->GetPdgCode() << endl;
        } else {
          fHistoTrueMergedPartConvMissedPDG[fiCut]->Fill(fMCStack->Particle(motherLab)->GetPdgCode(), fWeightJetJetMC);
          fHistoTrueMergedPartConvNonLeadingPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
        }  
      } else if (motherLab1 > -1){
          if (fMCStack->Particle(motherLab1)->GetPdgCode() == 22){
            fHistoTrueClusPartConvGammaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
          } else {
            fHistoTrueMergedPartConvMissedPDG[fiCut]->Fill(fMCStack->Particle(motherLab1)->GetPdgCode(), fWeightJetJetMC);
          }    
      } else {
        fHistoTrueMergedPartConvMissedPDG[fiCut]->Fill(-1, fWeightJetJetMC);
      } 
    }
    
    if (!(TrueClusterCandidate->IsLargestComponentElectron() || TrueClusterCandidate->IsLargestComponentPhoton())){
      if (fEnableDetailedPrintOut) cout << "BG" << endl;
      fHistoTrueClusBGPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
      if (fDoMesonQA > 1) fHistoTrueClusBGInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
      if (abs(pdgCodeParticle) == 211) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 0.5, fWeightJetJetMC); // pi+/-
      else if (abs(pdgCodeParticle) == 2212) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 1.5, fWeightJetJetMC); // p
      else if (abs(pdgCodeParticle) == 321) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 2.5, fWeightJetJetMC); // K+-
      else if (abs(pdgCodeParticle) == 2112) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 3.5, fWeightJetJetMC); // n
      else if (abs(pdgCodeParticle) == 310) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 4.5, fWeightJetJetMC); // K0s
      else if (abs(pdgCodeParticle) == 3122) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 5.5, fWeightJetJetMC); // Lambda
      else if (abs(pdgCodeParticle) == 13) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 6.5, fWeightJetJetMC); // mu+/-
      else if (abs(pdgCodeParticle) == 130) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 7.5, fWeightJetJetMC); // K0l
      else fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 8.5, fWeightJetJetMC); // Rest
    }
      
    if ((TrueClusterCandidate->IsLargestComponentPhoton() || TrueClusterCandidate->IsConversionFullyContained() ) && !filled){
      if (fEnableDetailedPrintOut) cout << "photon" << endl;
      fHistoTrueClusGammaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
      if (fDoMesonQA > 1) fHistoTrueClusGammaInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
      Int_t motherLab = Photon->GetMother(0);
      if (motherLab == -1){
        fHistoTrueClusGammaPtvsSource[fiCut]->Fill(mesoncand->Pt(), 0.5, fWeightJetJetMC); // direct photon
      } else {
        Int_t pdgCodeMother = fMCStack->Particle(motherLab)->GetPdgCode();
        if (abs(pdgCodeMother) == 111)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(mesoncand->Pt(), 1.5, fWeightJetJetMC); // pi0
        else if (abs(pdgCodeMother) == 221)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(mesoncand->Pt(), 2.5, fWeightJetJetMC); // eta
        else if (abs(pdgCodeMother) == 331)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(mesoncand->Pt(), 3.5, fWeightJetJetMC); // eta'
        else if (abs(pdgCodeMother) == 223)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(mesoncand->Pt(), 4.5, fWeightJetJetMC); // omega
        else if (abs(pdgCodeMother) == 333)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(mesoncand->Pt(), 5.5, fWeightJetJetMC); // phi
        else if (abs(pdgCodeMother) == 3122)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(mesoncand->Pt(), 6.5, fWeightJetJetMC); // Lambda
        else 
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(mesoncand->Pt(), 7.5, fWeightJetJetMC); // rest
      }        
    }
    if (TrueClusterCandidate->IsLargestComponentElectron() && !filled){
      if (fEnableDetailedPrintOut) cout << "electron" << endl;
      fHistoTrueClusElectronPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
      if (fDoMesonQA > 1) fHistoTrueClusElectronInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
      Int_t motherLab = Photon->GetMother(0);
      if (motherLab == -1){
        fHistoTrueClusElectronPtvsSource[fiCut]->Fill(mesoncand->Pt(), 0.5, fWeightJetJetMC); // direct photon
      } else {
        Int_t pdgCodeMother = fMCStack->Particle(motherLab)->GetPdgCode();
        if (abs(pdgCodeMother) == 22)
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(mesoncand->Pt(), 1.5, fWeightJetJetMC); // gamma
        else if (abs(pdgCodeMother) == 111)
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(mesoncand->Pt(), 2.5, fWeightJetJetMC); // pi0
        else if (abs(pdgCodeMother) == 221)
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(mesoncand->Pt(), 3.5, fWeightJetJetMC); // eta
        else if (abs(pdgCodeMother) == 411 || abs(pdgCodeMother) == 421 || abs(pdgCodeMother) == 431 ||
                 abs(pdgCodeMother) == 443 || abs(pdgCodeMother) == 100443 || abs(pdgCodeMother) == 200443 ) 
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(mesoncand->Pt(), 4.5, fWeightJetJetMC); // c
        else if (abs(pdgCodeMother) == 511 || abs(pdgCodeMother) == 521 || abs(pdgCodeMother) == 531 ||
                 abs(pdgCodeMother) == 553 || abs(pdgCodeMother) == 100553 || abs(pdgCodeMother) == 200553 ) 
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(mesoncand->Pt(), 5.5, fWeightJetJetMC); // b
        else if (abs(pdgCodeMother) == 23 || abs(pdgCodeMother) == 24)
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(mesoncand->Pt(), 6.5, fWeightJetJetMC); // W/Z
        else if (abs(pdgCodeMother) == 15) 
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(mesoncand->Pt(), 7.5, fWeightJetJetMC); // tau
        else 
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(mesoncand->Pt(), 8.5, fWeightJetJetMC); // rest
      }  
    }
  }
  delete mesoncand;
  return;
}



//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
  
  // Loop over all primary MC particle
  for(UInt_t i = 0; i < fMCStack->GetNtrack(); i++) {
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){ 

      TParticle* particle = (TParticle *)fMCStack->Particle(i);
      if (!particle) continue;
      
      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }
      
      // check if particle is pi0/eta from di-photon decay
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
        ->MesonIsSelectedMC(particle,fMCStack,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        TParticle* daughter0 = (TParticle*)fMCStack->Particle(particle->GetFirstDaughter());
        TParticle* daughter1 = (TParticle*)fMCStack->Particle(particle->GetLastDaughter());
        
        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCStack, fInputEvent);
          }
        }
    
        if(particle->GetPdgCode() == 111 && GetSelectedMesonID() != 2){
          fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // All MC Pi0
          fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
          if (fIsMC==2)fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0 ){  
            if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
          }
        } else if(particle->GetPdgCode() == 221 && GetSelectedMesonID() != 1){
          fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // All MC Eta
          fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
          if (fIsMC==2)fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0 ){  
            if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
          }
        }
      
        // Check the acceptance for both gammas & whether they are counted as primaries as well
        Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, particle->GetFirstDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, particle->GetLastDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        
        if( kDaughter0IsPrim && kDaughter1IsPrim &&
          ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCStack) ||
          ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCStack) ){                    
          if(particle->GetPdgCode() == 111 && GetSelectedMesonID() != 2){
            fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Pi0 with gamma in acc
          } else if(particle->GetPdgCode() == 221 && GetSelectedMesonID() != 1){
            fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Eta with gamma in acc
          }
        }
      }
      Int_t gammaLabel    = -1;
      Int_t electronLabel = -1;
      Int_t positronLabel = -1;
      // check if particle is pi0/eta from Dalitz decay
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
        ->MesonIsSelectedMCDalitz(particle,fMCStack, electronLabel, positronLabel, gammaLabel, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        TParticle* gamma    = (TParticle*)fMCStack->Particle(gammaLabel);
        TParticle* electron = (TParticle*)fMCStack->Particle(electronLabel);
        TParticle* positron = (TParticle*)fMCStack->Particle(positronLabel);

        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCStack, fInputEvent);
          }
        }

        if(particle->GetPdgCode() == 111 && GetSelectedMesonID() != 2){
          fHistoMCPi0DalitzPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // All MC Pi0
          fHistoMCPi0DalitzWOWeightPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
          if (fIsMC==2)fHistoMCPi0DalitzWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0 ){  
            if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
          }
        } else if(particle->GetPdgCode() == 221 && GetSelectedMesonID() != 1){
          fHistoMCEtaDalitzPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // All MC Eta
          fHistoMCEtaDalitzWOWeightPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
          if (fIsMC==2)fHistoMCEtaDalitzWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0 ){  
            if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
          }
        }

        // Check the acceptance for both gammas & whether they are counted as primaries as well
        Bool_t kGammaIsPrim     = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, gammaLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kElectronIsPrim  = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, electronLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kPositronIsPrim  = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, positronLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        if( kGammaIsPrim && kElectronIsPrim && kPositronIsPrim &&
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma,fMCStack) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedElecMC(electron,fMCStack) || 
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedElecMC(positron,fMCStack)
          ){                    
          if(particle->GetPdgCode() == 111 && GetSelectedMesonID() != 2){
            fHistoMCPi0DalitzInAccPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Pi0 with gamma in acc
          } else if(particle->GetPdgCode() == 221 && GetSelectedMesonID() != 1){
            fHistoMCEtaDalitzInAccPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Eta with gamma in acc
          }
        }
      }
    }
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::SetLogBinningXTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloMerged::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
{
  if(tobechecked > -1)
  {
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else{
      vec.push_back(tobechecked);
      return false;
    }
  }
  return false;
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}


//_________________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked){
  if( ma.find(tobechecked) != ma.end() ) ma[tobechecked] += 1;
  else ma[tobechecked] = 2;
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist){
  map<Int_t, Int_t>::iterator it;
  for (it = ma.begin(); it != ma.end(); it++){
    hist->Fill(it->second, fWeightJetJetMC);
  }
  ma.clear();
  return;
}
