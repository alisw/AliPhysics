/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Florian Jonas                                                 *
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

#include "AliAnalysisTaskElectronStudies.h"
#include "TChain.h"
#include "TRandom.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "TFile.h"
#include "AliESDtrackCuts.h"
#include "AliAODMCParticle.h"
#include "AliAODConversionPhoton.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODCaloCluster.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliEMCALRecoUtilsBase.h"
#include "AliAODConversionMother.h"
#include "TObjectTable.h"

ClassImp(AliAnalysisTaskElectronStudies)
//________________________________________________________________________
AliAnalysisTaskElectronStudies::AliAnalysisTaskElectronStudies() : AliAnalysisTaskSE(),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fWeightJetJetMC(1),
  fOutputList(NULL),
  fAnalysisTree(NULL),
  fIsMC(0),
  fIsHeavyIon(0),
  fV0Reader(NULL),
  fV0ReaderName(""),
  fPIDResponse(NULL),
  fReaderGammas(NULL),
  fAODMCTrackArray(NULL),
  fGeomEMCAL(NULL),
  fCorrTaskSetting(""),
  fEventCuts(NULL),
  fClusterCutsEMC(NULL),
  fTMCuts(NULL),
  fConvCuts(NULL),
  fCaloUtils(NULL),
  fMinClsTPC(0),
  fMinFracClsTPC(0),
  fChi2PerClsTPC(9999),
  fMinClsITS(0),
  fEtaCut(9999),
  fPtCut(0),
  fYMCCut(9999),
  fMinNsigmaElec(-1),
  fMaxNsigmaElec(3),
  fMaxDCAxy(9999),
  fMaxDCAz(9999),
  fMatchingParamsPhi(),
  fMatchingParamsEta(),
  fUseRTrackMatching(kFALSE),
  fRTrackMatching(999),

  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fPtElectronTrack(NULL),
  fPtElectronTrackInEmcalAcc(NULL),
  fTruePtCluster(NULL),
  fTruePtElectronCluster(NULL),
  fTruePtElectronClusterMatchedWithTrack(NULL),
  fTruePtElectronTrack(NULL),
  fTruePtElectronTrackInEmcalAcc(NULL),
  fGenPtElectrons(NULL),
  fGenPtElectronsInEmcalAcc(NULL),
  fTreeBuffSize(60*1024*1024),
  fMemCountAOD(0),
  fTrackMatcherRunningMode(0),
  fBuffer_ClusterE(0), 
  fBuffer_ClusterM02(0), 
  fBuffer_ClusterM20(0), 
  fBuffer_Track_Pt(0), 
  fBuffer_Track_P(0), 
  fBuffer_Track_dEta(0), 
  fBuffer_Track_dPhi(0), 
  fBuffer_Track_NSigmaElec(0), 
  fBuffer_Track_IsFromV0(0), 
  fBuffer_MC_True_Cluster_E(0), 
  fBuffer_MC_True_Track_E(0), 
  fBuffer_MC_True_Track_Pt(0), 
  fBuffer_MC_True_Track_P(0), 
  fBuffer_MC_Track_Is_Electron(kFALSE), 
  fBuffer_MC_Cluster_Is_Electron(kFALSE), 
  fBuffer_MC_ClusterTrack_Same_Electron(kFALSE), 
  fBuffer_MC_JetJetWeight(1),
  fBuffer_MatchType(0),
  fTrackMatcher(NULL), 
  fTrackMatcherName("")
{
  SetEtaMatching(0.010,4.07,-2.5);
  SetPhiMatching(0.015,3.65,3.65);
}

AliAnalysisTaskElectronStudies::AliAnalysisTaskElectronStudies(const char *name) : AliAnalysisTaskSE(name),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fWeightJetJetMC(1),
  fOutputList(NULL),
  fAnalysisTree(NULL),
  fIsMC(0),
  fIsHeavyIon(0),
  fV0Reader(NULL),
  fV0ReaderName(""),
  fPIDResponse(NULL),
  fReaderGammas(NULL),
  fAODMCTrackArray(NULL),
  fGeomEMCAL(NULL),
  fCorrTaskSetting(""),
  fEventCuts(NULL),
  fClusterCutsEMC(NULL),
  fTMCuts(NULL),
  fConvCuts(NULL),
  fCaloUtils(NULL),
  fMinClsTPC(0),
  fMinFracClsTPC(0),
  fChi2PerClsTPC(9999),
  fMinClsITS(0),
  fEtaCut(9999),
  fPtCut(0),
  fYMCCut(9999),
  fMinNsigmaElec(-1),
  fMaxNsigmaElec(3),
  fMaxDCAxy(9999),
  fMaxDCAz(9999),
  fMatchingParamsPhi(),
  fMatchingParamsEta(),
  fUseRTrackMatching(kFALSE),
  fRTrackMatching(999),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fPtElectronTrack(NULL),
  fPtElectronTrackInEmcalAcc(NULL),
  fTruePtCluster(NULL),
  fTruePtElectronCluster(NULL),
  fTruePtElectronClusterMatchedWithTrack(NULL),
  fTruePtElectronTrack(NULL),
  fTruePtElectronTrackInEmcalAcc(NULL),
  fGenPtElectrons(NULL),
  fGenPtElectronsInEmcalAcc(NULL),
  fTreeBuffSize(60*1024*1024),
  fMemCountAOD(0),
  fTrackMatcherRunningMode(0),
  fBuffer_ClusterE(0), 
  fBuffer_ClusterM02(0), 
  fBuffer_ClusterM20(0), 
  fBuffer_Track_Pt(0), 
  fBuffer_Track_P(0), 
  fBuffer_Track_dEta(0), 
  fBuffer_Track_dPhi(0), 
  fBuffer_Track_NSigmaElec(0), 
  fBuffer_Track_IsFromV0(0), 
  fBuffer_MC_True_Cluster_E(0), 
  fBuffer_MC_True_Track_E(0), 
  fBuffer_MC_True_Track_Pt(0), 
  fBuffer_MC_True_Track_P(0), 
  fBuffer_MC_Track_Is_Electron(kFALSE), 
  fBuffer_MC_Cluster_Is_Electron(kFALSE), 
  fBuffer_MC_ClusterTrack_Same_Electron(kFALSE), 
  fBuffer_MC_JetJetWeight(1.),
  fBuffer_MatchType(0),
  fTrackMatcher(NULL), 
  fTrackMatcherName("") 
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  SetEtaMatching(0.010,4.07,-2.5);
  SetPhiMatching(0.015,3.65,3.65);

}

//________________________________________________________________________
AliAnalysisTaskElectronStudies::~AliAnalysisTaskElectronStudies()
{
  // default deconstructor
}
//________________________________________________________________________
void AliAnalysisTaskElectronStudies::UserCreateOutputObjects()
{
  // Create User Output Objects
  fOutputList                         = new TList();
  fOutputList->SetOwner(kTRUE);

  if(((AliConvEventCuts*)fEventCuts)->GetCutHistograms()){
    fOutputList->Add(((AliConvEventCuts*)fEventCuts)->GetCutHistograms());
  }

  if(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms()){
    fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms());
  }

    if(((AliCaloPhotonCuts*)fTMCuts)->GetCutHistograms()){
    fOutputList->Add(((AliCaloPhotonCuts*)fTMCuts)->GetCutHistograms());
  }

  if(((AliConversionPhotonCuts*)fConvCuts)->GetCutHistograms()){
    fOutputList->Add(((AliConversionPhotonCuts*)fConvCuts)->GetCutHistograms());
  }

  // for(Int_t iMatcherTask = 0; iMatcherTask < 5; iMatcherTask++){
  //     temp = 0x0;
  //     if(!fCorrTaskSetting.CompareTo("")){
  //     temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcherSignal_%i_%i",iMatcherTask,fTrackMatcherRunningMode)));
  //     } else {
  //     temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcherSignal_%i_%i_%s",iMatcherTask,fTrackMatcherRunningMode,fCorrTaskSetting.Data())));
  //     }
  //     if(temp) fOutputList->Add(temp->GetCaloTrackMatcherHistograms());
  // }
  fTrackMatcher = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager())->GetTask(fTrackMatcherName);
  if(fTrackMatcher)fOutputList->Add(fTrackMatcher->GetCaloTrackMatcherHistograms());

  fHistoNEvents           = new TH1F("NEvents","NEvents",14,-0.5,13.5);
  fHistoNEvents->GetXaxis()->SetBinLabel(1,"Accepted");
  fHistoNEvents->GetXaxis()->SetBinLabel(2,"Centrality");
  fHistoNEvents->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
  if (((AliConvEventCuts*)fEventCuts)->IsSpecialTrigger() > 1 ){
    TString TriggerNames  = "Not Trigger: ";
    TriggerNames          = TriggerNames+ ( (AliConvEventCuts*)fEventCuts)->GetSpecialTriggerName();
    fHistoNEvents->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
  } else {
    fHistoNEvents->GetXaxis()->SetBinLabel(4,"Trigger");
  }
  fHistoNEvents->GetXaxis()->SetBinLabel(5,"Vertex Z");
  fHistoNEvents->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
  fHistoNEvents->GetXaxis()->SetBinLabel(7,"Pile-Up");
  fHistoNEvents->GetXaxis()->SetBinLabel(8,"no SDD");
  fHistoNEvents->GetXaxis()->SetBinLabel(9,"no V0AND");
  fHistoNEvents->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problem");
  fHistoNEvents->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
  fHistoNEvents->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
  fHistoNEvents->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
  fHistoNEvents->GetYaxis()->SetTitle("N_{events}");
  fOutputList->Add(fHistoNEvents);

  fPtElectronTrack = new TH1F("fPtElectronTrack","fPtElectronTrack;p_{T} (GeV/c; counts",200,0.,50.);
  fOutputList->Add(fPtElectronTrack);

  fPtElectronTrackInEmcalAcc = new TH1F("fPtElectronTrackInEmcalAcc","fPtElectronTrackInEmcalAcc;p_{T} (GeV/c; counts",200,0.,50.);
  fOutputList->Add(fPtElectronTrackInEmcalAcc);

  if(fIsMC > 1){
    fHistoNEventsWOWeight           = new TH1F("NEventsWOWeight","NEventsWOWeight",14,-0.5,13.5);
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    if (((AliConvEventCuts*)fEventCuts)->IsSpecialTrigger() > 1 ){
      TString TriggerNames  = "Not Trigger: ";
      TriggerNames          = TriggerNames+ ( (AliConvEventCuts*)fEventCuts)->GetSpecialTriggerName();
      fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    } else {
      fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(5,"Vertex Z");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(7,"Pile-Up");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(10,"EMCAL problem");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fHistoNEventsWOWeight->GetYaxis()->SetTitle("N_{events}");
    fOutputList->Add(fHistoNEventsWOWeight);
  }

  if(fIsMC){
    fTruePtCluster = new TH1F("fTruePtCluster","fTruePtCluster;p_{T} (GeV/c; counts",200,0.,50.);
    fOutputList->Add(fTruePtCluster);
    fTruePtElectronCluster = new TH1F("fTruePtElectronCluster","fTruePtElectronCluster;p_{T} (GeV/c; counts",200,0.,50.);
    fOutputList->Add(fTruePtElectronCluster);
    fTruePtElectronClusterMatchedWithTrack = new TH1F("fTruePtElectronClusterMatchedWithTrack","fTruePtElectronClusterMatchedWithTrack;p_{T} (GeV/c; counts",200,0.,50.);
    fOutputList->Add(fTruePtElectronClusterMatchedWithTrack);

    fTruePtElectronTrack = new TH1F("fTruePtElectronTrack","fTruePtElectronTrack;p_{T} (GeV/c; counts",200,0.,50.);
    fOutputList->Add(fTruePtElectronTrack);
    fTruePtElectronTrackInEmcalAcc = new TH1F("fTruePtElectronTrackInEmcalAcc","fTruePtElectronTrackInEmcalAcc;p_{T} (GeV/c; counts",200,0.,50.);
    fOutputList->Add(fTruePtElectronTrackInEmcalAcc);
    fGenPtElectrons = new TH1F("fGenPtElectrons","fGenPtElectrons;p_{T} (GeV/c; counts",200,0.,50.);
    fOutputList->Add(fGenPtElectrons);
    fGenPtElectronsInEmcalAcc = new TH1F("fGenPtElectronsInEmcalAcc","fGenPtElectronsInEmcalAcc;p_{T} (GeV/c; counts",200,0.,50.);
    fOutputList->Add(fGenPtElectronsInEmcalAcc);
  }

  
  
  PostData(1, fOutputList);
  TString eventCutString = fEventCuts->GetCutNumber();
  TString clusterCutString = fClusterCutsEMC->GetCutNumber();
  TString tmCutString = fTMCuts->GetCutNumber();
  OpenFile(2);
  fAnalysisTree = new TTree(Form("AnalysisTree_%s_%s_%s_%s",eventCutString.Data(),clusterCutString.Data(),tmCutString.Data(),fCorrTaskSetting.Data()),Form("AnalysisTree_%s_%s_%s_%s",eventCutString.Data(),clusterCutString.Data(),tmCutString.Data(),fCorrTaskSetting.Data()));
  fAnalysisTree->Branch("Cluster_E", &fBuffer_ClusterE, "Cluster_E/F");
  fAnalysisTree->Branch("Cluster_M02", &fBuffer_ClusterM02, "Cluster_M02/F");
  fAnalysisTree->Branch("Cluster_M20", &fBuffer_ClusterM20, "Cluster_M20/F");
  fAnalysisTree->Branch("Track_Pt", &fBuffer_Track_Pt, "Track_Pt/F");
  fAnalysisTree->Branch("Track_P", &fBuffer_Track_P, "Track_P/F");
  fAnalysisTree->Branch("Track_dEta", &fBuffer_Track_dEta, "Track_dEta/F");
  fAnalysisTree->Branch("Track_dPhi", &fBuffer_Track_dPhi, "Track_dPhi/F");
  fAnalysisTree->Branch("Track_NSigmaElec", &fBuffer_Track_NSigmaElec, "Track_NSigmaElec/F");
  fAnalysisTree->Branch("Track_IsFromV0", &fBuffer_Track_IsFromV0, "Track_IsFromV0/I");

  
  fAnalysisTree->Branch("MatchType", &fBuffer_MatchType, "MatchType/S");
  if(fIsMC>0){    
     fAnalysisTree->Branch("MC_True_Cluster_E", &fBuffer_MC_True_Cluster_E, "MC_True_Cluster_E/F");       
     fAnalysisTree->Branch("MC_True_Track_E", &fBuffer_MC_True_Track_E, "MC_True_Track_E/F");       
     fAnalysisTree->Branch("MC_True_Track_Pt", &fBuffer_MC_True_Track_Pt, "MC_True_Track_Pt/F");       
     fAnalysisTree->Branch("MC_True_Track_P", &fBuffer_MC_True_Track_P, "MC_True_Track_P/F");       
     fAnalysisTree->Branch("MC_Track_Is_Electron", &fBuffer_MC_Track_Is_Electron, "MC_Track_Is_Electron/O");       
     fAnalysisTree->Branch("MC_Cluster_Is_Electron", &fBuffer_MC_Cluster_Is_Electron, "MC_Cluster_Is_Electron/O");       
     fAnalysisTree->Branch("MC_ClusterTrack_Same_Electron", &fBuffer_MC_ClusterTrack_Same_Electron, "MC_ClusterTrack_Same_Electron/O");       
     fAnalysisTree->Branch("MC_JetJetWeight", &fBuffer_MC_JetJetWeight, "MC_JetJetWeight/F");       
  }
  PostData(2, fAnalysisTree);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskElectronStudies::Notify()
{
    return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskElectronStudies::UserExec(Option_t *){

  fInputEvent                         = InputEvent();
  ((AliCaloPhotonCuts*)fClusterCutsEMC)->InitializeEMCAL(fInputEvent);
  //((AliCaloPhotonCuts*)fTMCuts)->InitializeEMCAL(fInputEvent);
  if(fIsMC>0) fMCEvent                  = MCEvent();
  if((fIsMC) > 0 && (!fMCEvent)) {printf("Error: No MC Event");return;}
  
  // Get V0 reader
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
  
  if(fIsMC > 0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
     RelabelAODPhotonCandidates(kTRUE);    // In case of AODMC relabeling MC
     fV0Reader->RelabelAODs(kTRUE);
  }

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(InputEvent()->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    fHistoNEvents->Fill(eventQuality);
    if (fIsMC>1) fHistoNEventsWOWeight->Fill(eventQuality);
    return;
  }
  Int_t eventNotAccepted              = fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
  if(eventNotAccepted) return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1

  if (fIsMC > 1){
      fWeightJetJetMC       = 1;
      Float_t maxjetpt      = -1.;
      Float_t pthard = -1;
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCuts)->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC ,pthard, fInputEvent, maxjetpt);
      if (!isMCJet){
        fHistoNEvents->Fill(10,fWeightJetJetMC);
        if (fIsMC>1) fHistoNEventsWOWeight->Fill(10);
        return;
      }
  }
  
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  fPIDResponse = inputHandler->GetPIDResponse(); 
  if (!fPIDResponse){AliFatal("fPIDResponse does not exist!"); return;}
 
  if (fIsMC > 1){
      fWeightJetJetMC       = 1;
      Float_t maxjetpt      = -1.;
      Float_t pthard = -1;
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCuts)->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC ,pthard, fInputEvent, maxjetpt);
      if (fIsMC == 3){
        Double_t weightMult   = ((AliConvEventCuts*)fEventCuts)->GetWeightForMultiplicity(fV0Reader->GetNumberOfPrimaryTracks());
        fWeightJetJetMC       = fWeightJetJetMC*weightMult;
      }

      if (!isMCJet){
        fHistoNEvents->Fill(10,fWeightJetJetMC);
        if (fIsMC>1) fHistoNEventsWOWeight->Fill(10);
        return;
      }
  }

  Bool_t triggered = kTRUE;
  if(eventNotAccepted!=0){
      fHistoNEvents->Fill(eventNotAccepted,fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if (fIsMC>1) fHistoNEventsWOWeight->Fill(eventNotAccepted);
      if (eventNotAccepted==3 && fIsMC > 0){
        triggered = kFALSE;
      }else {
        return;
      }
  }

  if(eventQuality != 0 && triggered== kTRUE){// Event Not Accepted
    fHistoNEvents->Fill(eventQuality, fWeightJetJetMC);
    if (fIsMC>1) fHistoNEventsWOWeight->Fill(eventQuality); // Should be 0 here
    return;
  }

  if (triggered == kTRUE) {
    fHistoNEvents->Fill(eventQuality,fWeightJetJetMC);
    if (fIsMC>1) fHistoNEventsWOWeight->Fill(eventQuality); // Should be 0 here
  }

  fGeomEMCAL                          = AliEMCALGeometry::GetInstance();
  if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}

  fBuffer_MC_JetJetWeight = fWeightJetJetMC;
  //
  // ─── MAIN PROCESSING ────────────────────────────────────────────────────────────
  //
  if (triggered==kFALSE) return;

  ProcessCaloPhotons(); // track matching is done here as well
  ProcessTracks();
  if(fIsMC>0) ProcessMCParticles();
  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
 
  if( fIsMC > 0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }
  // fill output
  PostData(2, fAnalysisTree);
  
  ResetBuffer();
  //gObjectTable->Print();
  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskElectronStudies::Terminate(Option_t *){

}

//________________________________________________________________________
void AliAnalysisTaskElectronStudies::ResetBuffer(){

}

//________________________________________________________________________
void AliAnalysisTaskElectronStudies::ProcessCaloPhotons(){
   Int_t nclus                         = 0;
   Int_t nclusCorr                     = 0;

   if(fIsMC && !fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

   TClonesArray * arrClustersProcess   = NULL;
   if(!fCorrTaskSetting.CompareTo("")){
     nclus = fInputEvent->GetNumberOfCaloClusters();

     nclusCorr = nclus;
   } else {
    arrClustersProcess                = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclusCorr                            = arrClustersProcess->GetEntries();
    nclus = fInputEvent->GetNumberOfCaloClusters();
  }
  if(nclus == 0)  return;
  // ((AliCaloPhotonCuts*)fClusterCutsEMC)->FillHistogramsExtendedQA(fInputEvent,fIsMC);
  // ((AliCaloPhotonCuts*)fClusterCutsPHOS)->FillHistogramsExtendedQA(fInputEvent,fIsMC);
  // in case user wants to use default track matching
  fTMCuts->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
  AliAODCaloCluster* clus                       = NULL;   
  if(arrClustersProcess){ 
     // EMCal correction framework was used
     // we need to loop over this for EMCal clusters and 
     // over the others for PHOS cluster

      // Loop over EMCal clusters
      for(Long_t i = 0; i < nclusCorr; i++){
        Double_t tempClusterWeight              =  fWeightJetJetMC;
        clus                                = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(i));
        
        if(!clus) continue;
        if ( !clus->IsEMCAL()){ // for PHOS: cluster->GetType() == AliVCluster::kPHOSNeutral
          delete clus;
          continue;
        }
        // check if given EMC cuts are fulfilled
        if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
          delete clus;
          continue;
        }
        
        if(fIsMC){
          Int_t *mclabelsCluster = clus->GetLabels();
          if (clus->GetNLabels() > 0)
          {
            if(mclabelsCluster[0]!=-1){
              AliAODMCParticle* clusterMother = (AliAODMCParticle* )fAODMCTrackArray->At(mclabelsCluster[0]);
              fTruePtCluster->Fill(clusterMother->Pt(),fWeightJetJetMC);
              if (TMath::Abs(clusterMother->PdgCode()) == 11) {
                  fTruePtElectronCluster->Fill(clusterMother->Pt(),fWeightJetJetMC);
              }
            }
          }
        }
        AliAODTrack *aodt = NULL;
        //ProcessTrackMatching(clus);
        if(fTMCuts->CheckClusterForTrackMatch(clus)){
           Int_t labelTrackClosest = -1;
           if(fTMCuts->GetClosestMatchedTrackToCluster(fInputEvent,clus,labelTrackClosest)){
              Int_t properLabel = labelTrackClosest;
              // if(labelTrack<0) properLabel = (-1 * labelTrack) - 1; // conversion hybrid none hybrid
              aodt = dynamic_cast<AliAODTrack*> (fInputEvent->GetTrack(properLabel));
              if(TrackIsSelectedAOD(aodt)){
                 ProcessMatchedTrack(aodt,clus,kFALSE);
              }
           }
        }

        delete clus;
      } // end of initial cluster loop
  }

  
  // no need to loop over normal clusters as well
  if(arrClustersProcess) return; 

  // Loop over normal clusters
  for(Long_t i = 0; i < nclus; i++){
    Double_t tempClusterWeight              =  fWeightJetJetMC;                   
    clus                                = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    
    if(!clus) continue;

    if(clus->IsEMCAL()){ // if is was not saved already
      if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        delete clus;
        continue;
      }

      if(fIsMC){
          Int_t *mclabelsCluster = clus->GetLabels();
          if (clus->GetNLabels() > 0)
          {
            if(mclabelsCluster[0]!=-1){
              AliAODMCParticle* clusterMother = (AliAODMCParticle* )fAODMCTrackArray->At(mclabelsCluster[0]);
              fTruePtCluster->Fill(clusterMother->Pt(),fWeightJetJetMC);
              if (TMath::Abs(clusterMother->PdgCode()) == 11) {
                  fTruePtElectronCluster->Fill(clusterMother->Pt(),fWeightJetJetMC);
              }
            }
          }
        }
        AliAODTrack *aodt = NULL;
        //ProcessTrackMatching(clus);
        if(fTMCuts->CheckClusterForTrackMatch(clus)){
           Int_t labelTrack = -1;
           if(fTMCuts->GetClosestMatchedTrackToCluster(fInputEvent,clus,labelTrack)){
              Int_t properLabel = labelTrack;
              // if(labelTrack<0) properLabel = (-1 * labelTrack) - 1; // conversion hybrid none hybrid
              aodt = dynamic_cast<AliAODTrack*> (fInputEvent->GetTrack(properLabel));
              if(TrackIsSelectedAOD(aodt)) ProcessMatchedTrack(aodt,clus,kFALSE);
           }
        }

      //ProcessTrackMatching(clus); // tree filling done here too
      
      delete clus;
      continue;
    }
    delete clus;
  }  
}

//________________________________________________________________________
void AliAnalysisTaskElectronStudies::ProcessTracks(){
   for(Int_t t=0;t<fInputEvent->GetNumberOfTracks();t++){
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*> (fInputEvent->GetTrack(t));
      if(!aodt) continue;
      if(!TrackIsSelectedAOD(aodt)) continue;

      // apply electron PID cut
      Float_t nsigmaelec = fPIDResponse->NumberOfSigmasTPC(aodt,AliPID::kElectron); 
      if ((nsigmaelec > fMaxNsigmaElec) || (nsigmaelec < fMinNsigmaElec)) continue;
     
      fPtElectronTrack->Fill(aodt->Pt(),fWeightJetJetMC);
      if(IsInEMCalAcceptance(aodt)) fPtElectronTrackInEmcalAcc->Fill(aodt->Pt(),fWeightJetJetMC);

      if(fIsMC){
         Int_t mclabel = aodt->GetLabel();
         if(mclabel != -1){
            if(mclabel<0) mclabel = (-1 * mclabel) - 1;
            AliAODMCParticle* particle =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(mclabel));
            if(particle){
               if(TMath::Abs(particle->GetPdgCode())== 11){
                  fTruePtElectronTrack->Fill(particle->Pt(),fWeightJetJetMC);
                  if(IsInEMCalAcceptance(particle)) fTruePtElectronTrackInEmcalAcc->Fill(particle->Pt(),fWeightJetJetMC);
               }
            }
         }
      }

   }
}

///________________________________________________________________________
Bool_t AliAnalysisTaskElectronStudies::TrackIsSelectedAOD(AliAODTrack* lTrack) {
  // apply filter bits 
  if(! lTrack) return kFALSE;
  if( ! lTrack->IsHybridGlobalConstrainedGlobal()){
    return kFALSE;
  }

	// Absolute TPC Cluster cut
	if(lTrack->GetTPCNcls()<fMinClsTPC) return kFALSE;
	if(lTrack->GetTPCchi2perCluster()>fChi2PerClsTPC) return kFALSE;

  // Found / findable cluster in TPC cuts
  Double_t clsToF=0;
  if(lTrack->GetTPCNclsF()!=0){
    clsToF = (Double_t)lTrack->GetNcls(1)/(Double_t)lTrack->GetTPCNclsF();
  }
  if(clsToF < fMinFracClsTPC) return kFALSE;
  
  // DCA cut 
  Float_t b[2];
  Float_t bCov[3];
  lTrack->GetImpactParameters(b,bCov);
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }

  Float_t dcaToVertexXY = b[0];
  Float_t dcaToVertexZ = b[1];

  if(dcaToVertexXY > fMaxDCAxy) return kFALSE;
  if(dcaToVertexZ > fMaxDCAz) return kFALSE;


  // ITS Cluster Cut
	// SetClusterRequirementITS and SetRequireITSRefit can
	// not be set for AODs after filtering
 
	if(lTrack->GetITSNcls()<fMinClsITS) return kFALSE;

  if(  TMath::Abs(lTrack->Eta()) > fEtaCut ) {
      return kFALSE;
  }

  if( lTrack->Pt() < fPtCut ) {
    return kFALSE;
  }

  // n sigma preselection
  Float_t nsigmaelec = fPIDResponse->NumberOfSigmasTPC(lTrack,AliPID::kElectron); 
  if ((nsigmaelec > 10) || (nsigmaelec < -10))
    return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskElectronStudies::ProcessMatchedTrack(AliAODTrack* track, AliAODCaloCluster* clus, Bool_t isV0){
    Float_t clusPos[3]                      = { 0,0,0 };
    clus->GetPosition(clusPos);
    TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
    
    fBuffer_ClusterE = clus->E();
    fBuffer_ClusterM02 = clus->GetM02();
    fBuffer_ClusterM20 = clus->GetM20();
    fBuffer_Track_Pt = track->Pt();
    fBuffer_Track_P = track->P();
    fBuffer_Track_NSigmaElec = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron); 
    fBuffer_Track_IsFromV0 = isV0; 
    
    fBuffer_MC_True_Cluster_E = 0; 
    fBuffer_MC_True_Track_E = 0; 
    fBuffer_MC_True_Track_Pt = 0; 
    fBuffer_MC_True_Track_P = 0; 
    fBuffer_MC_Track_Is_Electron= kFALSE;
    fBuffer_MC_Cluster_Is_Electron= kFALSE;
    fBuffer_MC_ClusterTrack_Same_Electron= kFALSE;
    
    Float_t tempEta = -99999;
    Float_t tempPhi = -99999;
    // Get dEta and dPhi of track on EMCal surface!
    ((AliCaloTrackMatcher*) fTMCuts->GetCaloTrackMatcherInstance())->GetTrackClusterMatchingResidual(track->GetID(),clus->GetID(),tempEta,tempPhi);
    fBuffer_Track_dEta = tempEta; 
    fBuffer_Track_dPhi = tempPhi;
    
    if(fIsMC){
        // check if leading contribution is electron
        Int_t *mclabelsCluster = clus->GetLabels();

        AliAODMCParticle* clusterMother = NULL;
        if (clus->GetNLabels() > 0)
        {
        // for (Int_t k = 0; k < (Int_t)clusterE->GetNLabels(); k++)
        // {
          if(mclabelsCluster[0]!=-1){
            clusterMother = (AliAODMCParticle* )fAODMCTrackArray->At(mclabelsCluster[0]);
            if (TMath::Abs(clusterMother->PdgCode()) == 11) {
                fBuffer_MC_Cluster_Is_Electron = kTRUE;
                fBuffer_MC_True_Cluster_E = clusterMother->E(); 
                fTruePtElectronClusterMatchedWithTrack->Fill(clusterMother->Pt(),fWeightJetJetMC);
            }
          }
        }
  
            
        // Check Track
        Int_t trackMCLabel = track->GetLabel();
        AliAODMCParticle* trackMother = NULL;
        if(trackMCLabel>-1){
          trackMother = (AliAODMCParticle* )fAODMCTrackArray->At(trackMCLabel);
          if(TMath::Abs(trackMother->GetPdgCode()) == 11){
            fBuffer_MC_Track_Is_Electron = kTRUE;
            fBuffer_MC_True_Track_E = trackMother->E(); 
            fBuffer_MC_True_Track_Pt = trackMother->Pt(); 
            fBuffer_MC_True_Track_P = trackMother->P();
          }
        }

        if((clus->GetNLabels()>0) && (mclabelsCluster[0] == trackMCLabel) && (fBuffer_MC_Track_Is_Electron ==1)){
          fBuffer_MC_ClusterTrack_Same_Electron = kTRUE;
        }
          
      // }
    } // end is MC

    AliAODTrack* highTrack = NULL;
    // Check other case
    Int_t labelTrackHighest = -1;
    if(fTMCuts->GetHighestPtMatchedTrackToCluster(fInputEvent,clus,labelTrackHighest)){
          Int_t properLabel = labelTrackHighest;
          // if(labelTrack<0) properLabel = (-1 * labelTrack) - 1; // conversion hybrid none hybrid
          highTrack = dynamic_cast<AliAODTrack*> (fInputEvent->GetTrack(properLabel));
          if(TrackIsSelectedAOD(highTrack)){
            if(highTrack->GetID()!=track->GetID()){
              // found track that is different
              fBuffer_MatchType = 1;
              fAnalysisTree->Fill();

              // Fill high pT track stuff
              fBuffer_MatchType = 2;
              Float_t tempEtaHigh = -99999;
              Float_t tempPhiHigh = -99999;
              ((AliCaloTrackMatcher*) fTMCuts->GetCaloTrackMatcherInstance())->GetTrackClusterMatchingResidual(highTrack->GetID(),clus->GetID(),tempEtaHigh,tempPhiHigh);
              fBuffer_Track_dEta = tempEtaHigh; 
              fBuffer_Track_dPhi = tempPhiHigh;

              fBuffer_Track_Pt = highTrack->Pt();
              fBuffer_Track_P = highTrack->P();
              fBuffer_Track_NSigmaElec = fPIDResponse->NumberOfSigmasTPC(highTrack,AliPID::kElectron); 
              fBuffer_Track_IsFromV0 = isV0; 
              
              fBuffer_MC_True_Cluster_E = 0; 
              fBuffer_MC_True_Track_E = 0; 
              fBuffer_MC_True_Track_Pt = 0; 
              fBuffer_MC_True_Track_P = 0; 
              fBuffer_MC_Track_Is_Electron= kFALSE;
              fBuffer_MC_Cluster_Is_Electron= kFALSE;
              fBuffer_MC_ClusterTrack_Same_Electron= kFALSE;

              if(fIsMC){
                  // check if leading contribution is electron
                  Int_t *mclabelsCluster = clus->GetLabels();

                  AliAODMCParticle* clusterMother = NULL;
                  if (clus->GetNLabels() > 0)
                  {
                  // for (Int_t k = 0; k < (Int_t)clusterE->GetNLabels(); k++)
                  // {
                    if(mclabelsCluster[0]!=-1){
                      clusterMother = (AliAODMCParticle* )fAODMCTrackArray->At(mclabelsCluster[0]);
                      if (TMath::Abs(clusterMother->PdgCode()) == 11) {
                          fBuffer_MC_Cluster_Is_Electron = kTRUE;
                          fBuffer_MC_True_Cluster_E = clusterMother->E(); 
                          fTruePtElectronClusterMatchedWithTrack->Fill(clusterMother->Pt(),fWeightJetJetMC);
                      }
                    }
                  }
            
                      
                  // Check Track
                  Int_t trackMCLabel = track->GetLabel();
                  AliAODMCParticle* trackMother = NULL;
                  if(trackMCLabel>-1){
                    trackMother = (AliAODMCParticle* )fAODMCTrackArray->At(trackMCLabel);
                    if(TMath::Abs(trackMother->GetPdgCode()) == 11){
                      fBuffer_MC_Track_Is_Electron = kTRUE;
                      fBuffer_MC_True_Track_E = trackMother->E(); 
                      fBuffer_MC_True_Track_Pt = trackMother->Pt(); 
                      fBuffer_MC_True_Track_P = trackMother->P();
                    }
                  }

                  if((clus->GetNLabels()>0) && (mclabelsCluster[0] == trackMCLabel) && (fBuffer_MC_Track_Is_Electron ==1)){
                    fBuffer_MC_ClusterTrack_Same_Electron = kTRUE;
                  }
                    
                // }
              } // end is MC
              fAnalysisTree->Fill();
            } else{
              fBuffer_MatchType = 0;
              fAnalysisTree->Fill();
            } 
          }
    } else{
      // did not work or did not find
      fBuffer_MatchType = 0;
      fAnalysisTree->Fill();
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskElectronStudies::ProcessMCParticles(){
  // Loop over all primary MC particle
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray){
    for(Int_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {
      AliAODMCParticle* particle =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
      if(!particle) continue;
      if(particle->MCStatusCode() != 1 || !particle->IsPhysicalPrimary()) continue;   
      if(TMath::Abs(particle->PdgCode()) != 11 ) continue;
      fGenPtElectrons->Fill(particle->Pt(),fWeightJetJetMC);
      if(IsInEMCalAcceptance(particle)){
        fGenPtElectronsInEmcalAcc->Fill(particle->Pt(),fWeightJetJetMC);
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskElectronStudies::ProcessTrackMatching(AliAODCaloCluster* clus){
     Int_t nModules = fGeomEMCAL->GetNumberOfSuperModules();
     AliExternalTrackParam *trackParam = 0;
     for(Int_t t=0;t<fInputEvent->GetNumberOfTracks();t++){
        AliAODTrack *aodt = dynamic_cast<AliAODTrack*> (fInputEvent->GetTrack(t));
        if(!aodt) continue;
        if(!TrackIsSelectedAOD(aodt)) continue;

        Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
        aodt->GetPxPyPz(pxpypz);
        aodt->GetXYZ(xyz);
        aodt->GetCovarianceXYZPxPyPz(cv);

        trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,aodt->Charge());
     
        AliExternalTrackParam emcParam(*trackParam);
        Float_t eta, phi, pt;
        //propagate tracks to emc surfaces
        if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
          delete trackParam;
          continue;
        }
        if( TMath::Abs(eta) > 0.75 ) {
          delete trackParam;
          continue;
        }
        // Save some time and memory in case of no DCal present
        if( nModules < 13 && ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())){
          delete trackParam;
          continue;
        }
        // Save some time and memory in case of run2
        if( nModules > 12 ){
          if (( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad()) && ( phi < 250*TMath::DegToRad() || phi > 340*TMath::DegToRad())){
            delete trackParam;
            continue;
          }
        }
        Float_t dEta=-999, dPhi=-999;
        Double_t trkPos[3] = {0.,0.,0.};
        if (!emcParam.GetXYZ(trkPos)){
          delete trackParam;
          continue;
        }

        AliExternalTrackParam trackParamTmp(emcParam);//Retrieve the starting point every time before the extrapolation
        if(!AliEMCALRecoUtils::ExtrapolateTrackToCluster(&trackParamTmp, clus, 0.139, 5., dEta, dPhi)){
             delete trackParam;
             continue;
        }
        if(!fUseRTrackMatching){
          if(TMath::Abs(dEta) > (fMatchingParamsEta[0] + pow(aodt->Pt() + fMatchingParamsEta[1],fMatchingParamsEta[2]))){
              delete trackParam;
              continue;
          }
          if(TMath::Abs(dPhi) > (fMatchingParamsPhi[0] + pow(aodt->Pt() + fMatchingParamsPhi[1],fMatchingParamsPhi[2]))){
              delete trackParam;
              continue;
          }
        } else{ // use R track matching
            Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
            if(dR > fRTrackMatching){
              delete trackParam;
              continue;
            }

        }
        // track is matched

        ProcessMatchedTrack(aodt,clus,kFALSE);
        delete trackParam;

     }

    // check also conversion sample to be sure nothing from there is missing
     if(!fReaderGammas) fReaderGammas    = fV0Reader->GetReconstructedGammas();
     for (Int_t c = 0; c < fReaderGammas->GetEntriesFast(); c++)
     {
        AliAODConversionPhoton* photon = (AliAODConversionPhoton*) fReaderGammas->At(c);
        if(!photon) continue;
        if(!((AliConversionPhotonCuts*)fConvCuts)->PhotonIsSelected(photon,fInputEvent)) continue;
        
        for (Int_t iElec = 0;iElec < 2;iElec++){
          Int_t tracklabel = photon->GetLabel(iElec);
          AliAODTrack *convtrack = dynamic_cast<AliAODTrack*> (fInputEvent->GetTrack(tracklabel));
          if(!convtrack) continue;
          if(convtrack->IsHybridGlobalConstrainedGlobal()) continue; // that means we already treated it;

          Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
          convtrack->GetPxPyPz(pxpypz);
          convtrack->GetXYZ(xyz);
          convtrack->GetCovarianceXYZPxPyPz(cv);

          trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,convtrack->Charge());
      
          AliExternalTrackParam emcParam(*trackParam);
          Float_t eta, phi, pt;
          //propagate tracks to emc surfaces
          if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
            delete trackParam;
            continue;
          }
          if( TMath::Abs(eta) > 0.75 ) {
            delete trackParam;
            continue;
          }
          // Save some time and memory in case of no DCal present
          if( nModules < 13 && ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())){
            delete trackParam;
            continue;
          }
          // Save some time and memory in case of run2
          if( nModules > 12 ){
            if (( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad()) && ( phi < 250*TMath::DegToRad() || phi > 340*TMath::DegToRad())){
              delete trackParam;
              continue;
            }
          }
          Float_t dEta=-999, dPhi=-999;
          Double_t trkPos[3] = {0.,0.,0.};
          if (!emcParam.GetXYZ(trkPos)){
            delete trackParam;
            continue;
          }

          AliExternalTrackParam trackParamTmp(emcParam);//Retrieve the starting point every time before the extrapolation
          if(!AliEMCALRecoUtils::ExtrapolateTrackToCluster(&trackParamTmp, clus, 0.139, 5., dEta, dPhi)){
               delete trackParam;
               continue;
          }
          if(TMath::Abs(dEta) > (fMatchingParamsEta[0] + pow(convtrack->Pt() + fMatchingParamsEta[1],fMatchingParamsEta[2]))){
               delete trackParam;
               continue;
          }
          if(TMath::Abs(dPhi) > (fMatchingParamsPhi[0] + pow(convtrack->Pt() + fMatchingParamsPhi[1],fMatchingParamsPhi[2]))){
               delete trackParam;
               continue;
          }
          
          // track is matched
          ProcessMatchedTrack(convtrack,clus,kTRUE);
          delete trackParam;
        }
     }
    //  return highestMatchIndex; 
}
Bool_t AliAnalysisTaskElectronStudies::IsSameTrack(Int_t id1, Int_t id2){
    if((id1 == -999) || (id2 == -999)){
        cout << "ERROR: Track info is missing for one track!" << endl;
        return kFALSE;
    }
    Int_t esdID1 = id1; 
    Int_t esdID2 = id2; 
    if(id1<0) esdID1 = (-1 * id1) - 1;
    if(id2<0) esdID2 = (-1 * id2) - 1;
    if(esdID1 == esdID2){
        return kTRUE;
    } else{
        return kFALSE;
    }

}

Bool_t AliAnalysisTaskElectronStudies::IsInEMCalAcceptance(AliAODConversionPhoton *photon)
{
    Double_t eta = photon->GetPhotonEta();
    Double_t phi = photon->GetPhotonPhi();
    // cout << phi << endl;
    // cout << eta << endl;
    if (phi < 0)
        phi += 2 * TMath::Pi();
    if ((eta < -0.6687) || (eta > 0.66465))
        return kFALSE;
    if ((phi < 1.39626) || (phi > 3.15))
        return kFALSE;
    return kTRUE;
}
Bool_t AliAnalysisTaskElectronStudies::IsInEMCalAcceptance(AliAODMCParticle *part)
{
    Double_t eta = part->Eta();
    Double_t phi = part->Phi();

    if (phi < 0)
        phi += 2 * TMath::Pi();
    // cout << phi << endl;
    // cout << eta << endl;
    if ((eta < -0.6687) || (eta > 0.66465))
        return kFALSE;
    if ((phi < 1.39626) || (phi > 3.15))
        return kFALSE;
    return kTRUE;
}
Bool_t AliAnalysisTaskElectronStudies::IsInEMCalAcceptance(AliAODTrack *part)
{
    Double_t eta = part->Eta();
    Double_t phi = part->Phi();

    if (phi < 0)
        phi += 2 * TMath::Pi();
    // cout << phi << endl;
    // cout << eta << endl;
    if ((eta < -0.6687) || (eta > 0.66465))
        return kFALSE;
    if ((phi < 1.39626) || (phi > 3.15))
        return kFALSE;
    return kTRUE;
}

Int_t AliAnalysisTaskElectronStudies::CheckClustersForMCContribution(Int_t mclabel, TClonesArray *vclus)
{
    Int_t clusterLabel = -1; // position of cluster in array where mc label was found as contribution
    for (Int_t p = 0; p < vclus->GetEntriesFast(); p++)
    {
        AliAODCaloCluster *clus = (AliAODCaloCluster *)vclus->At(p);
        if (!clus)
            continue;
        Int_t *mclabelsCluster = clus->GetLabels();
        if (clus->GetNLabels() > 0)
        {
            for (Int_t k = 0; k < (Int_t)clus->GetNLabels(); k++)
            {
                if (mclabelsCluster[p] == mclabel)
                    clusterLabel = p;
            }
        }
    }

    return clusterLabel;
}
//________________________________________________________________________
void AliAnalysisTaskElectronStudies::RelabelAODPhotonCandidates(Bool_t mode){

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel

  Int_t* fMCEventPos = nullptr;
  Int_t* fMCEventNeg = nullptr;
  Int_t* fESDArrayPos = nullptr;
  Int_t* fESDArrayNeg = nullptr;
  if(mode){
    fMCEventPos = new Int_t[fReaderGammas->GetEntries()];
    fMCEventNeg = new Int_t[fReaderGammas->GetEntries()];
    fESDArrayPos = new Int_t[fReaderGammas->GetEntries()];
    fESDArrayNeg = new Int_t[fReaderGammas->GetEntries()];
  }

  for(Int_t iGamma = 0;iGamma<fReaderGammas->GetEntries();iGamma++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(iGamma);
    if(!PhotonCandidate) continue;
    if(!mode){// Back to ESD Labels
      PhotonCandidate->SetMCLabelPositive(fMCEventPos[iGamma]);
      PhotonCandidate->SetMCLabelNegative(fMCEventNeg[iGamma]);
      PhotonCandidate->SetLabelPositive(fESDArrayPos[iGamma]);
      PhotonCandidate->SetLabelNegative(fESDArrayNeg[iGamma]);
      continue;
    }
    fMCEventPos[iGamma] =  PhotonCandidate->GetMCLabelPositive();
    fMCEventNeg[iGamma] =  PhotonCandidate->GetMCLabelNegative();
    fESDArrayPos[iGamma] = PhotonCandidate->GetTrackLabelPositive();
    fESDArrayNeg[iGamma] = PhotonCandidate->GetTrackLabelNegative();

    Bool_t AODLabelPos = kFALSE;
    Bool_t AODLabelNeg = kFALSE;

    for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
      AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
      if(!AODLabelPos){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
        PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelPositive(i);
        AODLabelPos = kTRUE;
        }
      }
      if(!AODLabelNeg){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
        PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelNegative(i);
        AODLabelNeg = kTRUE;
        }
      }
      if(AODLabelNeg && AODLabelPos){
        break;
      }
    }
    if(!AODLabelPos || !AODLabelNeg){
      cout<<"WARNING!!! AOD TRACKS NOT FOUND FOR"<<endl;
      if(!AODLabelNeg){
        PhotonCandidate->SetMCLabelNegative(-999999);
        PhotonCandidate->SetLabelNegative(-999999);
      }
      if(!AODLabelPos){
        PhotonCandidate->SetMCLabelPositive(-999999);
        PhotonCandidate->SetLabelPositive(-999999);
      }
    }
  }


  if(!mode){
    delete[] fMCEventPos;
    delete[] fMCEventNeg;
    delete[] fESDArrayPos;
    delete[] fESDArrayNeg;
  }
}




