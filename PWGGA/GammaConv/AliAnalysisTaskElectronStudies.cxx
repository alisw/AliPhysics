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
// #include "TParticle.h"
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
#include "AliEMCALRecoUtils.h"
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
  fEMCalRecoUtils(NULL),
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
  fTrackPvsPOnSurface(NULL),
  fTrackPvsPOnSurfaceOwn(NULL),
  fTrackRefPvsR(NULL),
  fTrackPOnSurface(NULL),
  fTrackPOnSurfaceOwn(NULL),
  fTrackPOnSurfaceTrue(NULL),
  fTreeBuffSize(60*1024*1024),
  fMemCountAOD(0),
  fTrackMatcherRunningMode(0),
  fIsoMaxRadius(0.2),
  fConversionTrackMatchR(0.4),
  fBuffer_NPrimaryTracks(0), 
  fBuffer_NClus(0), 
  fBuffer_IsProblem(kFALSE), 
  fBuffer_ClusterE(0), 
  fBuffer_ClusterM02(0), 
  fBuffer_ClusterM20(0), 
  fBuffer_ClusterNCells(0), 
  fBuffer_Track_E(0), 
  fBuffer_Track_Px(0), 
  fBuffer_Track_Py(0), 
  fBuffer_Track_Pz(0), 
  fBuffer_Track_PonEMCal(0), 
  fBuffer_Track_Charge(0), 
  fBuffer_Track_dEta(0), 
  fBuffer_Track_dPhi(0), 
  fBuffer_Track_NSigmaElec(0), 
  fBuffer_Track_IsFromV0(0), 
  fBuffer_Track_ClosestR(0), 
  fBuffer_Track_ChargedIso(0), 
  fBuffer_MatchType(0), 
  fBuffer_MC_True_Cluster_E(0), 
  fBuffer_MC_True_Track_E(0), 
  fBuffer_MC_True_Track_Px(0), 
  fBuffer_MC_True_Track_Py(0), 
  fBuffer_MC_True_Track_Pz(0), 
  fBuffer_MC_True_Track_MotherPDG(0), 
  fBuffer_MC_Track_Is_Electron(0), 
  fBuffer_MC_Cluster_Is_Electron(0), 
  fBuffer_MC_ClusterTrack_Same_Electron(0), 
  fBuffer_MC_JetJetWeight(1),
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
  fEMCalRecoUtils(NULL),
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
  fTrackPvsPOnSurface(NULL),
  fTrackPvsPOnSurfaceOwn(NULL),
  fTrackRefPvsR(NULL),
   fTrackPOnSurface(NULL),
  fTrackPOnSurfaceOwn(NULL),
  fTrackPOnSurfaceTrue(NULL),
  fTreeBuffSize(60*1024*1024),
  fMemCountAOD(0),
  fTrackMatcherRunningMode(0),
  fIsoMaxRadius(0.2),
  fConversionTrackMatchR(0.4),
  fBuffer_NPrimaryTracks(0), 
  fBuffer_NClus(0), 
  fBuffer_IsProblem(kFALSE), 
  fBuffer_ClusterE(0), 
  fBuffer_ClusterM02(0), 
  fBuffer_ClusterM20(0), 
  fBuffer_ClusterNCells(0), 
  fBuffer_Track_E(0), 
  fBuffer_Track_Px(0), 
  fBuffer_Track_Py(0), 
  fBuffer_Track_Pz(0), 
  fBuffer_Track_PonEMCal(0), 
  fBuffer_Track_Charge(0), 
  fBuffer_Track_dEta(0), 
  fBuffer_Track_dPhi(0), 
  fBuffer_Track_NSigmaElec(0), 
  fBuffer_Track_IsFromV0(0), 
  fBuffer_Track_ClosestR(0), 
  fBuffer_Track_ChargedIso(0), 
  fBuffer_MatchType(0), 
  fBuffer_MC_True_Cluster_E(0), 
  fBuffer_MC_True_Track_E(0), 
  fBuffer_MC_True_Track_Px(0), 
  fBuffer_MC_True_Track_Py(0), 
  fBuffer_MC_True_Track_Pz(0), 
  fBuffer_MC_True_Track_MotherPDG(0), 
  fBuffer_MC_Track_Is_Electron(0), 
  fBuffer_MC_Cluster_Is_Electron(0), 
  fBuffer_MC_ClusterTrack_Same_Electron(0), 
  fBuffer_MC_JetJetWeight(1.),
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
  fEMCalRecoUtils  = new AliEMCALRecoUtils;
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

  fTrackPvsPOnSurface =  new TH2F("fTrackPvsPOnSurface","fTrackPvsPOnSurface;P at vertex (GeV/c); P on EMCal surface",100,0.,10.,100,0,10);
  fOutputList->Add(fTrackPvsPOnSurface);

  fTrackPvsPOnSurfaceOwn =  new TH2F("fTrackPvsPOnSurfaceOwn","fTrackPvsPOnSurfaceOwn;P at vertex (GeV/c); P on EMCal surface (with Bremsstrahlung)",100,0.,10.,100,0,10);
  fOutputList->Add(fTrackPvsPOnSurfaceOwn);

  fTrackRefPvsR =  new TH2F("fTrackRefPvsR","fTrackRefPvsR;R (cm); P (GeV/c)",48,0.,480,100,0,10);
  fOutputList->Add(fTrackRefPvsR);
  
  fTrackPOnSurface =  new TH1F("fTrackPOnSurface","fTrackPOnSurface;P on EMCal surface ",100,0,10);
  fTrackPOnSurfaceOwn =  new TH1F("fTrackPOnSurfaceOwn","fTrackPOnSurfaceOwn;P on EMCal surface ",100,0,10);
  fTrackPOnSurfaceTrue =  new TH1F("fTrackPOnSurfaceTrue","fTrackPOnSurfaceTrue;P on EMCal surface ",100,0,10);
  fOutputList->Add(fTrackPOnSurface);
  fOutputList->Add(fTrackPOnSurfaceOwn);
  fOutputList->Add(fTrackPOnSurfaceTrue);
  
  PostData(1, fOutputList);
  TString eventCutString = fEventCuts->GetCutNumber();
  TString clusterCutString = fClusterCutsEMC->GetCutNumber();
  TString tmCutString = fTMCuts->GetCutNumber();
  OpenFile(2);
  fAnalysisTree = new TTree(Form("AnalysisTree_%s_%s_%s_%s",eventCutString.Data(),clusterCutString.Data(),tmCutString.Data(),fCorrTaskSetting.Data()),Form("AnalysisTree_%s_%s_%s_%s",eventCutString.Data(),clusterCutString.Data(),tmCutString.Data(),fCorrTaskSetting.Data()));
  
  fAnalysisTree->Branch("Event_NPrimaryTracks", &fBuffer_NPrimaryTracks,"Event_NPrimaryTracks/s");
  fAnalysisTree->Branch("Event_NClus", &fBuffer_NClus,"Event_NClus/s");
  fAnalysisTree->Branch("Event_IsProblem", &fBuffer_IsProblem,"Event_IsProblem/O");
  fAnalysisTree->Branch("Cluster_E","std::vector<UShort_t>",&fBuffer_ClusterE);
  fAnalysisTree->Branch("Cluster_M02","std::vector<UShort_t>", &fBuffer_ClusterM02);
  fAnalysisTree->Branch("Cluster_M20","std::vector<UShort_t>", &fBuffer_ClusterM20);
  fAnalysisTree->Branch("Cluster_NCells","std::vector<UShort_t>", &fBuffer_ClusterNCells);
  fAnalysisTree->Branch("Track_E","std::vector<UShort_t>", &fBuffer_Track_E);
  fAnalysisTree->Branch("Track_Px","std::vector<Short_t>", &fBuffer_Track_Px);
  fAnalysisTree->Branch("Track_Py","std::vector<Short_t>", &fBuffer_Track_Py);
  fAnalysisTree->Branch("Track_Pz","std::vector<Short_t>", &fBuffer_Track_Pz);
  fAnalysisTree->Branch("Track_PonEMCal","std::vector<UShort_t>", &fBuffer_Track_PonEMCal);
  fAnalysisTree->Branch("Track_Charge","std::vector<Short_t>", &fBuffer_Track_Charge);
  fAnalysisTree->Branch("Track_dEta","std::vector<Short_t>", &fBuffer_Track_dEta);
  fAnalysisTree->Branch("Track_dPhi","std::vector<Short_t>", &fBuffer_Track_dPhi);
  fAnalysisTree->Branch("Track_NSigmaElec","std::vector<Short_t>", &fBuffer_Track_NSigmaElec);
  fAnalysisTree->Branch("Track_IsFromV0","std::vector<Bool_t>", &fBuffer_Track_IsFromV0);
  fAnalysisTree->Branch("Track_ClosestR","std::vector<UShort_t>", &fBuffer_Track_ClosestR);
  fAnalysisTree->Branch("Track_ChargedIso","std::vector<UShort_t>", &fBuffer_Track_ChargedIso);
  fAnalysisTree->Branch("Track_MatchType","std::vector<UShort_t>", &fBuffer_MatchType);

  
  if(fIsMC>0){    
     fAnalysisTree->Branch("MC_True_Cluster_E","std::vector<UShort_t>", &fBuffer_MC_True_Cluster_E);
     fAnalysisTree->Branch("MC_True_Track_E","std::vector<UShort_t>", &fBuffer_MC_True_Track_E);
     fAnalysisTree->Branch("MC_True_Track_Px","std::vector<Short_t>", &fBuffer_MC_True_Track_Px);
     fAnalysisTree->Branch("MC_True_Track_Py","std::vector<Short_t>", &fBuffer_MC_True_Track_Py);
     fAnalysisTree->Branch("MC_True_Track_Pz","std::vector<Short_t>", &fBuffer_MC_True_Track_Pz);
     fAnalysisTree->Branch("MC_True_Track_MotherPDG","std::vector<Int_t>", &fBuffer_MC_True_Track_MotherPDG);
     fAnalysisTree->Branch("MC_Track_Is_Electron","std::vector<Bool_t>", &fBuffer_MC_Track_Is_Electron);
     fAnalysisTree->Branch("MC_Cluster_Is_Electron","std::vector<Bool_t>", &fBuffer_MC_Cluster_Is_Electron);
     fAnalysisTree->Branch("MC_ClusterTrack_Same_Electron","std::vector<UShort_t>", &fBuffer_MC_ClusterTrack_Same_Electron); // 1 for leading label of clus electron, 2 for sub leading
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

  if(fInputEvent->IsA()==AliESDEvent::Class()){
     // do only tracking studies and quit
     ProcessTracksESD();
     return;
  }

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
  fAnalysisTree->Fill();

  PostData(2, fAnalysisTree);

  // cout << "-------- END OF EVENT" << endl;
  // cout << "Size vector = " << fBuffer_ClusterE.size() << endl;
  // for (Int_t i = 0; i < fBuffer_ClusterE.size(); i++)
  // {
  //   cout << "i=" << i << "   ClusterE=" << fBuffer_ClusterE.at(i) << endl;
  //   cout << "i=" << i << "   MatchType=" << fBuffer_MatchType.at(i) << endl;
  //   cout << "i=" << i << "   TrackPx=" << (Float_t)fBuffer_Track_Px.at(i)/kShortScaleLow << endl;
  //   cout << "i=" << i << "   TrackPy=" << (Float_t)fBuffer_Track_Py.at(i)/kShortScaleLow  << endl;
  //   cout << "i=" << i << "   TrackPz=" << (Float_t)fBuffer_Track_Pz.at(i)/kShortScaleLow  << endl;
  //   cout << "i=" << i << "   TrackdEta=" << (Float_t)fBuffer_Track_dEta.at(i)/kShortScaleHigh << endl;
  //   cout << "i=" << i << "   TrackdPhi=" << (Float_t)fBuffer_Track_dPhi.at(i)/kShortScaleHigh << endl;
  //   cout << "i=" << i << "   TrackNSigma=" << (Float_t)fBuffer_Track_NSigmaElec.at(i)/kShortScaleHigh << endl;
  
  // }
       
    // fBuffer_ClusterM02.push_back(input.ClusterM02); 
    // fBuffer_ClusterM20.push_back(input.ClusterM20); 
  
  ResetBuffer();
  //gObjectTable->Print();
  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskElectronStudies::Terminate(Option_t *){

}

//________________________________________________________________________
void AliAnalysisTaskElectronStudies::ResetBuffer(){
    fBuffer_ClusterE.clear();  
    fBuffer_ClusterM02.clear();
    fBuffer_ClusterM20.clear();
    fBuffer_ClusterNCells.clear();
    fBuffer_Track_E.clear();
    fBuffer_Track_Px.clear();
    fBuffer_Track_Py.clear();
    fBuffer_Track_Pz.clear();
    fBuffer_Track_PonEMCal.clear();
    fBuffer_Track_Charge.clear();
    fBuffer_Track_dEta.clear();
    fBuffer_Track_dPhi.clear();
    fBuffer_Track_NSigmaElec.clear();
    fBuffer_Track_IsFromV0.clear();
    fBuffer_Track_ClosestR.clear();
    fBuffer_Track_ChargedIso.clear();
    fBuffer_MatchType.clear();

    fBuffer_MC_True_Cluster_E.clear();
    fBuffer_MC_True_Track_E.clear();
    fBuffer_MC_True_Track_Px.clear();
    fBuffer_MC_True_Track_Py.clear();
    fBuffer_MC_True_Track_Pz.clear();
    fBuffer_MC_Track_Is_Electron.clear();
    fBuffer_MC_Cluster_Is_Electron.clear();
    fBuffer_MC_ClusterTrack_Same_Electron.clear();
    fBuffer_MC_True_Track_MotherPDG.clear();

    fBuffer_NClus = 0.;
    fBuffer_NPrimaryTracks = 0.;
    fBuffer_IsProblem = kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskElectronStudies::ProcessCaloPhotons(){
   Int_t nclus                         = 0;

   if(fIsMC && !fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

   TClonesArray * arrClustersProcess   = NULL;
   if(!fCorrTaskSetting.CompareTo("")){
     nclus = fInputEvent->GetNumberOfCaloClusters();
   } else {
    arrClustersProcess                = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus                        = arrClustersProcess->GetEntries();
  }
  if(nclus == 0)  return;
  fBuffer_NClus = nclus;
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
      for(Long_t i = 0; i < nclus; i++){
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
        // fTMCuts->CleanClusterLabels(clus,fAODMCTrackArray);
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

          } else{
            // check for conv matching
            if(!fReaderGammas) fReaderGammas    = fV0Reader->GetReconstructedGammas();
            AliAODTrack *inTrackClosest = 0x0;
            Float_t RClosest = 99999;
            Float_t dEtaClosest = 99999;
            Float_t dPhiClosest= 99999;
            for (Int_t conv = 0; conv < fReaderGammas->GetEntriesFast(); conv++)
            {
              AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(conv);
              for (Int_t i = 0;i < 2;i++){ // loop over daughters
                Int_t tracklabel = PhotonCandidate->GetLabel(i);
                AliAODTrack *inTrack = 0x0;
                if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled()){
                    inTrack = static_cast<AliAODTrack*>(fInputEvent->GetTrack(tracklabel));
                } else {
                    for(Int_t ii=0;ii<fInputEvent->GetNumberOfTracks();ii++) {
                      inTrack = static_cast<AliAODTrack*>(fInputEvent->GetTrack(ii));
                      if(inTrack){
                        if(inTrack->GetID() == tracklabel) {
                        break;
                      }
                    }
                  }
                } 
                Float_t dEtaConv = 0;
                Float_t dPhiConv = 0;
                Bool_t propagated = ((AliCaloTrackMatcher*) fTMCuts->GetCaloTrackMatcherInstance())->PropagateV0TrackToClusterAndGetMatchingResidual(inTrack,clus,fInputEvent,dEtaConv,dPhiConv);
                if (propagated){
                  Double_t r = TMath::Sqrt(dEtaConv*dEtaConv + dPhiConv*dPhiConv);
                  if(r < RClosest){
                     inTrackClosest = inTrack;
                     RClosest = r;
                     dEtaClosest = dEtaConv;
                     dPhiClosest = dPhiConv;
                  }
                }
              } // end daughter loop
            } // end conv photon loop 
            if((RClosest < fConversionTrackMatchR) && inTrackClosest){ 
             ProcessMatchedTrack(inTrackClosest,clus,kTRUE,dEtaClosest,dPhiClosest);
            }
          } // end conv matching case

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
           Int_t labelTrackClosest = -1;
          if(fTMCuts->GetClosestMatchedTrackToCluster(fInputEvent,clus,labelTrackClosest)){
              Int_t properLabel = labelTrackClosest;
              // if(labelTrack<0) properLabel = (-1 * labelTrack) - 1; // conversion hybrid none hybrid
              aodt = dynamic_cast<AliAODTrack*> (fInputEvent->GetTrack(properLabel));
              if(TrackIsSelectedAOD(aodt)){
                ProcessMatchedTrack(aodt,clus,kFALSE);
              }
          }
                  } else{
            // check for conv matching
            if(!fReaderGammas) fReaderGammas    = fV0Reader->GetReconstructedGammas();
            AliAODTrack *inTrackClosest = 0x0;
            Float_t RClosest = 99999;
            Float_t dEtaClosest = 99999;
            Float_t dPhiClosest= 99999;
            for (Int_t conv = 0; conv < fReaderGammas->GetEntriesFast(); conv++)
            {
              AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(conv);
              for (Int_t i = 0;i < 2;i++){ // loop over daughters
                Int_t tracklabel = PhotonCandidate->GetLabel(i);
                AliAODTrack *inTrack = 0x0;
                if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled()){
                    inTrack = static_cast<AliAODTrack*>(fInputEvent->GetTrack(tracklabel));
                } else {
                    for(Int_t ii=0;ii<fInputEvent->GetNumberOfTracks();ii++) {
                      inTrack = static_cast<AliAODTrack*>(fInputEvent->GetTrack(ii));
                      if(inTrack){
                        if(inTrack->GetID() == tracklabel) {
                        break;
                      }
                    }
                  }
                } 
                Float_t dEtaConv = 0;
                Float_t dPhiConv = 0;
                Bool_t propagated = ((AliCaloTrackMatcher*) fTMCuts->GetCaloTrackMatcherInstance())->PropagateV0TrackToClusterAndGetMatchingResidual(inTrack,clus,fInputEvent,dEtaConv,dPhiConv);
                if (propagated){
                  Double_t r = TMath::Sqrt(dEtaConv*dEtaConv + dPhiConv*dPhiConv);
                  if(r < RClosest){
                     inTrackClosest = inTrack;
                     RClosest = r;
                     dEtaClosest = dEtaConv;
                     dPhiClosest = dPhiConv;
                  }
                }
              } // end daughter loop
            } // end conv photon loop 
            if((RClosest < fConversionTrackMatchR) && inTrackClosest){ 
             ProcessMatchedTrack(inTrackClosest,clus,kTRUE,dEtaClosest,dPhiClosest);
            }
          } // end conv matching case

      //ProcessTrackMatching(clus); // tree filling done here too
      
      delete clus;
      continue;
    }
    delete clus;
  }  
}

//________________________________________________________________________
void AliAnalysisTaskElectronStudies::ProcessTracks(){
   UShort_t ntracks = 0;
   for(Int_t t=0;t<fInputEvent->GetNumberOfTracks();t++){
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*> (fInputEvent->GetTrack(t));
      if(!aodt) continue;
      if(!TrackIsSelectedAOD(aodt)) continue;
      ntracks++;
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
   fBuffer_NPrimaryTracks = ntracks;
}

//________________________________________________________________________
void AliAnalysisTaskElectronStudies::ProcessTracksESD(){
  AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);
	if(fESDEvent){
		for(Int_t i=0;i<fESDEvent->GetNumberOfTracks();i++){
            AliESDtrack *track = dynamic_cast<AliESDtrack*> (fESDEvent->GetTrack(i));
            cout << "P" << track->P() << endl;

            // Calculate P on emcal surface

          // If the esdFriend is available, use the TPCOuter point as the starting point of extrapolation
          // Otherwise use the TPCInner point. Ignored special case for ITS standalone

          Float_t fStep =  20;
          Float_t fEMCalSurfaceDistance = 440; //cm
           AliExternalTrackParam *trkParam = 0;

           const AliESDfriendTrack*  friendTrack = track->GetFriendTrack();
     
          if (friendTrack && friendTrack->GetTPCOut())
               trkParam = const_cast<AliExternalTrackParam*>(friendTrack->GetTPCOut());
          else if (track->GetInnerParam())
               trkParam = const_cast<AliExternalTrackParam*>(track->GetInnerParam());
          if(!trkParam) {
            AliWarning("Could not find track params ... Ignoring track ...");
            continue;
          }
          AliExternalTrackParam trkParamTmp(*trkParam);
          Float_t eta, phi, pt;
          if (AliEMCALRecoUtilsBase::ExtrapolateTrackToEMCalSurface(&trkParamTmp, fEMCalSurfaceDistance, track->GetMass(kTRUE), fStep, eta, phi, pt))
          {
             track->SetTrackPhiEtaPtOnEMCal(phi,eta,pt);
          } else{
            AliWarning("Could not propagate!");
            continue;
          }
          fTrackPvsPOnSurface->Fill(track->P(),track->GetTrackPOnEMCal());

          if(track->P()>=4.5 && track->P()<=5.5) fTrackPOnSurface->Fill(track->GetTrackPOnEMCal());

        
          cout << "POnEMCal Step 20 (MeV) " << track->GetTrackPOnEMCal()*1000 << endl;


           // Do own propagation now
          if (fEMCalRecoUtils->ExtrapolateTrackToEMCalSurfaceExperimental(&trkParamTmp, fEMCalSurfaceDistance, track->GetMass(kTRUE), fStep, eta, phi, pt))
          {
            cout << "Success" << endl;
             track->SetTrackPhiEtaPtOnEMCal(phi,eta,pt);
          } else{
            AliWarning("Could not propagate!");
            continue;
          }
          cout << "POnEMCal Own (MeV) " << track->GetTrackPOnEMCal()*1000 << endl;
           fTrackPvsPOnSurfaceOwn->Fill(track->P(),track->GetTrackPOnEMCal());
           if(track->P()>=4.5 && track->P()<=5.5) fTrackPOnSurfaceOwn->Fill(track->GetTrackPOnEMCal());

          
    }
  }


  // do a loop over MC particles
  AliMCEventHandler* mcinfo = (AliMCEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());  	     
  AliMCEvent* mcevent = mcinfo->MCEvent();
 
 Double_t truePOnSurface = 0;
 Double_t truePOnSurfaceBefore = 0; // last value before EMCal surface was passed
 for (Int_t ipart=0; ipart < mcevent->GetNumberOfTracks(); ipart++) {
    AliMCParticle *mcPart = (AliMCParticle*)mcevent->GetTrack(ipart);
    if (TMath::Abs(mcPart->PdgCode()) != 11) continue;
    if(mcPart->E()>5.5 || mcPart->E()<4.5) continue;
    Int_t nTrackRefs = mcPart->GetNumberOfTrackReferences();
    cout << "nTrackRefs " << nTrackRefs << endl;
    for (Int_t iTrackRef1 = 0; iTrackRef1 < nTrackRefs; iTrackRef1++){
	      AliTrackReference * trackRef1 = mcPart->GetTrackReference(iTrackRef1);
        if((iTrackRef1==0) && (trackRef1->R()>10)) break; // only consider tracks that start close to the center

            cout << trackRef1->R() << "\t" << trackRef1->P() << endl;

            fTrackRefPvsR->Fill(trackRef1->R(),trackRef1->P());
        if(trackRef1->R()>=430){
           truePOnSurface = truePOnSurfaceBefore;
            cout << "R() "<< trackRef1->R() << endl;
            cout << "P() "<< trackRef1->P() << endl;

           break;
        }
        truePOnSurfaceBefore = trackRef1->P();
    }
    if(truePOnSurface ==0) truePOnSurface =truePOnSurfaceBefore;
    cout << "------" << endl;
 }
 fTrackPOnSurfaceTrue->Fill(truePOnSurface);
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
void AliAnalysisTaskElectronStudies::ProcessMatchedTrack(AliAODTrack* track, AliAODCaloCluster* clus, Bool_t isV0, Float_t dEtaV0, Float_t dPhiV0){
    Float_t clusPos[3]                      = { 0,0,0 };
    clus->GetPosition(clusPos);
    TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
    
    treeWriteContainer output;
    
    output.ClusterE = ConvertToUShort(clus->E(),kShortScaleLow);
    output.ClusterM02 =  ConvertToUShort(clus->GetM02(),kShortScaleMiddle);
    output.ClusterNCells = (UShort_t) clus->GetNCells();
   
    // fix for rounding issue causing it to sometimes be just below 0
    Double_t m20 = clus->GetM20();
    if(m20<0) m20 =0;
    output.ClusterM20 = ConvertToUShort(m20,kShortScaleMiddle);
    output.Track_E = ConvertToUShort(track->E(),kShortScaleLow);
    output.Track_Px = ConvertToShort(track->Px(),kShortScaleLow);
    output.Track_Py = ConvertToShort(track->Py(),kShortScaleLow);
    output.Track_Pz = ConvertToShort(track->Pz(),kShortScaleLow);
    output.Track_PonEMCal = ConvertToUShort(track->GetTrackPOnEMCal(),kShortScaleLow);

    // if(clus->E() < 0) cout  << "ClusterE" << endl;
    // if(clus->GetM02() < 0) cout  << "ClusterM02" << endl;
    // if(clus->GetM20() < 0) cout  << "ClusterM20" << endl;
    // if(track->E() < 0) cout  << "TrackE" << endl;
    // if(track->GetTrackPOnEMCal() < 0) cout  << "TrackPOnEmcal" << endl;

    // Float_t LowTrack_E = Float_t(output.Track_E) / kShortScaleLow;
    // Float_t LowTrack_Px = Float_t( output.Track_Px) / kShortScaleLow;
    // Float_t LowTrack_Py = Float_t(output.Track_Py ) / kShortScaleLow;
    // Float_t LowTrack_Pz = Float_t(output.Track_Pz) / kShortScaleLow;
    // TLorentzVector vecLowResolution;
    // vecLowResolution.SetPxPyPzE(LowTrack_Px,LowTrack_Py,LowTrack_Pz,LowTrack_E);

    // if(track->Pt()>10.){
    //   cout << "-------" << endl;
    // cout << "Track pt = " << track->Pt() << " LowRes Pt = " << vecLowResolution.Pt() << "Ratio = " <<vecLowResolution.Pt()/track->Pt() << endl;
    // cout << "Track P = " << track->P() << " LowRes P = " << vecLowResolution.P() << "Ratio =" <<vecLowResolution.P()/track->P() << endl;
    // cout << "Track P = " << track->P() << " POnEmcal = " << track->GetTrackPOnEMCal() << "Ratio =" <<track->GetTrackPOnEMCal()/track->P() << endl;
    // }

    output.Track_NSigmaElec =  ConvertToShort(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron),kShortScaleMiddle); 
    output.Track_IsFromV0 = isV0; 
    output.Track_Charge = track->Charge();

    output.MC_True_Cluster_E = 0; 
    output.MC_True_Track_E = 0; 
    output.MC_True_Track_Px = 0; 
    output.MC_True_Track_Py = 0; 
    output.MC_True_Track_Pz = 0; 
    output.MC_Track_Is_Electron= kFALSE;
    output.MC_Cluster_Is_Electron= kFALSE;
    output.MC_ClusterTrack_Same_Electron= 0;
    output.MC_True_Track_MotherPDG= 0;
    output.matchType = 0;

    std::pair<Double_t,Double_t> isoStudy = ProcessChargedIsolation(track); // R, E
    output.minR = ConvertToUShort(isoStudy.first,kShortScaleLow);
    output.isoE = ConvertToUShort(isoStudy.second,kShortScaleLow);

    // if(isoStudy.first < 0) cout  << "MinR" << endl;
    // if(isoStudy.second < 0) cout  << "IsoE" << endl;
    
    Float_t tempEta = -9;
    Float_t tempPhi = -9;
    // Get dEta and dPhi of track on EMCal surface!
    if(!isV0){
      ((AliCaloTrackMatcher*) fTMCuts->GetCaloTrackMatcherInstance())->GetTrackClusterMatchingResidual(track->GetID(),clus->GetID(),tempEta,tempPhi);
      output.Track_dEta = ConvertToShort(tempEta,kShortScaleHigh); 
      output.Track_dPhi = ConvertToShort(tempPhi,kShortScaleHigh);
    } else{ // matched a V0
      output.Track_dEta = ConvertToShort(dEtaV0,kShortScaleHigh); 
      output.Track_dPhi = ConvertToShort(dPhiV0,kShortScaleHigh);
    }

    
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
                output.MC_Cluster_Is_Electron = kTRUE;
                output.MC_True_Cluster_E = ConvertToUShort(clusterMother->E(),kShortScaleLow); 
                fTruePtElectronClusterMatchedWithTrack->Fill(clusterMother->Pt(),fWeightJetJetMC);
            }
          }
        }
             
        // Check Track
        Int_t trackMCLabel = track->GetLabel();
  
        AliAODMCParticle* trackMC = NULL;
        
        if(trackMCLabel>-1){
          trackMC = (AliAODMCParticle* )fAODMCTrackArray->At(trackMCLabel);
          if(TMath::Abs(trackMC->GetPdgCode()) == 11){
            output.MC_Track_Is_Electron = kTRUE;
            output.MC_True_Track_E = ConvertToUShort(trackMC->E(),kShortScaleLow); 
            output.MC_True_Track_Px = ConvertToShort(trackMC->Px(),kShortScaleLow); 
            output.MC_True_Track_Py = ConvertToShort(trackMC->Py(),kShortScaleLow);
            output.MC_True_Track_Pz = ConvertToShort(trackMC->Pz(),kShortScaleLow);
            Int_t motherLabel = trackMC->GetMother();
            if(motherLabel>=0){
               output.MC_True_Track_MotherPDG = ((AliAODMCParticle* )fAODMCTrackArray->At(motherLabel))->PdgCode();
            }
          }
        }

        if((clus->GetNLabels()>0) && (mclabelsCluster[0] == trackMCLabel) && (output.MC_Track_Is_Electron ==1)){
          output.MC_ClusterTrack_Same_Electron = 1;
        } else if((clus->GetNLabels()>0) && (mclabelsCluster[1] == trackMCLabel) && (output.MC_Track_Is_Electron ==1)){
          output.MC_ClusterTrack_Same_Electron = 2;
        }
          
      // }
    } // end is MC

    // TODO remove high PT track match stuff, can't think of any clever way of doing this 
    // on event basis right now

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
              output.matchType = 1;
              
              PushToVectors(output);

              // Fill high pT track stuff
              output.matchType = 2;
              Float_t tempEtaHigh = -99999;
              Float_t tempPhiHigh = -99999;
              ((AliCaloTrackMatcher*) fTMCuts->GetCaloTrackMatcherInstance())->GetTrackClusterMatchingResidual(highTrack->GetID(),clus->GetID(),tempEtaHigh,tempPhiHigh);
              output.Track_dEta = ConvertToShort(tempEtaHigh,kShortScaleHigh); 
              output.Track_dPhi = ConvertToShort(tempPhiHigh,kShortScaleHigh);

              output.Track_E = ConvertToUShort(highTrack->E(),kShortScaleLow);
              output.Track_Px = ConvertToShort(highTrack->Px(),kShortScaleLow);
              output.Track_Py = ConvertToShort(highTrack->Py(),kShortScaleLow);
              output.Track_Pz = ConvertToShort(highTrack->Pz(),kShortScaleLow);
              output.Track_NSigmaElec = ConvertToShort(fPIDResponse->NumberOfSigmasTPC(highTrack,AliPID::kElectron),kShortScaleMiddle); 
              output.Track_IsFromV0 = isV0; 
              
              output.MC_True_Cluster_E = 0; 
              output.MC_True_Track_E = 0; 
              output.MC_True_Track_Px = 0; 
              output.MC_True_Track_Py = 0; 
              output.MC_True_Track_Pz = 0; 
              output.MC_True_Track_MotherPDG = 0;
              output.MC_Track_Is_Electron= kFALSE;
              output.MC_Cluster_Is_Electron= kFALSE;
              output.MC_ClusterTrack_Same_Electron= 0;

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
                          output.MC_Cluster_Is_Electron = kTRUE;
                          output.MC_True_Cluster_E = ConvertToUShort(clusterMother->E(),kShortScaleLow); 
                          fTruePtElectronClusterMatchedWithTrack->Fill(clusterMother->Pt(),fWeightJetJetMC);
                      }
                    }
                  }
            
                      
                  // Check Track
                  Int_t trackMCLabel = track->GetLabel();
                  AliAODMCParticle* trackMC = NULL;
                  if(trackMCLabel>-1){
                    trackMC = (AliAODMCParticle* )fAODMCTrackArray->At(trackMCLabel);
                    if(TMath::Abs(trackMC->GetPdgCode()) == 11){
                      output.MC_Track_Is_Electron = kTRUE;
                      output.MC_True_Track_E = ConvertToUShort(trackMC->E(),kShortScaleLow); 
                      output.MC_True_Track_Px = ConvertToShort(trackMC->Px(),kShortScaleLow); 
                      output.MC_True_Track_Py = ConvertToShort(trackMC->Py(),kShortScaleLow);
                      output.MC_True_Track_Pz = ConvertToShort(trackMC->Pz(),kShortScaleLow);
                      Int_t motherLabel = trackMC->GetMother();
                      if(motherLabel>=0){
                        output.MC_True_Track_MotherPDG = ((AliAODMCParticle* )fAODMCTrackArray->At(motherLabel))->PdgCode();
                      }
                    }
                  }

                  if((clus->GetNLabels()>0) && (mclabelsCluster[0] == trackMCLabel) && (output.MC_Track_Is_Electron ==1)){
                    output.MC_ClusterTrack_Same_Electron = 1;
                  } else if((clus->GetNLabels()>0) && (mclabelsCluster[1] == trackMCLabel) && (output.MC_Track_Is_Electron ==1)){
                    output.MC_ClusterTrack_Same_Electron = 2;
                  }
                    
                // }
              } // end is MC
              PushToVectors(output);
            } else{
              output.matchType = 0;
              PushToVectors(output);
            } 
          }
    } else{
      // did not work or did not find
      output.matchType = 0;
      PushToVectors(output);
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
// Charged isolation for tracks 
//_____________________________________________________________________________
std::pair<Double_t,Double_t> AliAnalysisTaskElectronStudies::ProcessChargedIsolation(AliAODTrack* track){
    TLorentzVector v4ThisTrack;
    v4ThisTrack.SetPxPyPzE(track->Px(),track->Py(),track->Pz(),track->E());

    Double_t rMin =9999;
    Double_t EIso = 0;
    for (Int_t t = 0; t < fInputEvent->GetNumberOfTracks(); t++)
    {
        AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
        if(!aodt) continue;
        if(!TrackIsSelectedAOD(aodt)) continue;
        if(IsSameTrack(track->GetID(),aodt->GetID())) continue;

        TLorentzVector v4OtherTrack;
        v4OtherTrack.SetPxPyPzE(aodt->Px(),aodt->Py(),aodt->Pz(),aodt->E());
        
        Double_t thistrackEta = v4ThisTrack.Eta();
        Double_t thistrackPhi = v4ThisTrack.Phi();
        if (thistrackPhi < 0) thistrackPhi += 2*TMath::Pi();

        Double_t othertrackEta = v4OtherTrack.Eta();
        Double_t othertrackPhi = v4OtherTrack.Phi();
        if (othertrackPhi < 0) othertrackPhi += 2*TMath::Pi();

        Double_t dEta = thistrackEta - othertrackEta;
        Double_t dPhi = thistrackPhi - othertrackPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

        if(dR <= rMin) rMin = dR; // set as new minimum
        
        //if((dR > fTrackIsolationR[0]) && (dR > fTrackIsolationR[1])) continue;

        // track is in cone
        // check if track comes from pi0
        if(dR<= fIsoMaxRadius) EIso += aodt->E();
    }
    return std::make_pair(rMin,EIso);
}
 Short_t AliAnalysisTaskElectronStudies::ConvertToShort(Float_t input, Int_t scale){
    Short_t value = Short_t(input*scale);
    Int_t min = std::numeric_limits<short>::min();
    Int_t max = std::numeric_limits<short>::max();
    if((input*scale) > max){
      fBuffer_IsProblem = kTRUE;
      value = Short_t(max);
    } else if((input*scale) < min){
      fBuffer_IsProblem = kTRUE;
      value = Short_t(min);
    }
    return value;
 }
  Short_t AliAnalysisTaskElectronStudies::ConvertToShort(Double_t input, Int_t scale){
    Short_t value = Short_t(input*scale);
    Int_t min = std::numeric_limits<short>::min();
    Int_t max = std::numeric_limits<short>::max();
    if((input*scale) > max){
      fBuffer_IsProblem = kTRUE;
      value = Short_t(max);
    } else if((input*scale) < min){
      fBuffer_IsProblem = kTRUE;
      value = Short_t(min);
    }
    return value;
 }
UShort_t AliAnalysisTaskElectronStudies::ConvertToUShort(Float_t input, Int_t scale){
    UShort_t value = UShort_t(input*scale);
    Int_t min = std::numeric_limits<ushort>::min();
    Int_t max = std::numeric_limits<ushort>::max();
    if((input*scale) > max){
      fBuffer_IsProblem = kTRUE;
      value = UShort_t(max);
    } else if((input*scale) < min){
      fBuffer_IsProblem = kTRUE;
      value = UShort_t(min);
    }
    return value;
 }
 UShort_t AliAnalysisTaskElectronStudies::ConvertToUShort(Double_t input, Int_t scale){
    UShort_t value = UShort_t(input*scale);
    Int_t min = std::numeric_limits<ushort>::min();
    Int_t max = std::numeric_limits<ushort>::max();
    if((input*scale) > max){
      fBuffer_IsProblem = kTRUE;
      value = UShort_t(max);
    } else if((input*scale) < min){
      fBuffer_IsProblem = kTRUE;
      value = UShort_t(min);
    }
    return value;
 }
void AliAnalysisTaskElectronStudies::PushToVectors(treeWriteContainer input){
    // Push back values to vectors and set buffers
    fBuffer_ClusterE.push_back(input.ClusterE);     
    fBuffer_ClusterM02.push_back(input.ClusterM02); 
    fBuffer_ClusterM20.push_back(input.ClusterM20); 
    fBuffer_ClusterNCells.push_back(input.ClusterNCells); 
    fBuffer_Track_E.push_back(input.Track_E);
    fBuffer_Track_Px.push_back(input.Track_Px);
    fBuffer_Track_Py.push_back(input.Track_Py); 
    fBuffer_Track_Pz.push_back(input.Track_Pz); 
    fBuffer_Track_PonEMCal.push_back(input.Track_PonEMCal); 
    fBuffer_Track_ClosestR.push_back(input.minR); 
    fBuffer_Track_ChargedIso.push_back(input.isoE); 
    fBuffer_MatchType.push_back(input.matchType); 
    fBuffer_Track_Charge.push_back(input.Track_Charge);
    fBuffer_Track_dEta.push_back(input.Track_dEta); 
    fBuffer_Track_dPhi.push_back(input.Track_dPhi); 
    fBuffer_Track_NSigmaElec.push_back(input.Track_NSigmaElec); 
    fBuffer_Track_IsFromV0.push_back(input.Track_IsFromV0); 

    fBuffer_MC_True_Cluster_E.push_back(input.MC_True_Cluster_E); 
    fBuffer_MC_True_Track_E.push_back(input.MC_True_Track_E);
    fBuffer_MC_True_Track_Px.push_back(input.MC_True_Track_Px); 
    fBuffer_MC_True_Track_Py.push_back(input.MC_True_Track_Py); 
    fBuffer_MC_True_Track_Pz.push_back(input.MC_True_Track_Pz); 
    fBuffer_MC_Track_Is_Electron.push_back(input.MC_Track_Is_Electron); 
    fBuffer_MC_Cluster_Is_Electron.push_back(input.MC_Cluster_Is_Electron); 
    fBuffer_MC_ClusterTrack_Same_Electron.push_back(input.MC_ClusterTrack_Same_Electron);
    fBuffer_MC_True_Track_MotherPDG.push_back(input.MC_True_Track_MotherPDG);
 }







