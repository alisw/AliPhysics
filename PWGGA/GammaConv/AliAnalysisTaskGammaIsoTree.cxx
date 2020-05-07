/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Nicolas Schmidt                                               *
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

#include "AliAnalysisTaskGammaIsoTree.h"
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

ClassImp(AliAnalysisTaskGammaIsoTree)
//________________________________________________________________________
AliAnalysisTaskGammaIsoTree::AliAnalysisTaskGammaIsoTree() : AliAnalysisTaskSE(),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fWeightJetJetMC(1),
  fOutputList(NULL),
  fAnalysisTree(NULL),
  fIsMC(0),
  fIsHeavyIon(0),
  fV0Reader(NULL),
  fV0ReaderName(""),
  fReaderGammas(NULL),
  fConversionCandidates(0),
  fClusterEMCalCandidates(0),
  fClusterEMCalCandidatesNoCuts(0),
  fClusterPHOSCandidates(0),
  fTracks(0),
  fMCParticles(0),
  fIDMatchedEMCTrack(0),
  fIDMatchedEMCTrackNoCuts(0),
  fDataEvtHeader(),
  fMCEvtHeader(),
  fConvIsoInfo(0),
  fCaloIsoInfo(0),
  //fConvChargedTrackISORaw(0),
  //fEMCalClusterChargedTrackISORaw(0),
  fGeomEMCAL(NULL),
  fCorrTaskSetting(""),
  fEventCuts(NULL),
  fClusterCutsEMC(NULL),
  fClusterCutsPHOS(NULL),
  fConvCuts(NULL),
  fMinClsTPC(0),
  fChi2PerClsTPC(9999),
  fMinClsITS(0),
  fEtaCut(9999),
  fPtCut(0),
  fYMCCut(9999),
  fMatchingParamsPhi(),
  fMatchingParamsEta(),
  fMatchingEOverP(),
  fDoTrackIsolation(kFALSE),
  fTrackIsolationR(0.4),
  fDoNeutralIsolation(kFALSE),
  fNeutralIsolationR(0.4),
  fPi0TaggingWindow(),
  fEtaTaggingWindow(),
  fSaveConversions(kTRUE),
  fSaveEMCClusters(kTRUE),
  fSavePHOSClusters(kTRUE),
  fSaveTracks(kTRUE),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fRhoOutName("Rho"),
  fTreeBuffSize(100000),
  fMemCountAOD(0)
{
  fDataEvtHeader.pVtxX = -9999;
  fDataEvtHeader.pVtxY = -9999;
  fDataEvtHeader.pVtxZ = -9999;
  fDataEvtHeader.runnumber = -1;
  fDataEvtHeader.numberESDtracks = -1;
  
  fMCEvtHeader.pVtxX = -9999;
  fMCEvtHeader.pVtxY = -9999;
  fMCEvtHeader.pVtxZ = -9999;
  fMCEvtHeader.runnumber = -1;
  fMCEvtHeader.numberESDtracks = -1;

  SetEtaMatching(0.010,4.07,-2.5);
  SetPhiMatching(0.015,3.65,3.65);
  SetEOverP(1.75);

  SetPi0TaggingWindow(0.120,0.145);
  SetEtaTaggingWindow(0.5,0.6);

}

AliAnalysisTaskGammaIsoTree::AliAnalysisTaskGammaIsoTree(const char *name) : AliAnalysisTaskSE(name),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fWeightJetJetMC(1),
  fOutputList(NULL),
  fAnalysisTree(NULL),
  fIsMC(0),
  fIsHeavyIon(0),
  fV0Reader(NULL),
  fV0ReaderName(""),
  fReaderGammas(NULL),
  fConversionCandidates(0),
  fClusterEMCalCandidates(0),
  fClusterEMCalCandidatesNoCuts(0),
  fClusterPHOSCandidates(0),
  fTracks(0),
  fMCParticles(0),
  fIDMatchedEMCTrack(0),
  fIDMatchedEMCTrackNoCuts(0),
  fDataEvtHeader(),
  fMCEvtHeader(),
  fConvIsoInfo(0),
  fCaloIsoInfo(0),
  //fConvChargedTrackISORaw(0),
  //fEMCalClusterChargedTrackISORaw(0),
  fGeomEMCAL(NULL),
  fCorrTaskSetting(""),
  fEventCuts(NULL),
  fClusterCutsEMC(NULL),
  fClusterCutsPHOS(NULL),
  fConvCuts(NULL),
  fMinClsTPC(0),
  fChi2PerClsTPC(9999),
  fMinClsITS(0),
  fEtaCut(9999),
  fPtCut(0),
  fYMCCut(99999),
  fMatchingParamsPhi(),
  fMatchingParamsEta(),
  fMatchingEOverP(),
  fDoTrackIsolation(kFALSE),
  fTrackIsolationR(0.4),
  fDoNeutralIsolation(kFALSE),
  fNeutralIsolationR(0.4),
  fPi0TaggingWindow(),
  fEtaTaggingWindow(),
  fSaveConversions(kTRUE),
  fSaveEMCClusters(kTRUE),
  fSavePHOSClusters(kTRUE),
  fSaveTracks(kTRUE),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fRhoOutName("Rho"),
  fTreeBuffSize(100000),
  fMemCountAOD(0)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  
  fDataEvtHeader.pVtxX = -9999;
  fDataEvtHeader.pVtxY = -9999;
  fDataEvtHeader.pVtxZ = -9999;
  fDataEvtHeader.runnumber = -1;
  fDataEvtHeader.numberESDtracks = -1;
  
  fMCEvtHeader.pVtxX = -9999;
  fMCEvtHeader.pVtxY = -9999;
  fMCEvtHeader.pVtxZ = -9999;
  fMCEvtHeader.runnumber = -1;
  fMCEvtHeader.numberESDtracks = -1;

  SetEtaMatching(0.010,4.07,-2.5);
  SetPhiMatching(0.015,3.65,3.65);
  SetEOverP(1.75);

  SetPi0TaggingWindow(0.120,0.145);
  SetEtaTaggingWindow(0.5,0.6);
}

//________________________________________________________________________
AliAnalysisTaskGammaIsoTree::~AliAnalysisTaskGammaIsoTree()
{
  // default deconstructor

}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::UserCreateOutputObjects()
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

  if(((AliCaloPhotonCuts*)fClusterCutsPHOS)->GetCutHistograms()){
    fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsPHOS)->GetCutHistograms());
  }

  if(((AliConversionPhotonCuts*)fConvCuts)->GetCutHistograms()){
    fOutputList->Add(((AliConversionPhotonCuts*)fConvCuts)->GetCutHistograms());
  }


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
  fHistoNEvents->GetXaxis()->SetBinLabel(10,"EMCAL problem");
  fHistoNEvents->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
  fHistoNEvents->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
  fHistoNEvents->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
  fHistoNEvents->GetYaxis()->SetTitle("N_{events}");
  fOutputList->Add(fHistoNEvents);

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
  PostData(1, fOutputList);

  fAnalysisTree = new TTree("AnalysisTree","AnalysisTree");
  fAnalysisTree->Branch("fDataEvtHeader",&fDataEvtHeader,"pVtxX/D:pVtxY/D:pVtxZ/D:runnumber/I:numberESDtracks/I:rho/D");
  fAnalysisTree->Branch("fConversionCandidates","vector<AliAODConversionPhoton*>",&fConversionCandidates,101);
  fAnalysisTree->Branch("fClusterEMCalCandidates","vector<AliAODCaloCluster*>",&fClusterEMCalCandidates,101);
  fAnalysisTree->Branch("fClusterPHOSCandidates","vector<AliAODCaloCluster*>",&fClusterPHOSCandidates,101);
  fAnalysisTree->Branch("fTracks","vector<AliAODTrack*>",&fTracks,101);
  fAnalysisTree->Branch("fIDMatchedEMCTrack","vector<Short_t>",&fIDMatchedEMCTrack,101);
  fAnalysisTree->Branch("fConvIsoInfo","vector<isoInfo>",&fConvIsoInfo,101);
  fAnalysisTree->Branch("fCaloIsoInfo","vector<isoInfo>",&fCaloIsoInfo,101);
  if(fIsMC>0){
    fAnalysisTree->Branch("fMCParticles","vector<AliAODMCParticle*>",&fMCParticles,101);
    fAnalysisTree->Branch("fMCEvtHeader",&fMCEvtHeader,"pVtxX/D:pVtxY/D:pVtxZ/D:runnumber/I:numberESDtracks/I:weightJJ/F:rho/D");
  }

  OpenFile(2);
  PostData(2,fAnalysisTree);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaIsoTree::Notify()
{
    return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::UserExec(Option_t *){

  fInputEvent                         = InputEvent();
  if(fIsMC>0) fMCEvent                  = MCEvent();
  // Get V0 reader
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
  
  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(InputEvent()->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    fHistoNEvents->Fill(eventQuality);
    if (fIsMC>1) fHistoNEventsWOWeight->Fill(eventQuality);
    return;
  }
  Int_t eventNotAccepted              = fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
  if(eventNotAccepted) return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1



  AliRhoParameter* outrho= (AliRhoParameter*) InputEvent()->FindListObject(fRhoOutName.Data());
  if(!outrho) AliFatal("could not find rho container!");
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


  //
  // ─── MAIN PROCESSING ────────────────────────────────────────────────────────────
  //

  ProcessTracks(); // always run ProcessTracks before calo photons! (even if save tracks is false)
  ProcessCaloPhotons(); // track matching is done here as well
  if(fSaveConversions)
    ProcessConversionPhotons();
  ReduceTrackInfo(); // track matching is done, we can remove cov matrix etc now
  if(fIsMC>0) ProcessMCParticles();
  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  fDataEvtHeader.pVtxX = vertex[0];
  fDataEvtHeader.pVtxY = vertex[1];
  fDataEvtHeader.pVtxZ = vertex[2];
  fDataEvtHeader.rho = outrho->GetVal();
  fDataEvtHeader.runnumber = InputEvent()->GetRunNumber();
  fDataEvtHeader.numberESDtracks = InputEvent()->GetNumberOfESDTracks();

  if(fIsMC>0){
    fMCEvtHeader.pVtxX = fMCEvent->GetPrimaryVertex()->GetX();
    fMCEvtHeader.pVtxY = fMCEvent->GetPrimaryVertex()->GetY();
    fMCEvtHeader.pVtxZ = fMCEvent->GetPrimaryVertex()->GetZ();
    fMCEvtHeader.runnumber = fMCEvent->GetRunNumber();
    fMCEvtHeader.numberESDtracks = fMCEvent->GetNumberOfESDTracks();
    fMCEvtHeader.weightJJ = fWeightJetJetMC;
    fMCEvtHeader.rho = outrho->GetVal();
  }

  // do not fill these in tree if user wants it
  // processing was needed anyways because of track matching
  // and isolation
  if(!fSaveTracks) fTracks.clear();
  if(!fSaveEMCClusters) fClusterEMCalCandidates.clear();
  if(!fSavePHOSClusters) fClusterPHOSCandidates.clear();

  // fill output
  fAnalysisTree->Fill();
  ResetBuffer();
  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::Terminate(Option_t *){

}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ResetBuffer(){
  // for vectors of pointers, memory needs to be freed manually
  for (auto p : fConversionCandidates){ delete p;} 
  fConversionCandidates.clear();
  for (auto p : fClusterEMCalCandidates){ delete p;} 
  fClusterEMCalCandidates.clear();
  for (auto p : fClusterEMCalCandidatesNoCuts){ delete p;} 
  fClusterEMCalCandidatesNoCuts.clear();
  for (auto p : fClusterPHOSCandidates){ delete p;} 
  fClusterPHOSCandidates.clear();
  for (auto p : fTracks){ delete p;} 
  fTracks.clear();
  for (auto p : fMCParticles){ delete p;} 
  fMCParticles.clear();
  fIDMatchedEMCTrack.clear();
  fIDMatchedEMCTrackNoCuts.clear();
  fConvIsoInfo.clear();
  fCaloIsoInfo.clear();

  fConversionCandidates.resize(0);
  fClusterEMCalCandidates.resize(0);
  fClusterPHOSCandidates.resize(0);
  fTracks.resize(0);
  fMCParticles.resize(0);
  fIDMatchedEMCTrack.resize(0);
  fConvIsoInfo.resize(0);
  fCaloIsoInfo.resize(0);

  fDataEvtHeader.pVtxX = -9999;
  fDataEvtHeader.pVtxY = -9999;
  fDataEvtHeader.pVtxZ = -9999;
  fDataEvtHeader.runnumber = -1;
  fDataEvtHeader.numberESDtracks = -1;
  
  fMCEvtHeader.pVtxX = -9999;
  fMCEvtHeader.pVtxY = -9999;
  fMCEvtHeader.pVtxZ = -9999;
  fMCEvtHeader.runnumber = -1;
  fMCEvtHeader.numberESDtracks = -1;
}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessConversionPhotons(){
   fReaderGammas    = fV0Reader->GetReconstructedGammas();
   for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    if(!((AliConversionPhotonCuts*)fConvCuts)->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    fConversionCandidates.push_back(new AliAODConversionPhoton(PhotonCandidate));
    
    if(fDoTrackIsolation){
      isoInfo convIso;
      convIso.isoRawCharged = ProcessChargedIsolation(PhotonCandidate);
      convIso.isoRawNeutral = ProcessNeutralIsolation(PhotonCandidate);
      convIso.isTagged = ProcessTagging(PhotonCandidate); 
      fConvIsoInfo.push_back(convIso);
    }
    
    if(fMCEvent){
    }

   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessCaloPhotons(){
   Int_t nclus                         = 0;
   Int_t nclusCorr                     = 0;
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

        Int_t matchIndex = ProcessTrackMatching(clus,fTracks);
        fIDMatchedEMCTrackNoCuts.push_back(matchIndex);
        fClusterEMCalCandidatesNoCuts.push_back(new AliAODCaloCluster(*clus));

        // check if given EMC cuts are fulfilled
        if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
          delete clus;
          continue;
        }
        
        fIDMatchedEMCTrack.push_back(matchIndex);
        fClusterEMCalCandidates.push_back(new AliAODCaloCluster(*clus));
        if(fDoTrackIsolation){
          isoInfo caloIso;
          caloIso.isoRawCharged = ProcessChargedIsolation(clus);
          caloIso.isoRawNeutral = ProcessNeutralIsolation(clus);
          caloIso.isTagged = ProcessTagging(clus);  // TODO
          fCaloIsoInfo.push_back(caloIso);
        }
        delete clus;
      }
  }

  // Loop over normal clusters
  for(Long_t i = 0; i < nclus; i++){
    Double_t tempClusterWeight              =  fWeightJetJetMC;                   
    clus                                = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    
    if(!clus) continue;

    if(!arrClustersProcess && clus->IsEMCAL()){ // if is was not saved already
      // check if given EMC cuts are fulfilled
      Int_t matchIndex = ProcessTrackMatching(clus,fTracks);
      fIDMatchedEMCTrackNoCuts.push_back(matchIndex);
      fClusterEMCalCandidatesNoCuts.push_back(new AliAODCaloCluster(*clus));
      if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        delete clus;
        continue;
      }
      fIDMatchedEMCTrack.push_back(matchIndex);
      fClusterEMCalCandidates.push_back(new AliAODCaloCluster(*clus));
      if(fDoTrackIsolation){
          isoInfo caloIso;
          caloIso.isoRawCharged = ProcessChargedIsolation(clus);
          caloIso.isoRawNeutral = ProcessNeutralIsolation(clus);
          caloIso.isTagged = ProcessTagging(clus);  // TODO
          fCaloIsoInfo.push_back(caloIso);
      }
      delete clus;
      continue;
    }
    if(clus->IsPHOS()){
    // if(clus->GetType() == AliVCluster::kPHOSNeutral){
      if(!((AliCaloPhotonCuts*)fClusterCutsPHOS)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        delete clus;
        continue;
      }

      fClusterPHOSCandidates.push_back(new AliAODCaloCluster(*clus));
      delete clus;
      continue;
    }
    delete clus;
  }
}

///________________________________________________________________________
Bool_t AliAnalysisTaskGammaIsoTree::TrackIsSelectedAOD(AliAODTrack* lTrack) {
  // apply filter bits 
  if( ! lTrack->IsHybridGlobalConstrainedGlobal()){
    return kFALSE;
  }

	// Absolute TPC Cluster cut
	if(lTrack->GetTPCNcls()<fMinClsTPC) return kFALSE;
	if(lTrack->GetTPCchi2perCluster()>fChi2PerClsTPC) return kFALSE;
  // DCA cut 
  //if(!IsDCACutAccepted(lTrack)) return kFALSE;

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
  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessTracks(){
  for(Int_t t=0;t<fInputEvent->GetNumberOfTracks();t++){
      AliAODTrack *fCurrentTrack = dynamic_cast<AliAODTrack*> (fInputEvent->GetTrack(t));
      if(!TrackIsSelectedAOD(fCurrentTrack)) continue;

      fTracks.push_back(new AliAODTrack(*fCurrentTrack));
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessMCParticles(){
  // Loop over all primary MC particle
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  if (AODMCTrackArray){
    for(Int_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
      AliAODMCParticle* particle =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
      if(TMath::Abs(particle->Y())< fYMCCut){
        fMCParticles.push_back(new AliAODMCParticle(*particle));
      } else
      {
        fMCParticles.push_back(new AliAODMCParticle());
      }
    }
  }
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaIsoTree::ProcessTrackMatching(AliAODCaloCluster* clus, std::vector<AliAODTrack*> tracks){
     Int_t nModules = fGeomEMCAL->GetNumberOfSuperModules();
     Int_t highestMatchIndex = -1;
     for (UInt_t t = 0; t < tracks.size(); t++)
     {
      AliExternalTrackParam *trackParam = 0;
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(tracks.at(t));
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
        if(!AliEMCALRecoUtils::ExtrapolateTrackToCluster(&trackParamTmp, clus, 0.139, 5., dEta, dPhi)) continue;

        if(dEta > (fMatchingParamsEta[0] + pow(aodt->Pt() + fMatchingParamsEta[1],fMatchingParamsEta[2]))) continue;
        if(dPhi > (fMatchingParamsPhi[0] + pow(aodt->Pt() + fMatchingParamsPhi[1],fMatchingParamsPhi[2]))) continue;
        if((clus->E()/aodt->P()) > fMatchingEOverP) continue;
        if(highestMatchIndex == -1){ // this is the first match
           highestMatchIndex = t;
        } else{
           if(aodt->P()>tracks.at(highestMatchIndex)->P()){
             highestMatchIndex = t;
           }
        }

        delete trackParam;

     } 
     return highestMatchIndex; 
}

// Charged isolation for conversion photons
//_____________________________________________________________________________
Double32_t AliAnalysisTaskGammaIsoTree::ProcessChargedIsolation(AliAODConversionPhoton* photon){
    TLorentzVector* v4photon = new TLorentzVector();
    Double32_t isoRaw = 0.;
    

    v4photon->SetPxPyPzE(photon->Px(),photon->Py(),photon->Pz(),photon->E());
    for (Int_t t = 0; t < fInputEvent->GetNumberOfTracks(); t++)
    {
        AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
        if(!aodt) continue;
        if(!TrackIsSelectedAOD(aodt)) continue;
        TLorentzVector v4track;
        v4track.SetPxPyPzE(aodt->Px(),aodt->Py(),aodt->Pz(),aodt->E());
        Float_t dR = v4photon->DeltaR(v4track);
        
        if(dR > fTrackIsolationR) continue;

        // track is in cone
        // check if track comes from pi0
         Bool_t trackIsFromV0 = kFALSE;
        // check that track is not from conversion
        for (Int_t iElec = 0;iElec < 2;iElec++){
          Int_t tracklabel = photon->GetLabel(iElec);
          if(tracklabel==t){
            trackIsFromV0 = kTRUE;
          } else {
            trackIsFromV0 = kFALSE;
          }
        }
        if(trackIsFromV0) continue;

        isoRaw += v4track.Pt();
    }
    delete v4photon;
    return isoRaw;
}

// charged isolation for clusters
//_____________________________________________________________________________
Double32_t AliAnalysisTaskGammaIsoTree::ProcessChargedIsolation(AliAODCaloCluster* cluster){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    
    TLorentzVector v4cluster;
    cluster->GetMomentum(v4cluster,vertex);
    Double32_t isoRaw = 0.;
    
    for (Int_t t = 0; t < fInputEvent->GetNumberOfTracks(); t++)
    {
        AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
        if(!aodt) continue;
        if(!TrackIsSelectedAOD(aodt)) continue;
        TLorentzVector v4track;
        v4track.SetPxPyPzE(aodt->Px(),aodt->Py(),aodt->Pz(),aodt->E());
        Float_t dR = v4cluster.DeltaR(v4track);
        
        if(dR > fTrackIsolationR) continue;
        isoRaw += v4track.Pt();
    }
    return isoRaw;
}

//_____________________________________________________________________________
Double32_t AliAnalysisTaskGammaIsoTree::ProcessNeutralIsolation(AliAODConversionPhoton* photon){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
   
    TLorentzVector* v4photon = new TLorentzVector();
    Double32_t isoRaw = 0.;
  
    v4photon->SetPxPyPzE(photon->Px(),photon->Py(),photon->Pz(),photon->E());
    for (UInt_t c = 0; c < fClusterEMCalCandidatesNoCuts.size(); c++)
    {
        AliAODCaloCluster* clusterE = (AliAODCaloCluster*) fClusterEMCalCandidatesNoCuts.at(c);
        if(!clusterE) continue;
        // check if cluster is neutral
        if(fIDMatchedEMCTrackNoCuts.at(c) != -1) continue;

    
        TLorentzVector v4cluster;
        clusterE->GetMomentum(v4cluster,vertex);
        
        Double_t dR = v4photon->DeltaR(v4cluster);
        if(dR > fNeutralIsolationR) continue;

        isoRaw += clusterE->E();
    }
    delete v4photon;
    return isoRaw;
}
//_____________________________________________________________________________
Double32_t AliAnalysisTaskGammaIsoTree::ProcessNeutralIsolation(AliAODCaloCluster* cluster){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    Double32_t isoRaw = 0.;

    TLorentzVector tmp;
    cluster->GetMomentum(tmp,vertex);
    
    TLorentzVector* v4thiscluster = new TLorentzVector(tmp);
    for (UInt_t c = 0; c < fClusterEMCalCandidatesNoCuts.size(); c++)
    {
        AliAODCaloCluster* clusterE = (AliAODCaloCluster*) fClusterEMCalCandidatesNoCuts.at(c);
        if(!clusterE) continue;
        if(clusterE->GetID() == cluster->GetID()) continue;
        // check if cluster is neutral
        if(fIDMatchedEMCTrackNoCuts.at(c) != -1) continue;

        TLorentzVector v4othercluster;
        clusterE->GetMomentum(v4othercluster,vertex);
        
        Double_t dR = v4thiscluster->DeltaR(v4othercluster);
        if(dR > fNeutralIsolationR) continue;

        isoRaw += v4othercluster.Et();
    }
    delete v4thiscluster;
    return isoRaw;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaIsoTree::ProcessTagging(AliAODConversionPhoton* photon){
  Int_t taggedConv = 0;
  Int_t taggedClus = 0;
  AliAODConversionMother *pi0cand = NULL;
  fReaderGammas    = fV0Reader->GetReconstructedGammas();
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* otherPhoton = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!otherPhoton) continue; // bla
    if(!((AliConversionPhotonCuts*)fConvCuts)->PhotonIsSelected(otherPhoton,fInputEvent)) continue;
    
    pi0cand = new AliAODConversionMother(photon,otherPhoton);
    
    // check mass window
    Double_t mass = pi0cand->M();
    if((mass > fPi0TaggingWindow[0]) && (mass < fPi0TaggingWindow[1])){ // pi0
        taggedConv = 1;
        delete pi0cand;
        break;
    } else if((mass > fEtaTaggingWindow[0]) && (mass < fEtaTaggingWindow[1])){ // eta
        taggedConv = 1;
        delete pi0cand;
        break;
    }   
    delete pi0cand; 
  }

  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  for (UInt_t c = 0; c < fClusterEMCalCandidates.size(); c++)
  {
    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    fClusterEMCalCandidatesNoCuts.at(c)->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCluster = new AliAODConversionPhoton(tmpvec);
    if(!PhotonCluster){ delete tmpvec; continue;}

    pi0cand = new AliAODConversionMother(photon,PhotonCluster);
    
    // check mass window
    Double_t mass = pi0cand->M();
    if((mass > fPi0TaggingWindow[0]) && (mass < fPi0TaggingWindow[1])){ // pi0
        taggedClus = 2;
        delete pi0cand;
        delete PhotonCluster; 
        delete tmpvec;
        break;
    } else  if((mass > fEtaTaggingWindow[0]) && (mass < fEtaTaggingWindow[1])){ // eta
        delete pi0cand;
        delete PhotonCluster;
        delete tmpvec; 
        taggedClus = 2;
        break;
    }  
    delete pi0cand;
    delete PhotonCluster;
    delete tmpvec; 
  }
  return taggedConv+taggedClus;
}

// tagging of calo clusters
Int_t AliAnalysisTaskGammaIsoTree::ProcessTagging(AliAODCaloCluster* cluster){
  Int_t taggedConv = 0;
  Int_t taggedClus = 0;
  
  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  TLorentzVector tmp;
  cluster->GetMomentum(tmp,vertex);

  TLorentzVector* thisclustervec = new TLorentzVector();
  thisclustervec->SetPxPyPzE(tmp.Px(),tmp.Py(),tmp.Pz(),tmp.E());

  // convert to AODConversionPhoton
  AliAODConversionPhoton *ThisPhotonCluster = new AliAODConversionPhoton(thisclustervec);
  if(!ThisPhotonCluster){ delete thisclustervec; return -1;}
  
  AliAODConversionMother *pi0cand = NULL;
  fReaderGammas    = fV0Reader->GetReconstructedGammas();
  
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* otherPhoton = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!otherPhoton) continue; // bla
    if(!((AliConversionPhotonCuts*)fConvCuts)->PhotonIsSelected(otherPhoton,fInputEvent)) continue;
    
    pi0cand = new AliAODConversionMother(ThisPhotonCluster,otherPhoton);
    
    // check mass window
    Double_t mass = pi0cand->M();
    if((mass > fPi0TaggingWindow[0]) && (mass < fPi0TaggingWindow[1])){ // pi0
        taggedConv = 1;
        delete pi0cand;
        break;
    } else if((mass > fEtaTaggingWindow[0]) && (mass < fEtaTaggingWindow[1])){ // eta
        taggedConv = 1;
        delete pi0cand;
        break;
    }  
    delete pi0cand; 
  }


  for (UInt_t c = 0; c < fClusterEMCalCandidates.size(); c++)
  {
    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    
    if(fClusterEMCalCandidatesNoCuts.at(c)->GetID() == cluster->GetID()) continue;
    fClusterEMCalCandidatesNoCuts.at(c)->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCluster = new AliAODConversionPhoton(tmpvec);
    if(!PhotonCluster){ delete tmpvec; continue;}

    pi0cand = new AliAODConversionMother(ThisPhotonCluster,PhotonCluster);
    
    // check mass window
    Double_t mass = pi0cand->M();
    if((mass > fPi0TaggingWindow[0]) && (mass < fPi0TaggingWindow[1])){ // pi0
        taggedClus = 2;
        delete pi0cand;
        delete PhotonCluster; 
        delete tmpvec;
        break;
    } else  if((mass > fEtaTaggingWindow[0]) && (mass < fEtaTaggingWindow[1])){ // eta
        delete pi0cand;
        delete PhotonCluster;
        delete tmpvec; 
        taggedClus = 2;
        break;
    }  
    delete pi0cand;
    delete PhotonCluster;
    delete tmpvec; 
  }
  delete ThisPhotonCluster;
  delete thisclustervec;
  return taggedConv+taggedClus;
}

// delete covariance matrix etc. to save space
void AliAnalysisTaskGammaIsoTree::ReduceTrackInfo(){
    for (UInt_t i = 0; i < fTracks.size(); i++)
    {
      // fTracks.at(i)->SetTOFchi2(0);
      // fTracks.at(i)->SetTOFLabel(NULL);
      // fTracks.at(i)->SetTOFsignalDx(0);
      // fTracks.at(i)->SetTOFsignalDz(0);
      fTracks.at(i)->Clear();
    }
}
