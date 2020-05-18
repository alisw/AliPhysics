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
  fClusterEMCalCandidatesBackground(0),
  fClusterPHOSCandidates(0),
  fTracks(0),
  fMCParticles(0),
  fExtraClusterInfo(),
  fExtraClusterInfoBackground(),
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
  fClusterCutsBackgroundEMC(NULL),
  fClusterCutsPHOS(NULL),
  fConvCuts(NULL),
  fCaloUtils(NULL),
  fMinClsTPC(0),
  fChi2PerClsTPC(9999),
  fMinClsITS(0),
  fEtaCut(9999),
  fPtCut(0),
  fYMCCut(9999),
  fMatchingParamsPhi(),
  fMatchingParamsEta(),
  fMatchingEOverP(),
  fDoBackgroundTrackMatching(kFALSE),
  fDoTrackIsolation(kFALSE),
  fTrackIsolationR(),
  fDoNeutralIsolation(kFALSE),
  fNeutralIsolationR(),
  fDoCellIsolation(kTRUE),
  fDoTagging(kFALSE),
  fPi0TaggingWindow(),
  fEtaTaggingWindow(),
  fSaveConversions(kTRUE),
  fSaveEMCClusters(kTRUE),
  fSavePHOSClusters(kTRUE),
  fSaveTracks(kTRUE),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoChargedIso(NULL),
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
  fClusterEMCalCandidatesBackground(0),
  fClusterPHOSCandidates(0),
  fTracks(0),
  fMCParticles(0),
  fExtraClusterInfo(),
  fExtraClusterInfoBackground(),
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
  fClusterCutsBackgroundEMC(NULL),
  fClusterCutsPHOS(NULL),
  fConvCuts(NULL),
  fCaloUtils(NULL),
  fMinClsTPC(0),
  fChi2PerClsTPC(9999),
  fMinClsITS(0),
  fEtaCut(9999),
  fPtCut(0),
  fYMCCut(99999),
  fMatchingParamsPhi(),
  fMatchingParamsEta(),
  fMatchingEOverP(),
  fDoBackgroundTrackMatching(kFALSE),
  fDoTrackIsolation(kFALSE),
  fTrackIsolationR(),
  fDoNeutralIsolation(kFALSE),
  fNeutralIsolationR(),
  fDoCellIsolation(kTRUE),
  fDoTagging(kFALSE),
  fPi0TaggingWindow(),
  fEtaTaggingWindow(),
  fSaveConversions(kTRUE),
  fSaveEMCClusters(kTRUE),
  fSavePHOSClusters(kTRUE),
  fSaveTracks(kTRUE),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoChargedIso(NULL),
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

  if(((AliCaloPhotonCuts*)fClusterCutsBackgroundEMC)->GetCutHistograms()){
    fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsBackgroundEMC)->GetCutHistograms());
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
  fHistoNEvents->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problem");
  fHistoNEvents->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
  fHistoNEvents->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
  fHistoNEvents->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
  fHistoNEvents->GetYaxis()->SetTitle("N_{events}");
  fOutputList->Add(fHistoNEvents);

  fHistoChargedIso           = new TH1F("fHistoChargedIso","fHistoChargedIso",500,-0.5,50);
  fOutputList->Add(fHistoChargedIso);
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
  fAnalysisTree->Branch("fExtraClusterInfo","vector<extraClusterInfo>",&fExtraClusterInfo,101);
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
  ((AliCaloPhotonCuts*)fClusterCutsEMC)->InitializeEMCAL(fInputEvent);
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
  for (auto p : fClusterEMCalCandidatesBackground){ delete p;} 
  fClusterEMCalCandidatesBackground.clear();
  for (auto p : fClusterPHOSCandidates){ delete p;} 
  fClusterPHOSCandidates.clear();
  for (auto p : fTracks){ delete p;} 
  fTracks.clear();
  for (auto p : fMCParticles){ delete p;} 
  fMCParticles.clear();
  fExtraClusterInfo.clear();
  fExtraClusterInfoBackground.clear();
  fConvIsoInfo.clear();
  fCaloIsoInfo.clear();

  fConversionCandidates.resize(0);
  fClusterEMCalCandidates.resize(0);
  fClusterPHOSCandidates.resize(0);
  fTracks.resize(0);
  fMCParticles.resize(0);
  fExtraClusterInfo.resize(0);
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
   if(fIsMC > 0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
     RelabelAODPhotonCandidates(kTRUE);    // In case of AODMC relabeling MC
     fV0Reader->RelabelAODs(kTRUE);
   }
   for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    if(!((AliConversionPhotonCuts*)fConvCuts)->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    fConversionCandidates.push_back(new AliAODConversionPhoton(*PhotonCandidate));
    
    isoInfo convIso;
    if(fDoTrackIsolation) ProcessChargedIsolation(PhotonCandidate,convIso.isoRawCharged);
    if(fDoNeutralIsolation) ProcessNeutralIsolation(PhotonCandidate,convIso.isoRawNeutral);
    if(fDoCellIsolation) ProcessCellIsolation(PhotonCandidate,convIso.isoCell);
    if(fDoTagging) convIso.isTagged = ProcessTagging(PhotonCandidate); 
    fConvIsoInfo.push_back(convIso);
    
    
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

        // get additional cluster info
        Short_t nLM = ((AliCaloPhotonCuts*)fClusterCutsEMC)->GetNumberOfLocalMaxima(clus, fInputEvent);
        Short_t matchIndex = ProcessTrackMatching(clus,fTracks);
        Float_t eFrac = GetExoticEnergyFraction(clus,fInputEvent);
        if(((AliCaloPhotonCuts*)fClusterCutsBackgroundEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
          
          extraClusterInfo clusInfo;  
          clusInfo.nLM = nLM;
          clusInfo.matchedTrackIndex = matchIndex;
          clusInfo.exoticEFrac = eFrac;
          
          fExtraClusterInfoBackground.push_back(clusInfo);
          fClusterEMCalCandidatesBackground.push_back(new AliAODCaloCluster(*clus));
        }

        // check if given EMC cuts are fulfilled
        if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
          delete clus;
          continue;
        }
        // get and set nLM  
        extraClusterInfo clusInfo;  
        clusInfo.nLM = nLM;
        clusInfo.matchedTrackIndex = matchIndex;
        clusInfo.exoticEFrac = eFrac;

        fExtraClusterInfo.push_back(clusInfo);
        fClusterEMCalCandidates.push_back(new AliAODCaloCluster(*clus));
        isoInfo caloIso;
        if(fDoTrackIsolation) ProcessChargedIsolation(clus,caloIso.isoRawCharged);
        if(fDoNeutralIsolation) ProcessNeutralIsolation(clus,caloIso.isoRawNeutral);
        if(fDoCellIsolation) ProcessCellIsolation(clus,caloIso.isoCell);
        if(fDoTagging) caloIso.isTagged = ProcessTagging(clus);  // TODO
        fCaloIsoInfo.push_back(caloIso);
        delete clus;
      }
  }

  // Loop over normal clusters
  for(Long_t i = 0; i < nclus; i++){
    Double_t tempClusterWeight              =  fWeightJetJetMC;                   
    clus                                = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    
    if(!clus) continue;

    if(!arrClustersProcess && clus->IsEMCAL()){ // if is was not saved already
      // get additional cluster info
      Short_t nLM = ((AliCaloPhotonCuts*)fClusterCutsEMC)->GetNumberOfLocalMaxima(clus, fInputEvent);
      Short_t matchIndex = ProcessTrackMatching(clus,fTracks);
      Float_t eFrac = GetExoticEnergyFraction(clus,fInputEvent);

      extraClusterInfo clusInfo;  
      clusInfo.nLM = nLM;
      clusInfo.matchedTrackIndex = matchIndex;
      clusInfo.exoticEFrac = eFrac;
      if(((AliCaloPhotonCuts*)fClusterCutsBackgroundEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        
        fExtraClusterInfoBackground.push_back(clusInfo);
        fClusterEMCalCandidatesBackground.push_back(new AliAODCaloCluster(*clus));
      }
      if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        delete clus;
        continue;
      }
      fExtraClusterInfo.push_back(clusInfo);
      fClusterEMCalCandidates.push_back(new AliAODCaloCluster(*clus));
      
      isoInfo caloIso;
      if(fDoTrackIsolation) ProcessChargedIsolation(clus,caloIso.isoRawCharged);
      if(fDoNeutralIsolation)ProcessNeutralIsolation(clus,caloIso.isoRawNeutral);
      if(fDoCellIsolation) ProcessCellIsolation(clus,caloIso.isoCell);
      if(fDoTagging) caloIso.isTagged = ProcessTagging(clus);  // TODO
      fCaloIsoInfo.push_back(caloIso);
      
      delete clus;
      continue;
    }
    if(clus->IsPHOS() && fSavePHOSClusters){
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
        if(TMath::Abs(dEta) > (fMatchingParamsEta[0] + pow(aodt->Pt() + fMatchingParamsEta[1],fMatchingParamsEta[2]))) continue;
        if(TMath::Abs(dPhi) > (fMatchingParamsPhi[0] + pow(aodt->Pt() + fMatchingParamsPhi[1],fMatchingParamsPhi[2]))) continue;
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
void AliAnalysisTaskGammaIsoTree::ProcessChargedIsolation(AliAODConversionPhoton* photon, Double32_t arrIso[]){
    TLorentzVector* v4photon = new TLorentzVector();
    arrIso[0] = 0.;
    arrIso[1] = 0.;
    

    v4photon->SetPxPyPzE(photon->Px(),photon->Py(),photon->Pz(),photon->E());
    for (Int_t t = 0; t < fInputEvent->GetNumberOfTracks(); t++)
    {
        AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
        if(!aodt) continue;
        if(!TrackIsSelectedAOD(aodt)) continue;
        TLorentzVector v4track;
        v4track.SetPxPyPzE(aodt->Px(),aodt->Py(),aodt->Pz(),aodt->E());
        
        Double_t trackEta = v4track.Eta();
        Double_t trackPhi = v4track.Phi();
        if (trackPhi < 0) trackPhi += 2*TMath::Pi();

        Double_t photonEta = v4photon->Eta();
        Double_t photonPhi = v4photon->Phi();
        if (photonPhi < 0) photonPhi += 2*TMath::Pi();

        Double_t dEta = trackEta - photonEta;
        Double_t dPhi = trackPhi - photonPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
        
        if((dR > fTrackIsolationR[0]) && (dR > fTrackIsolationR[1])) continue;

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
        if(dR <= fTrackIsolationR[0]) arrIso[0] += v4track.Pt();
        if(dR <= fTrackIsolationR[1]) arrIso[1] += v4track.Pt();
    }
    delete v4photon;
    return;
}

// charged isolation for clusters
//_____________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessChargedIsolation(AliAODCaloCluster* cluster, Double32_t arrIso[]){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    
    TLorentzVector v4cluster;
    cluster->GetMomentum(v4cluster,vertex);
    arrIso[0] = 0.;
    arrIso[1] = 0.;
    // cout << fInputEvent->GetNumberOfTracks() << endl;
    for (Int_t t = 0; t < fInputEvent->GetNumberOfTracks(); t++)
    {
        AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
        if(!aodt) continue;
        if(!TrackIsSelectedAOD(aodt)) continue;
        TLorentzVector v4track;
        v4track.SetPxPyPzE(aodt->Px(),aodt->Py(),aodt->Pz(),aodt->E());        
        Double_t trackEta = v4track.Eta();
        Double_t trackPhi = v4track.Phi();
        if (trackPhi < 0) trackPhi += 2*TMath::Pi();

        Double_t clusterEta = v4cluster.Eta();
        Double_t clusterPhi = v4cluster.Phi();
        if (clusterPhi < 0) clusterPhi += 2*TMath::Pi();

        Double_t dEta = trackEta - clusterEta;
        Double_t dPhi = trackPhi - clusterPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
        if(dR <= fTrackIsolationR[0]) arrIso[0] += v4track.Pt();
        if(dR <= fTrackIsolationR[1]) arrIso[1] += v4track.Pt();
    }
    fHistoChargedIso->Fill(arrIso[0]); // debug only
    return;
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessNeutralIsolation(AliAODConversionPhoton* photon, Double32_t arrIso[]){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
   
    TLorentzVector* v4photon = new TLorentzVector();
    arrIso[0] = 0.;
    arrIso[1] = 0.;
  
    v4photon->SetPxPyPzE(photon->Px(),photon->Py(),photon->Pz(),photon->E());
    for (UInt_t c = 0; c < fClusterEMCalCandidatesBackground.size(); c++)
    {
        AliAODCaloCluster* clusterE = (AliAODCaloCluster*) fClusterEMCalCandidatesBackground.at(c);
        if(!clusterE) continue;

        extraClusterInfo clusInfo = fExtraClusterInfoBackground.at(c);
        // check if cluster is neutral
        if(fDoBackgroundTrackMatching){
           if(clusInfo.matchedTrackIndex != -1) continue;
        }
     
    
        TLorentzVector v4cluster;
        clusterE->GetMomentum(v4cluster,vertex);
        Double_t photonEta = v4photon->Eta();
        Double_t photonPhi = v4photon->Phi();
        if (photonPhi < 0) photonPhi += 2*TMath::Pi();

        Double_t clusterEta = v4cluster.Eta();
        Double_t clusterPhi = v4cluster.Phi();
        if (clusterPhi < 0) clusterPhi += 2*TMath::Pi();

        Double_t dEta = photonEta - clusterEta;
        Double_t dPhi = photonPhi - clusterPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
        
        if(dR <= fNeutralIsolationR[0]) arrIso[0] += v4cluster.Et();
        if(dR <= fNeutralIsolationR[1]) arrIso[1] += v4cluster.Et();
    }
    delete v4photon;
    return;
}

//____Experimental isolation using EMC cells
void AliAnalysisTaskGammaIsoTree::ProcessCellIsolation(AliAODConversionPhoton* photon, Double32_t arrIso[]){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
   
    TLorentzVector* v4photon = new TLorentzVector();
    v4photon->SetPxPyPzE(photon->Px(),photon->Py(),photon->Pz(),photon->E());

    arrIso[0] = 0.;
    arrIso[1] = 0.;
  
    AliVCaloCells *cells = InputEvent()->GetEMCALCells();
    const Short_t nCells = cells->GetNumberOfCells();
  
    // count cells above threshold per sm
    Int_t bunchCrossNo = InputEvent()->GetBunchCrossNumber();
    for(Int_t iCell=0; iCell<nCells; ++iCell) {

              // Define necessary variables
        Short_t cellId                   = 0;
        Double_t cellE = 0,  cellTime = 0, cellEFrac = 0;
        Int_t cellMCLabel = 0;

        // Get Cell 
        cells->GetCell(iCell,cellId,cellE,cellTime,cellMCLabel,cellEFrac);


        UShort_t cellMax[]  = {(UShort_t) cellId};
        Bool_t   badCell    = GetCaloUtils()->GetEMCALRecoUtils()->ClusterContainsBadChannel(GetCaloUtils()->GetEMCALGeometry(),cellMax,1);
        if(badCell) continue;
        Bool_t   exoticCell = GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(cellId,cells,bunchCrossNo);
        if(exoticCell) continue;
        // Int_t sm       = cellId / (24*48);
        

        // energy cut
        if(cellE < 0.1) continue;

        Float_t cellEta = 0;
        Float_t cellPhi = 0;

        fGeomEMCAL->EtaPhiFromIndex(cellId,cellEta,cellPhi);
         
        Double_t photonEta = v4photon->Eta();
        Double_t photonPhi = v4photon->Phi();
        if (photonPhi < 0) photonPhi += 2*TMath::Pi();

        if (cellPhi < 0) cellPhi += 2*TMath::Pi();

        Double_t dEta = photonEta - cellEta;
        Double_t dPhi = photonPhi - cellPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
        
        if(dR <= fNeutralIsolationR[0]) arrIso[0] += cellE;
        if(dR <= fNeutralIsolationR[1]) arrIso[1] += cellE;
    }
    delete v4photon;
    return;
}

//____Experimental isolation using EMC cells
void AliAnalysisTaskGammaIsoTree::ProcessCellIsolation(AliAODCaloCluster* cluster, Double32_t arrIso[]){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
   
    TLorentzVector tmp;
    cluster->GetMomentum(tmp,vertex);
    
    TLorentzVector* v4thiscluster = new TLorentzVector(tmp);

    const Short_t nClusterCells = cluster->GetNCells();
    UShort_t* idClusterCells = cluster->GetCellsAbsId();

    arrIso[0] = 0.;
    arrIso[1] = 0.;
  
    AliVCaloCells *cells = InputEvent()->GetEMCALCells();
    const Short_t nCells = cells->GetNumberOfCells();
  
    // count cells above threshold per sm
    Int_t bunchCrossNo = InputEvent()->GetBunchCrossNumber();
   

    for(Int_t iCell=0; iCell<nCells; ++iCell) {
      
        Short_t cellId                   = 0;
        Double_t cellE = 0,  cellTime = 0, cellEFrac = 0;
        Int_t cellMCLabel = 0;

        // Get Cell 
        cells->GetCell(iCell,cellId,cellE,cellTime,cellMCLabel,cellEFrac);

        UShort_t cellMax[]  = {(UShort_t) cellId};
        Bool_t   badCell    = GetCaloUtils()->GetEMCALRecoUtils()->ClusterContainsBadChannel(GetCaloUtils()->GetEMCALGeometry(),cellMax,1);
        if(badCell) continue;
        Bool_t   exoticCell = GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(cellId,cells,bunchCrossNo);
        if(exoticCell) continue;
        // Int_t sm       = cellId / (24*48);

        // energy cut
        if(cellE < 0.1) continue;

        // check that cell is not contained in cluster
        Bool_t cellInCluster = kFALSE;
        for (Short_t c = 0; c < nClusterCells; c++)
        {
            if(idClusterCells[c] == cellId) cellInCluster = kTRUE;
        }

        if(cellInCluster) continue;

        Float_t cellEta = 0.;
        Float_t cellPhi = 0.;
        fGeomEMCAL->EtaPhiFromIndex(cellId,cellEta,cellPhi);

        Double_t clusterEta = v4thiscluster->Eta();
        Double_t clusterPhi = v4thiscluster->Phi();
        if (clusterPhi < 0) clusterPhi += 2*TMath::Pi();

        if (cellPhi < 0) cellPhi += 2*TMath::Pi();

        Double_t dEta = clusterEta - cellEta;
        Double_t dPhi = clusterPhi - cellPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
        
        if(dR <= fNeutralIsolationR[0]) arrIso[0] += cellE;
        if(dR <= fNeutralIsolationR[1]) arrIso[1] += cellE;
    }
    delete v4thiscluster;
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessNeutralIsolation(AliAODCaloCluster* cluster, Double32_t arrIso[]){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    arrIso[0] = 0.;
    arrIso[1] = 0.;

    TLorentzVector tmp;
    cluster->GetMomentum(tmp,vertex);
    
    TLorentzVector* v4thiscluster = new TLorentzVector(tmp);
    for (UInt_t c = 0; c < fClusterEMCalCandidatesBackground.size(); c++)
    {
        AliAODCaloCluster* clusterE = (AliAODCaloCluster*) fClusterEMCalCandidatesBackground.at(c);
        if(!clusterE) continue;
        if(clusterE->GetID() == cluster->GetID()) continue;
        extraClusterInfo clusInfo = fExtraClusterInfoBackground.at(c);
        if(fDoBackgroundTrackMatching){
           if(clusInfo.matchedTrackIndex != -1) continue;
        }
 
        TLorentzVector v4othercluster;
        clusterE->GetMomentum(v4othercluster,vertex);
        
        Double_t thisclusterEta = v4thiscluster->Eta();
        Double_t thisclusterPhi = v4thiscluster->Phi();
        if (thisclusterPhi < 0) thisclusterPhi += 2*TMath::Pi();

        Double_t otherclusterEta = v4othercluster.Eta();
        Double_t otherclusterPhi = v4othercluster.Phi();
        if (otherclusterPhi < 0) otherclusterPhi += 2*TMath::Pi();

        Double_t dEta = thisclusterEta - otherclusterEta;
        Double_t dPhi = thisclusterPhi - otherclusterPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

        if(dR <= fNeutralIsolationR[0]) arrIso[0] += v4othercluster.Et();
        if(dR <= fNeutralIsolationR[1]) arrIso[1] += v4othercluster.Et();
    }
    delete v4thiscluster;
    return;
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

  for (UInt_t c = 0; c < fClusterEMCalCandidatesBackground.size(); c++)
  {
    extraClusterInfo clusInfo = fExtraClusterInfoBackground.at(c);
    if(fDoBackgroundTrackMatching){
       if(clusInfo.matchedTrackIndex != -1) continue;
    }
    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    fClusterEMCalCandidatesBackground.at(c)->GetMomentum(clusterVector,vertex);

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


  for (UInt_t c = 0; c < fClusterEMCalCandidatesBackground.size(); c++)
  {
    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    
    if(fClusterEMCalCandidatesBackground.at(c)->GetID() == cluster->GetID()) continue;
    
    extraClusterInfo clusInfo = fExtraClusterInfoBackground.at(c);
    if(fDoBackgroundTrackMatching){
       if(clusInfo.matchedTrackIndex != -1) continue;
    }
    fClusterEMCalCandidatesBackground.at(c)->GetMomentum(clusterVector,vertex);

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

//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::RelabelAODPhotonCandidates(Bool_t mode){

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

Float_t AliAnalysisTaskGammaIsoTree::GetExoticEnergyFraction(AliVCluster *cluster, AliVEvent *event){
    Float_t exoticEnergyFrac =-999;   
    if (!cluster) {
      AliInfo("Cluster pointer null!");
      return 999;
    }
    AliVCaloCells* cells    = event->GetEMCALCells();

    Int_t largestCellID     = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindLargestCellInCluster(cluster,event);
    Float_t ecell1          = cells->GetCellAmplitude(largestCellID);
    Float_t eCross          = ((AliCaloPhotonCuts*)fClusterCutsEMC)->GetECross(largestCellID,cells);
    exoticEnergyFrac = 1-eCross/ecell1;
    
    return exoticEnergyFrac;

}

