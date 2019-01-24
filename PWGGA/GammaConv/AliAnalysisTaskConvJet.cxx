/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
 **************************************************************************/

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskConvJet.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskConvJet);
/// \endcond


AliAnalysisTaskConvJet::AliAnalysisTaskConvJet() :
  AliAnalysisTaskEmcalJet(),
  fIsMC(0),
  fNJetContainers(0),
  fNTrueJetContainers(0),
  fJetContainersAdded(0),
  fTrueJetContainersAdded(0),
  fJetNameArray(NULL),
  fTrueJetNameArray(NULL),
  fTrainconfigArray(NULL),
  fTrueTrainconfigArray(NULL),
  fListNJets(0),
  fListJetPt(0),
  fListJetPx(0),
  fListJetPy(0),
  fListJetPz(0),
  fListJetEta(0),
  fListJetPhi(0),
  fListJetArea(0),
  fListTrueNJets(0),
  fListTrueJetPt(0),
  fListTrueJetPx(0),
  fListTrueJetPy(0),
  fListTrueJetPz(0),
  fListTrueJetEta(0),
  fListTrueJetPhi(0),
  fListTrueJetArea(0)
{
}

AliAnalysisTaskConvJet::AliAnalysisTaskConvJet(const char *name, Int_t IsMC) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fIsMC(IsMC),
  fNJetContainers(0),
  fNTrueJetContainers(0),
  fJetContainersAdded(0),
  fTrueJetContainersAdded(0),
  fJetNameArray(NULL),
  fTrueJetNameArray(NULL),
  fTrainconfigArray(NULL),
  fTrueTrainconfigArray(NULL),
  fListNJets(0),
  fListJetPt(0),
  fListJetPx(0),
  fListJetPy(0),
  fListJetPz(0),
  fListJetEta(0),
  fListJetPhi(0),
  fListJetArea(0),
  fListTrueNJets(0),
  fListTrueJetPt(0),
  fListTrueJetPx(0),
  fListTrueJetPy(0),
  fListTrueJetPz(0),
  fListTrueJetEta(0),
  fListTrueJetPhi(0),
  fListTrueJetArea(0)
{
  SetMakeGeneralHistograms(kTRUE);
}
/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskConvJet::UserCreateOutputObjects()
{
   AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/**
 * The body of this function should contain instructions to fill the vectors.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskConvJet::FillHistograms()
{
  DoJetLoop();

  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the vectors.
 */
void AliAnalysisTaskConvJet::DoJetLoop()
{
  AliJetContainer* jetCont = 0;
  vector<Int_t> MatchRec;
  vector<Int_t> MatchTrue;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    MatchRec.clear();
    MatchTrue.clear();
    TString JetName = jetCont->GetTitle();
    TObjArray *arr = JetName.Tokenize("__");
    TObjString* testObjString = (TObjString*)arr->At(2);
    if(testObjString->GetString() != "mcparticles"){
      for(Int_t i = 0; i < fNJetContainers; i++){
        if(JetName == fJetNameArray[i]){
          MatchRec.push_back(i);
        }
      }
      for(UInt_t k = 0; k < MatchRec.size(); k++){
        Int_t j = MatchRec.at(k);
        UInt_t count = 0;
        fListJetPt.at(j).clear();
        fListJetPx.at(j).clear();
        fListJetPy.at(j).clear();
        fListJetPz.at(j).clear();
        fListJetEta.at(j).clear();
        fListJetPhi.at(j).clear();
        fListJetArea.at(j).clear();
        for(auto jet : jetCont->accepted()) {
          if (!jet) continue;
          count++;
          fListJetPt.at(j).push_back(jet->Pt());
          fListJetPx.at(j).push_back(jet->Px());
          fListJetPy.at(j).push_back(jet->Py());
          fListJetPz.at(j).push_back(jet->Pz());
          fListJetEta.at(j).push_back(jet->Eta());
          fListJetPhi.at(j).push_back(jet->Phi());
          fListJetArea.at(j).push_back(jet->Area());
        }
        fListNJets.at(j) = count;
      }
    }else{
      for(Int_t i = 0; i < fNTrueJetContainers; i++){
        if(JetName == fTrueJetNameArray[i]){
            MatchTrue.push_back(i);
        }
      }
      for(UInt_t k = 0; k < MatchTrue.size(); k++){
        Int_t j = MatchTrue.at(k);
        UInt_t count = 0;
        fListTrueJetPt.at(j).clear();
        fListTrueJetPx.at(j).clear();
        fListTrueJetPy.at(j).clear();
        fListTrueJetPz.at(j).clear();
        fListTrueJetEta.at(j).clear();
        fListTrueJetPhi.at(j).clear();
        fListTrueJetArea.at(j).clear();
        for(auto jet : jetCont->accepted()) {
          if (!jet) continue;
          count++;
          fListTrueJetPt.at(j).push_back(jet->Pt());
          fListTrueJetPx.at(j).push_back(jet->Px());
          fListTrueJetPy.at(j).push_back(jet->Py());
          fListTrueJetPz.at(j).push_back(jet->Pz());
          fListTrueJetEta.at(j).push_back(jet->Eta());
          fListTrueJetPhi.at(j).push_back(jet->Phi());
          fListTrueJetArea.at(j).push_back(jet->Area());
        }
        fListTrueNJets.at(j) = count;
      }
    }
  }
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskConvJet::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskConvJet::Run()
{
  return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskConvJet::Terminate(Option_t *) 
{
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskConvJet * AliAnalysisTaskConvJet::AddTask_GammaConvJet(
  const char *ntracks,
  const char *nclusters,
  const char *ncells,
  const char *suffix,
  Int_t IsMC,
  Int_t NContainers)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTask_GammaConvJet", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTask_GammaConvJet", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(ntracks);
  TString clusName(nclusters);
  TString cellName(ncells);

  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  if (cellName == "usedefault") {
    if (dataType == kESD) {
      cellName = "EMCALCells";
    }
    else if (dataType == kAOD) {
      cellName = "emcalCells";
    }
    else {
      cellName = "";
    }
  }

  TString name("AliAnalysisTaskConvJet");

  AliAnalysisTaskConvJet* sampleTask = new AliAnalysisTaskConvJet(name, IsMC);
  sampleTask->SetCaloCellsName(cellName);
  sampleTask->SetVzRange(-10,10);

  sampleTask->SetNumberOfContainers(NContainers);

  for(Int_t i = 0; i < NContainers; i++){

    if (trackName == "mcparticles") {
      sampleTask->AddMCParticleContainer(trackName);
    }
    else if (trackName == "tracks" || trackName == "Tracks") {
      sampleTask->AddTrackContainer(trackName);
    }
    else if (!trackName.IsNull()) {
     sampleTask->AddParticleContainer(trackName);
    }
    sampleTask->AddClusterContainer(clusName);
  
    sampleTask->GetClusterContainer(i)->SetClusECut(0.);
    sampleTask->GetClusterContainer(i)->SetClusPtCut(0.);
    sampleTask->GetClusterContainer(i)->SetClusNonLinCorrEnergyCut(0.);
    sampleTask->GetClusterContainer(i)->SetClusHadCorrEnergyCut(0.30);
    sampleTask->GetClusterContainer(i)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    sampleTask->GetParticleContainer(i)->SetParticlePtCut(0.15);
    sampleTask->GetParticleContainer(i)->SetParticleEtaLimits(-0.8,0.8);

    if(trackName != "mcparticles"){
      sampleTask->GetTrackContainer(i)->SetFilterHybridTracks(kTRUE);
      sampleTask->GetTrackContainer(i)->SetParticlePtCut(0.15);
      sampleTask->GetTrackContainer(i)->SetParticleEtaLimits(-0.8,0.8);
    }
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(sampleTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(trackName);
  contname += suffix;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (sampleTask, 0,  cinput1 );
  mgr->ConnectOutput (sampleTask, 1, coutput1 );

  return sampleTask;
}
