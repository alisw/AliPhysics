/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <vector>

#include <TClonesArray.h>
#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliVEvent.h>
#include "AliEmcalJet.h"
#include "AliEmcalJetFinderKernel.h"
#include "AliEmcalJetTaskV1.h"

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::AliEmcalJetTaskV1);
/// \endcond

namespace PWGJE {

namespace EMCALJetTasks {

/**
 * Default constructor. This constructor is only for ROOT I/O and
 * not to be used by users.
 */
AliEmcalJetTaskV1::AliEmcalJetTaskV1() :
  AliAnalysisTaskEmcal(),
  fJetFinder(nullptr),
  fJets(nullptr),
  fJetsTag(),
  fJetsName()
{
}

/**
 * Standard named constructor.
 * @param name Name of the task.
 */
AliEmcalJetTaskV1::AliEmcalJetTaskV1(const char *name) :
  AliAnalysisTaskEmcal(name),
  fJetFinder(nullptr),
  fJets(nullptr),
  fJetsTag("Jets"),
  fJetsName()
{
}

/**
 * Destructor
 */
AliEmcalJetTaskV1::~AliEmcalJetTaskV1()
{
  if(fJetFinder) delete fJetFinder;
}


/**
 * This method is called for each event.
 * @return Always kTRUE
 */
Bool_t AliEmcalJetTaskV1::Run()
{
  // clear the jet array (normally a null operation)
  fJets->Delete();
  fJetFinder->RunJetFinder(fJets);
  Int_t n = fJets->GetEntriesFast();

  if (n == 0) return kFALSE;

  return kTRUE;
}

/**
 * This method is called once before analzying the first event.
 * It generates the output jet branch name, initializes the FastJet wrapper
 * and the utilities (FJ contribs).
 */
void AliEmcalJetTaskV1::ExecOnce()
{
  fJetsName = AliJetContainer::GenerateJetName(fJetFinder->GetJetType(), fJetFinder->GetJetAlgo(), fJetFinder->GetRecombScheme(), fJetFinder->GetRadius(), fJetFinder->GetParticleContainer(0), fJetFinder->GetClusterContainer(0), fJetsTag);
  std::cout << GetName() << ": Name of the jet container: " << fJetsName << std::endl;
  std::cout << "Use this name in order to connect jet containers in your task to connect to the collection of jets found by this jet finder" << std::endl;
  if(auto partcont = GetParticleContainer(0)) {
    std::cout << "Found particle container with name " << partcont->GetName() << std::endl;
  } else {
    std::cout << "Not particle container found for task" << std::endl;
  }
  if(auto clustcont = GetClusterContainer(0)){
    std::cout << "Found cluster container with name " << clustcont->GetName() << std::endl;
  } else {
    std::cout << "Not cluster container found for task" << std::endl;
  }

  // add jets to event if not yet there
  if (!(InputEvent()->FindListObject(fJetsName))) {
    fJets = new TClonesArray("AliEmcalJet");
    fJets->SetName(fJetsName);
    ::Info("AliEmcalJetTask::ExecOnce", "Jet collection with name '%s' has been added to the event.", fJetsName.Data());
    InputEvent()->AddObject(fJets);
  }
  else {
    AliError(Form("%s: Object with name %s already in event! Returning", GetName(), fJetsName.Data()));
    return;
  }

  AliAnalysisTaskEmcal::ExecOnce();

  fJetFinder->Init();

  // Set the particle containers
  for(auto pcont : fParticleCollArray)
    fJetFinder->AppendParticleContainer(static_cast<AliParticleContainer *>(pcont));
  for(auto ccont : fClusterCollArray)
    fJetFinder->AppendClusterContainer(static_cast<AliClusterContainer *>(ccont));

  fJetFinder->InitIndexMaps();
}

void AliEmcalJetTaskV1::RunChanged(Int_t runnumber) {
  fJetFinder->SetRunNumber(runnumber);
}

/**
 * Add an instance of this class to the analysis manager
 * @param nTracks name of the track collection
 * @param nClusters name of the calorimeter cluster collection
 * @param jetAlgo jet finding algorithm (anti-kt, kt, etc.)
 * @param radius jet resolution parameter (0.2, 0.4, tyc.)
 * @param jetType full, charged or neutral
 * @param minTrPt cut on the minimum transverse momentum of tracks
 * @param minClPt cut on the minimum transverse momentum of calorimeter clusters
 * @param ghostArea area of ghost particles (determines the jet area resolution)
 * @param reco recombination scheme
 * @param tag addtional information to be appended at the end of the output jet collection name
 * @param minJetPt cut on the minimum jet pt
 * @param lockTask lock the task - no further changes are possible if kTRUE
 * @param bFillGhosts add ghosts particles among the jet constituents in the output
 * @return a pointer to the new AliEmcalJetTask instance
 */
AliEmcalJetTaskV1* AliEmcalJetTaskV1::AddTaskEmcalJet(
  const TString nTracks, const TString nClusters,
  const AliJetContainer::EJetAlgo_t jetAlgo, const Double_t radius, const AliJetContainer::EJetType_t jetType,
  const Double_t minTrPt, const Double_t minClPt,
  const Double_t ghostArea, const AliJetContainer::ERecoScheme_t reco,
  const TString tag, const Double_t minJetPt,
  const Bool_t lockTask, const Bool_t bFillGhosts
)
{
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalJet", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AddTaskEmcalJet", "This task requires an input event handler");
    return 0;
  }

  EDataType_t dataType = kUnknownDataType;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(nTracks);
  TString clusName(nClusters);

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

  AliParticleContainer* partCont = 0;
  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = new AliMCParticleContainer(trackName);
    partCont = mcpartCont;
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = new AliTrackContainer(trackName);
    partCont = trackCont;
  }
  else if (!trackName.IsNull()) {
    partCont = new AliParticleContainer(trackName);
  }
  if (partCont) partCont->SetParticlePtCut(minTrPt);

  AliClusterContainer* clusCont = 0;
  if (!clusName.IsNull()) {
    clusCont = new AliClusterContainer(clusName);
    clusCont->SetClusECut(0.);
    clusCont->SetClusPtCut(0.);
    clusCont->SetClusHadCorrEnergyCut(minClPt);
    clusCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  switch (jetType) {
  case AliJetContainer::kChargedJet:
    if (partCont) partCont->SetCharge(AliParticleContainer::kCharged);
    break;
  case AliJetContainer::kNeutralJet:
    if (partCont) partCont->SetCharge(AliParticleContainer::kNeutral);
    break;
  default:
    break;
  }

  TString name = AliJetContainer::GenerateJetName(jetType, jetAlgo, reco, radius, partCont, clusCont, tag);

  Printf("Jet task name: %s", name.Data());

  auto mgrTask = static_cast<AliEmcalJetTaskV1 *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  auto jetTask = new AliEmcalJetTaskV1(name);
  jetTask->GetJetFinder()->SetJetType(jetType);
  jetTask->GetJetFinder()->SetJetAlgo(jetAlgo);
  jetTask->GetJetFinder()->SetRecombScheme(reco);
  jetTask->GetJetFinder()->SetRadius(radius);
  if (partCont) jetTask->AdoptParticleContainer(partCont);
  if (clusCont) jetTask->AdoptClusterContainer(clusCont);
  jetTask->SetJetsTag(tag);
  jetTask->GetJetFinder()->SetMinJetPt(minJetPt);
  jetTask->GetJetFinder()->SetGhostArea(ghostArea);

  if (bFillGhosts) jetTask->GetJetFinder()->SetFillGhost();
  if (lockTask) jetTask->GetJetFinder()->SetLocked();

  // Final settings, pass to manager and set the containers

  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(jetTask, 0, cinput);

  return jetTask;
}

}

}
