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
#include "AliAnalysisTaskRhoDev.h"

#include <TClonesArray.h>
#include <TMath.h>

#include <AliLog.h>
#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>

#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"

ClassImp(AliAnalysisTaskRhoDev);

AliAnalysisTaskRhoDev::AliAnalysisTaskRhoDev() :
  AliAnalysisTaskRhoBaseDev(),
  fNExclLeadJets(0),
  fRhoSparse(kFALSE),
  fExclJetOverlap(),
  fOccupancyFactor(0),
  fHistOccCorrvsCent(nullptr)
{
}

AliAnalysisTaskRhoDev::AliAnalysisTaskRhoDev(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBaseDev(name, histo),
  fNExclLeadJets(0),
  fRhoSparse(kFALSE),
  fExclJetOverlap(),
  fOccupancyFactor(0),
  fHistOccCorrvsCent(nullptr)
{
}

void AliAnalysisTaskRhoDev::UserCreateOutputObjects()
{
  if (!fCreateHisto) return;

  AliAnalysisTaskRhoBaseDev::UserCreateOutputObjects();

  fHistOccCorrvsCent = new TH2F("fHistOccCorrvsCent", "fHistOccCorrvsCent;Centrality (%);#it{C}", 100, 0, 100, 2000, 0 , 2);
  fOutput->Add(fHistOccCorrvsCent);
}

std::pair<AliEmcalJet*, AliEmcalJet*> AliAnalysisTaskRhoDev::GetLeadingJets()
{
  std::pair<AliEmcalJet*, AliEmcalJet*> maxJets = {nullptr, nullptr};
  if (fNExclLeadJets <= 0) return maxJets;

  auto itJet = fSortedJets["Background"].begin();

  maxJets.first = *itJet;
  if (fNExclLeadJets > 1) {
    itJet++;
    if (itJet != fSortedJets["Background"].end()) maxJets.second = *itJet;
  }

  return maxJets;
}

void AliAnalysisTaskRhoDev::CalculateRho()
{
  if (fJetCollArray.empty()) return;

  auto maxJets = GetLeadingJets();

  static Double_t rhovec[999];
  Int_t NjetAcc = 0;
  Double_t TotaljetArea = 0; // Total area of background jets (including ghost jets)
  Double_t TotaljetAreaPhys = 0; // Total area of physical background jets (excluding ghost jets)
  // Ghost jet is a jet made only of ghost particles

  AliJetContainer* bkgJetCont = fJetCollArray["Background"];
  AliJetContainer* sigJetCont = nullptr;
  if (!fExclJetOverlap.IsNull()) {
    auto sigJetContIt = fJetCollArray.find(fExclJetOverlap.Data());
    if (sigJetContIt != fJetCollArray.end()) sigJetCont = sigJetContIt->second;
  }

  // push all jets within selected acceptance into stack
  for (auto jet : bkgJetCont->accepted()) {

    TotaljetArea += jet->Area();

    if (jet->IsGhost()) continue;

    TotaljetAreaPhys += jet->Area();

    // excluding leading jets
    if (jet == maxJets.first || jet == maxJets.second) continue;

    Bool_t overlapsWithSignal = kFALSE;
    if (sigJetCont) {
      for (auto sigJet : sigJetCont->accepted()) {
        if (AreJetsOverlapping(jet, sigJet)) {
          overlapsWithSignal = kTRUE;
          break;
        }
      }
    }

    if (overlapsWithSignal) continue;

    rhovec[NjetAcc] = jet->Pt() / jet->Area();
    ++NjetAcc;
  }

  // Occupancy correction for sparse event described in https://arxiv.org/abs/1207.2392
  if (TotaljetArea > 0) {
    fOccupancyFactor = TotaljetAreaPhys / TotaljetArea;
  }
  else {
    fOccupancyFactor = 0;
  }

  if (NjetAcc > 0) {
    //find median value
    Double_t rho = TMath::Median(NjetAcc, rhovec);

    if (fRhoSparse) rho = rho * fOccupancyFactor;

    fOutRho->SetVal(rho);
  }
}

Bool_t AliAnalysisTaskRhoDev::FillHistograms()
{
  Bool_t r = AliAnalysisTaskRhoBaseDev::FillHistograms();
  if (!r) return kFALSE;

  fHistOccCorrvsCent->Fill(fCent, fOccupancyFactor);

  return kTRUE;
}

Bool_t AliAnalysisTaskRhoDev::VerifyContainers()
{
  if (fJetCollArray.count("Background") == 0) {
    AliError("No signal jet collection found. Task will not run!");
    return kFALSE;
  }

  return kTRUE;
}

AliAnalysisTaskRhoDev* AliAnalysisTaskRhoDev::AddTaskRhoDev(TString trackName, Double_t trackPtCut, TString clusName, Double_t clusECut, TString nRho, Double_t jetradius, UInt_t acceptance, AliJetContainer::EJetType_t jetType, AliJetContainer::ERecoScheme_t rscheme, Bool_t histo, TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisTaskRhoDev::AddTaskRhoDev", "No analysis manager to connect to.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AliAnalysisTaskRhoDev::AddTaskRhoDev", "This task requires an input event handler");
    return nullptr;
  }

  EDataType_t dataType = kUnknownDataType;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  // Init the task and do settings
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

  TString name(TString::Format("AliAnalysisTaskRhoDev_%s", nRho.Data()));
  if (!suffix.IsNull()) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskRhoDev* mgrTask = dynamic_cast<AliAnalysisTaskRhoDev*>(mgr->GetTask(name.Data()));
  if (mgrTask) {
    ::Warning("AliAnalysisTaskRhoDev::AddTaskRhoDev", "Not adding the task again, since a task with the same name '%s' already exists", name.Data());
    return mgrTask;
  }

  AliAnalysisTaskRhoDev* rhotask = new AliAnalysisTaskRhoDev(name, histo);
  rhotask->SetOutRhoName(nRho);

  AliParticleContainer* partCont = rhotask->AddParticleContainer(trackName.Data());
  partCont->SetMinPt(trackPtCut);
  AliClusterContainer *clusterCont = rhotask->AddClusterContainer(clusName.Data());
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(clusECut);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  AliJetContainer *jetCont = new AliJetContainer(jetType, AliJetContainer::kt_algorithm, rscheme, jetradius, partCont, clusterCont);
  if (jetCont) {
    jetCont->SetJetPtCut(0);
    jetCont->SetJetAcceptanceType(acceptance);
    jetCont->SetName("Background");
    rhotask->AdoptJetContainer(jetCont);
  }

  // Final settings, pass to manager and set the containers
  mgr->AddTask(rhotask);

  // Create containers for input/output
  mgr->ConnectInput(rhotask, 0, mgr->GetCommonInputContainer());
  if (histo) {
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
        TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(rhotask, 1, coutput1);
  }

  return rhotask;
}
