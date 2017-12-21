/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisTaskRhoDev.h"

#include <TClonesArray.h>
#include <TMath.h>

#include <AliLog.h>
#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>

#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskRhoDev);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskRhoDev::AliAnalysisTaskRhoDev() :
  AliAnalysisTaskRhoBaseDev(),
  fNExclLeadJets(0),
  fRhoSparse(kFALSE),
  fExclJetOverlap(),
  fOccupancyFactor(0),
  fHistOccCorrvsCent(nullptr)
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name  Name of the task
 * @param[in] histo If kTRUE, the task will also produce QA histograms
 */
AliAnalysisTaskRhoDev::AliAnalysisTaskRhoDev(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBaseDev(name, histo),
  fNExclLeadJets(0),
  fRhoSparse(kFALSE),
  fExclJetOverlap(),
  fOccupancyFactor(0),
  fHistOccCorrvsCent(nullptr)
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskRhoDev::UserCreateOutputObjects()
{
  if (!fCreateHisto) return;

  AliAnalysisTaskRhoBaseDev::UserCreateOutputObjects();

  fHistOccCorrvsCent = new TH2F("fHistOccCorrvsCent", "fHistOccCorrvsCent;Centrality (%);#it{C}", 100, 0, 100, 2000, 0 , 2);
  fOutput->Add(fHistOccCorrvsCent);
}

/**
 * Finds the first two leading jets
 * @return A pair with the leading and sub-leading jets respectively as first and second element
 */
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

/**
 * Calculates the average background using the median approach
 * as proposed in https://arxiv.org/pdf/0707.1378.pdf.
 * Rho is stored in fOutRho.
 */
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

/**
 * Fill histograms.
 */
Bool_t AliAnalysisTaskRhoDev::FillHistograms()
{
  Bool_t r = AliAnalysisTaskRhoBaseDev::FillHistograms();
  if (!r) return kFALSE;

  fHistOccCorrvsCent->Fill(fCent, fOccupancyFactor);

  return kTRUE;
}

/**
 * Verify that the required particle, cluster and jet containers were provided.
 * @return kTRUE if all requirements are satisfied, kFALSE otherwise
 */
Bool_t AliAnalysisTaskRhoDev::VerifyContainers()
{
  if (fJetCollArray.count("Background") == 0) {
    AliError("No signal jet collection found. Task will not run!");
    return kFALSE;
  }

  return kTRUE;
}

/**
 * Create an instance of this class and add it to the analysis manager
 * @param trackName name of the track collection
 * @param trackPtCut minimum pt of the tracks
 * @param clusName name of the calorimeter cluster collection
 * @param clusECut minimum energy of the calorimeter clustuers
 * @param nRho name of the output rho object
 * @param jetradius Radius of the kt jets used to calculate the background
 * @param acceptance Fiducial acceptance of the kt jets
 * @param jetType Jet type (full/charged)
 * @param rscheme Recombination scheme
 * @param histo If kTRUE the task will also produce QA histograms
 * @param suffix additional suffix that can be added at the end of the task name
 * @return pointer to the new AliAnalysisTaskRhoDev task
 */
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
