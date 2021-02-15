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
#include <array>
#include <bitset>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <TArrayD.h>
#include <TClonesArray.h>
#include <TGrid.h>
#include <THashList.h>
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParameter.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliClusterContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalList.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerDecision.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalTriggerSelectionCuts.h"
#include "AliEmcalTriggerStringDecoder.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliVVertex.h"

#include "AliAnalysisTaskEmcalClustersRef.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalClustersRef)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef() :
    AliAnalysisTaskEmcalTriggerBase(),
    AliAnalysisEmcalTriggerSelectionHelperImpl(),
    fCentralityRange(-999., 999.),
    fRequestCentrality(false),
    fEventCentrality(-1),
    fCentralityEstimator("V0M"),
    fBunchCrossingIndex(-1),
    fEnergyDefinition(kDefaultEnergy),
    fEnableSumw2(false),
    fDoFillMultiplicityHistograms(false),
    fUseFiredTriggers(false),
    fUseExclusiveTriggers(true),
    fFillTriggerClusters(true),
    fMonitorEtaPhi(false),
    fClusterTimeRange(-50e-6, 50e-6),
    fTriggerClusters(),
    fRequiredOverlaps(),
    fExcludedOverlaps()
{
}

AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef(const char *name) :
    AliAnalysisTaskEmcalTriggerBase(name),
    fCentralityRange(-999., 999.),
    fRequestCentrality(false),
    fEventCentrality(-1),
    fCentralityEstimator("V0M"),
    fBunchCrossingIndex(-1),
    fEnergyDefinition(kDefaultEnergy),
    fEnableSumw2(false),
    fDoFillMultiplicityHistograms(false),
    fUseFiredTriggers(false),
    fUseExclusiveTriggers(true),
    fFillTriggerClusters(true),
    fMonitorEtaPhi(false),
    fClusterTimeRange(-50e-6, 50e-6),
    fTriggerClusters(),
    fRequiredOverlaps(),
    fExcludedOverlaps()
{
}

AliAnalysisTaskEmcalClustersRef::~AliAnalysisTaskEmcalClustersRef() {
}


void AliAnalysisTaskEmcalClustersRef::CreateUserHistos(){

  EnergyBinning energybinning;
  TLinearBinning smbinning(21, -0.5, 20.5), smbinningUF(22, -1.5, 20.5), detbinning(3, -0.5, 2.5), etabinning(100, -0.7, 0.7), phibinning(200, 0., TMath::TwoPi()), timebinning(1000, -500e-9, 500e-9), ncellbinning(101, -0.5, 100.5); 
  TLinearBinning trgclustbinning(kTrgClusterN, -0.5, kTrgClusterN - 0.5);
  TString optionstring = fEnableSumw2 ? "s" : "";

  /*
   * Exclusive classes are defined as triggered events
   * in a class without lower threshold classes firing.
   * This is needed to make the triggers statistically
   * independent.
   */

  // Binnings for Multiplicity correlation
  TLinearBinning v0abinning(1000, 0., 1000.), centralitybinning(100, 0., 100.), trackletbinning(500, 0., 500.), itsclustbinning(500, 0., 500.), 
                 emcclustbinning(100, 0., 100.), emccellbinning(3000, 0., 3000.), adcbinning(2000, 0., 2000.);
  const TBinning *multbinning[6] = {&v0abinning, &trackletbinning, &trackletbinning, &itsclustbinning, &emcclustbinning, &emccellbinning};
  std::vector<const TBinning *> clusterallbinning, clustermaxbinning;
  clusterallbinning.push_back(&smbinning);
  clusterallbinning.push_back(&centralitybinning);
  clusterallbinning.push_back(&energybinning);
  if(fMonitorEtaPhi) {
    clusterallbinning.push_back(&etabinning);
    clusterallbinning.push_back(&phibinning);
  }
  // Max cluster: Must be separated between EMCAL and DCAL,
  // meaning for each event we fill the THnSparse with one
  // entry for EMCAL and one entry for DCAL (2 clusters per event).
  // It can be that there events where there is no cluster 
  // found in either of the detectors (or both). Those are
  // still counted as 0 but cannot be asigned to any detector,
  // so the SM ID is -1. In order to be able to distinguish between
  // EMCAL and DCAL a further axis with the detector index 
  // (0-EMCAL, 1-DCAL) is added.
  clusterallbinning.push_back(&trgclustbinning);
  clustermaxbinning.push_back(&smbinningUF);      // Add bin for -1 for events without clusters in detector
  clustermaxbinning.push_back(&detbinning);       // Add bin for Detector (EMCAL/DCAL), needed in order to separate events without clusters in only one of the detectors
  clustermaxbinning.push_back(&centralitybinning);
  clustermaxbinning.push_back(&energybinning);
  if(fMonitorEtaPhi) {
    clustermaxbinning.push_back(&etabinning);
    clustermaxbinning.push_back(&phibinning);
  }
  clustermaxbinning.push_back(&trgclustbinning);
  AliDebugStream(1) << "Using exclusive triggers: " << (fUseExclusiveTriggers ? "yes" : "no") << std::endl;
  for(auto trg : GetSupportedTriggers(fUseExclusiveTriggers)){
    AliDebugStream(1) << "Creating histograms for trigger " << trg << std::endl;
    fHistos->CreateTH1("hTrgClustCounter" + trg, "Event counter in trigger cluster " + trg, trgclustbinning, optionstring);
    fHistos->CreateTH1("hEventCentrality" + trg, "Event centrality for trigger class " + trg, 103, -2., 101., optionstring);
    fHistos->CreateTH1("hVertexZ" + trg, "z-position of the primary vertex for trigger class " + trg, 200, -40., 40., optionstring);
    if(this->fDoFillMultiplicityHistograms) fHistos->CreateTHnSparse("hMultiplicityCorrelation" + trg, "Multiplicity correlation for trigger" + trg, 6, multbinning, optionstring);
    fHistos->CreateTHnSparse("hClusterTHnSparseAll" + trg, "Cluster THnSparse (all) for trigger" + trg, clusterallbinning.size(), clusterallbinning.data(), optionstring);
    fHistos->CreateTHnSparse("hClusterTHnSparseMax" + trg, "Cluster THnSparse (max) for trigger" + trg, clustermaxbinning.size(), clustermaxbinning.data(), optionstring);
    if(fUseFiredTriggers) fHistos->CreateTHnSparse("hClusterTHnSparseFired" + trg, "Cluster THnSparse (firing) for trigger" + trg, clusterallbinning.size(), clusterallbinning.data(), optionstring);
    fHistos->CreateTH2("hTimeEnergy" + trg, "Cluster time vs. energy for trigger class " + trg, timebinning, energybinning, optionstring);
    fHistos->CreateTH2("hNCellEnergy" + trg, "Cluster number of cells vs energy for trigger class " + trg, ncellbinning, energybinning, optionstring);
    if(fUseFiredTriggers){
      fHistos->CreateTH2("hCorrClusterEPatchADC" + trg, "Correlation between cluster E and patch ADC for trigger " + trg, energybinning, adcbinning, optionstring);
      fHistos->CreateTH2("hCorrClusterEPatchE" + trg, "Correlation between cluster E and patch E for trigger " + trg, energybinning, energybinning, optionstring);
    }
  }
}

bool AliAnalysisTaskEmcalClustersRef::IsUserEventSelected(){
  fEventCentrality = 99;   // without centrality put everything in the peripheral bin
  if(fRequestCentrality){
    AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
    if(!mult){
      AliErrorStream() << GetName() << ": Centrality selection enabled but no centrality estimator found" << std::endl;
      return false;
    }
    //if(mult->IsEventSelected()) return false;
    fEventCentrality = mult->GetEstimator(fCentralityEstimator)->GetPercentile();
    AliDebugStream(1) << GetName() << ": Centrality " <<  fEventCentrality << std::endl;
    if(!fCentralityRange.IsInRange(fEventCentrality)){
      AliDebugStream(1) << GetName() << ": reject centrality: " << fEventCentrality << std::endl;
      return false;
    } else {
      AliDebugStream(1) << GetName() << ": select centrality " << fEventCentrality << std::endl;
    }
  } else {
    AliDebugStream(1) << GetName() << ": No centrality selection applied" << std::endl;
  }

  if(fBunchCrossingIndex > -1){
    int bcindex = fInputEvent->GetHeader()->GetBunchCrossNumber() % 4;
    if(bcindex != fBunchCrossingIndex) return false;
  }

  // handle overlaps
  if(fRequiredOverlaps.GetEntries() || fExcludedOverlaps.GetEntries()){
    if(fRequiredOverlaps.GetEntries()){
      Bool_t allFound(true);
      for(auto t : fRequiredOverlaps){
        auto trgstr = static_cast<TObjString *>(t)->String();
        if(std::find_if(fSelectedTriggers.begin(), fSelectedTriggers.end(), [&trgstr](const TString &seltrigger) -> bool { return seltrigger.Contains(trgstr); }) == fSelectedTriggers.end()) {
          allFound = false;
          break;
        }
      }
      if(!allFound) return false;
    }

    if(fExcludedOverlaps.GetEntries()){
      Bool_t oneFound(false);
      for(auto t : fExcludedOverlaps){
        auto trgstr = static_cast<TObjString *>(t)->String();
        if(std::find_if(fSelectedTriggers.begin(), fSelectedTriggers.end(), [&trgstr](const TString &seltrigger) -> bool { return seltrigger.Contains(trgstr); }) != fSelectedTriggers.end()) {
          oneFound = true;
          break;
        }
      }
      if(oneFound) return false;
    }
  }

  // determine trigger clusters (ANY - 0 always included)
  // Only meaningfull in data
  fTriggerClusters.clear();
  fTriggerClusters.emplace_back(kTrgClusterANY);
  if(fFillTriggerClusters) {
    auto triggerclusters =  GetTriggerClusterIndices(fInputEvent->GetFiredTriggerClasses().Data());
    for(auto en : triggerclusters) {
      if(std::find(fTriggerClusters.begin(), fTriggerClusters.end(), en) == fTriggerClusters.end()) fTriggerClusters.emplace_back(en);
    }
  }

  return true;
}

bool AliAnalysisTaskEmcalClustersRef::Run(){
  AliDebugStream(1) << GetName() << ": UserExec start" << std::endl;

  std::map<TString, const TList *> patchhandlers;
  const std::vector<TString> l1triggers = {"EJ1", "EJ2", "EG1", "EG2", "DJ1", "DJ2", "DG1", "DG2"};
  int energycomp = -1;
  // get fired trigger patches from the trigger selection task
  if(auto trgsel = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject("EmcalTriggerDecision"))){
    /*
    for(auto t : *(trgsel->GetListOfTriggerDecisions())){
      auto dec  = static_cast<PWG::EMCAL::AliEmcalTriggerDecision *>(t);
      std::cout << "Found trigger decision " << dec->GetName() << std::endl;
    }
    */
    for(auto t : l1triggers){
      auto decision = trgsel->FindTriggerDecision(t.Data());
      if(decision){
        patchhandlers[t] = decision->GetAcceptedPatches();
        if(energycomp < 0) {
          switch((decision->GetSelectionCuts()->GetSelectionMethod())){
            case PWG::EMCAL::AliEmcalTriggerSelectionCuts::kADC: energycomp = 0; break;
            case PWG::EMCAL::AliEmcalTriggerSelectionCuts::kEnergyOffline: energycomp = 1; break;
            case PWG::EMCAL::AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared: energycomp = 2; break;
          }
        }
      }
    }
  }

  auto supportedTriggers = GetSupportedTriggers(fUseExclusiveTriggers);
  Double_t energy, eta, phi, energyMaxEMCAL(0.), energyMaxDCAL(0.);
  const TList *selpatches(nullptr);
  AliVCluster *maxclusterEMCAL = nullptr,
              *maxclusterDCAL = nullptr;
  for(auto clust : GetClusterContainer(fNameClusterContainer.Data())->all()){
    //AliVCluster *clust = static_cast<AliVCluster *>(*clustIter);
    if(!clust->IsEMCAL()) continue;
    if(clust->GetIsExotic()) continue;
    if(!fClusterTimeRange.IsInRange(clust->GetTOF())) continue;

    // Distinguish energy definition
    switch(fEnergyDefinition){
    case kDefaultEnergy:
    	AliDebugStream(2) << GetName() << ": Using cluster energy definition: default" << std::endl;
    	energy = clust->E();
    	break;
    case kNonLinCorrEnergy:
    	AliDebugStream(2) << GetName() << ": Using cluster energy definition: corrected for non-linearity" << std::endl;
    	energy = clust->GetNonLinCorrEnergy();
    	break;
    case kHadCorrEnergy:
    	AliDebugStream(2) << GetName() << ": Using cluster energy definition: corrected for hadronic contribution" << std::endl;
    	energy = clust->GetHadCorrEnergy();
    	break;
    };

    AliDebugStream(2) << GetName() << ": Using energy " << energy << " (def: " << clust->E()
    		<< " | NL: " << clust->GetNonLinCorrEnergy()
			<< " | HD: " << clust->GetHadCorrEnergy()
			<< ")" << std::endl;

    TLorentzVector posvec;
    clust->GetMomentum(posvec, fVertex);
    posvec.SetE(energy);        // use energy definition as selected in the task
    eta = posvec.Eta();
    phi = posvec.Phi();
    if(phi < 0) phi += TMath::TwoPi();

    bool isEMCAL = phi < 3.8;
    if(isEMCAL) {
      if(!maxclusterEMCAL || (energy > energyMaxEMCAL)) {
        maxclusterEMCAL = clust;
        energyMaxEMCAL = energy;
      }
    } else {
      if(!maxclusterDCAL || (energy > energyMaxDCAL)){
        maxclusterDCAL = clust;
        energyMaxDCAL = energy;
      }
    }

    // fill histograms allEta
    for(const auto & trg : fSelectedTriggers){
      if(std::find(supportedTriggers.begin(), supportedTriggers.end(), trg) == supportedTriggers.end()) continue;
      selpatches = nullptr;
      for(auto t : l1triggers) {
        if(trg.Contains(t)) {
          auto patchdata = patchhandlers.find(t);
          if(patchdata != patchhandlers.end()){
            selpatches = patchdata->second;
          }
        }
      }
      for(auto trgclust : fTriggerClusters) {
        FillClusterHistograms(trg.Data(), energy, eta, phi, clust->GetTOF(), clust->GetNCells(), trgclust, selpatches, energycomp);
      }
    }
  }

  // Fill max cluster histogram 
  // in case not found fill also 0
  // Select a max. cluster in EMCAL and DCAL separately
  // and monitor both individually
  // In case of the combined triggers select the larger of the two
  // EMCAL
  double maxpointFull[6] = {-1., 0., fEventCentrality, 0., -1. -1.};
  std::vector<TString> combinedtriggers;
  if(maxclusterEMCAL) {
    maxpointFull[1] = 0;
    maxpointFull[3] = energyMaxEMCAL;
    TLorentzVector maxvector;
    maxclusterEMCAL->GetMomentum(maxvector, fVertex);
    maxpointFull[4] = maxvector.Eta();
    maxpointFull[5] = maxvector.Phi();
    if(maxpointFull[5] < 0) maxpointFull[5] += TMath::TwoPi();
    Int_t supermoduleID = -1;
    fGeom->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
    maxpointFull[0] = supermoduleID;
  }
  std::vector<double> maxpoint = {maxpointFull[0], maxpointFull[1], maxpointFull[2], maxpointFull[3]};
  if(fMonitorEtaPhi) {
    maxpoint.push_back(maxpointFull[4]);
    maxpoint.push_back(maxpointFull[5]);
  }
  // prepare trigger cluster
  maxpoint.push_back(0);
  int indexTrgCluster = maxpoint.size() - 1;
  for(const auto & trg : fSelectedTriggers){
    if(trg.Contains("ED")){
      combinedtriggers.push_back(trg);
      continue;
    } 
    if(std::find(supportedTriggers.begin(), supportedTriggers.end(), trg) == supportedTriggers.end()) continue;
    auto weight = GetTriggerWeight(trg.Data());
    for(auto trgclust : fTriggerClusters) {
      maxpoint[indexTrgCluster] = trgclust;
      fHistos->FillTHnSparse("hClusterTHnSparseMax" + trg, maxpoint.data(), weight);
    }
  }
  // DCAL
  if(maxclusterDCAL) {
    maxpointFull[1] = 1;
    maxpointFull[3] = energyMaxDCAL;
    TLorentzVector maxvector;
    maxclusterDCAL->GetMomentum(maxvector, fVertex);
    maxpointFull[4] = maxvector.Eta();
    maxpointFull[5] = maxvector.Phi();
    if(maxpointFull[5] < 0) maxpointFull[5] += TMath::TwoPi();
    Int_t supermoduleID = -1;
    fGeom->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
    maxpointFull[0] = supermoduleID;
  } else {
    // Reset max point
    maxpointFull[0] = -1.;
    maxpointFull[1] = 1.;
    maxpointFull[3] = 0.;
    maxpointFull[4] = -1.;
    maxpointFull[5] = -1.;
  }
  maxpoint[0] = maxpointFull[0]; 
  maxpoint[1] = maxpointFull[1];
  maxpoint[2] = maxpointFull[2];
  maxpoint[3] = maxpointFull[3];
  if(fMonitorEtaPhi) {
    maxpoint[4] = maxpointFull[4];
    maxpoint[5] = maxpointFull[5];
  }
  for(const auto & trg : fSelectedTriggers){
    if(trg.Contains("ED")) continue;
    if(std::find(supportedTriggers.begin(), supportedTriggers.end(), trg) == supportedTriggers.end()) continue;
    auto weight = GetTriggerWeight(trg.Data());
    for(auto trgclust : fTriggerClusters) {
      maxpoint[indexTrgCluster] = trgclust;
      fHistos->FillTHnSparse("hClusterTHnSparseMax" + trg, maxpoint.data(), weight);
    }
  }
  // handle combined trigger as the larger of the max. EMCAL or DCAL cluster
  if(combinedtriggers.size()) {
    AliVCluster *maxcluster = nullptr;
    Double_t energyMax = 0.;
    maxpointFull[1] = 2; // No selected cluster in event
    if(maxclusterEMCAL && maxclusterDCAL) {
      if(energyMaxEMCAL > energyMaxDCAL) {
        maxcluster = maxclusterEMCAL;
        energyMax = energyMaxEMCAL;
        maxpointFull[1] = 0;
      } else {
        maxcluster = maxclusterDCAL;
        energyMax = energyMaxDCAL;
        maxpointFull[1] = 1;
      }
    } else if(maxclusterEMCAL){
      maxcluster = maxclusterEMCAL;
      energyMax = energyMaxEMCAL;
      maxpointFull[1] = 0;
    } 
    else if(maxclusterDCAL) {
      maxcluster = maxclusterDCAL;
      energyMax = energyMaxDCAL;
      maxpointFull[1] = 1;
    }
    if(maxcluster) {
      maxpointFull[1] = 1;
      maxpointFull[3] = energyMax;
      TLorentzVector maxvector;
      maxcluster->GetMomentum(maxvector, fVertex);
      maxpointFull[4] = maxvector.Eta();
      maxpointFull[5] = maxvector.Phi();
      if(maxpointFull[5] < 0) maxpointFull[5] += TMath::TwoPi();
      Int_t supermoduleID = -1;
      fGeom->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
      maxpointFull[0] = supermoduleID;
    } else {
      // Reset max point
      maxpointFull[0] = -1.;
      maxpointFull[3] = 0.;
      maxpointFull[4] = -1.;
      maxpointFull[5] = -1.;
    }
    maxpoint[0] = maxpointFull[0]; 
    maxpoint[1] = maxpointFull[1];
    maxpoint[2] = maxpointFull[2];
    maxpoint[3] = maxpointFull[3];
    if(fMonitorEtaPhi) {
      maxpoint[4] = maxpointFull[4];
      maxpoint[5] = maxpointFull[5];
    }
    for(const auto & trg : combinedtriggers){
      if(std::find(supportedTriggers.begin(), supportedTriggers.end(), trg) == supportedTriggers.end()) continue;
      auto weight = GetTriggerWeight(trg.Data());
      for(auto trgclust : fTriggerClusters) {
        maxpoint[indexTrgCluster] = trgclust;
        fHistos->FillTHnSparse("hClusterTHnSparseMax" + trg, maxpoint.data(), weight);
      }
    }
  }
  return true;
}

void AliAnalysisTaskEmcalClustersRef::FillClusterHistograms(const TString &triggerclass, double energy, double eta, double phi, double clustertime, int ncell, int trgcluster, const TList *triggerPatches, int energycomp){
  std::vector<AliEMCALTriggerPatchInfo *> matchedPatches;
  if(fUseFiredTriggers && triggerPatches) {
    matchedPatches = CorrelateToTrigger(eta, phi, *triggerPatches);
  }
  auto hasTriggerPatch = matchedPatches.size() > 0;
  Int_t supermoduleID = -1;
  Double_t weight = GetTriggerWeight(triggerclass.Data());
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << triggerclass << std::endl;

  fGeom->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
  std::vector<double> point;
  point.push_back(static_cast<double>(supermoduleID));
  point.push_back(fEventCentrality);
  point.push_back(energy);
  if(fMonitorEtaPhi){
    point.push_back(eta);
    point.push_back(phi);
  }
  point.push_back(static_cast<double>(trgcluster));
  fHistos->FillTHnSparse("hClusterTHnSparseAll" + triggerclass, point.data(), weight);

  fHistos->FillTH2("hTimeEnergy" + triggerclass, clustertime, energy, weight);
  fHistos->FillTH2("hNCellEnergy" + triggerclass, ncell, energy, weight);
  if(fUseFiredTriggers && hasTriggerPatch){
    // find maximum trigger patch
    AliEMCALTriggerPatchInfo *maxpatch(nullptr);
    double maxenergy = 0;
    for(auto patch : matchedPatches) {
      double patche = 0;
      switch(energycomp){
        case 0: patche = patch->GetADCAmp(); break;
        case 1: patche = patch->GetPatchE(); break;
        case 2: patche = patch->GetSmearedEnergy(); break;
      };
      if(patche > maxenergy) {
        maxpatch = patch;
        maxenergy = patche;
      }
    }
    fHistos->FillTH2("hCorrClusterEPatchADC" + triggerclass, energy, maxpatch->GetADCAmp());
    fHistos->FillTH2("hCorrClusterEPatchE" + triggerclass, energy, maxpatch->GetPatchE());
    fHistos->FillTHnSparse("hClusterTHnSparseFired" + triggerclass, point.data(), weight);
  }
}

void AliAnalysisTaskEmcalClustersRef::UserFillHistosAfterEventSelection(){
  double v0amult = fInputEvent->GetVZEROData()->GetMTotV0A(),
         trackletmult = static_cast<double>(CountTracklets(-0.8, 0.8, 0., TMath::TwoPi())),
         emctrackletmult = static_cast<double>(CountTracklets(-0.8, 0.8, 1.4, TMath::Pi())),
         itsclustermult = fInputEvent->GetMultiplicity()->GetNumberOfSPDClusters(),
         emcclustermult = static_cast<double>(CountEmcalClusters(0.5)),
         emccellocc = static_cast<double>(this->GetEMCALCellOccupancy(0.1));

  auto supportedTriggers = GetSupportedTriggers(fUseExclusiveTriggers);

  for(const auto &t : fSelectedTriggers){
    if(std::find(supportedTriggers.begin(), supportedTriggers.end(), t) == supportedTriggers.end()) continue;
    Double_t weight = GetTriggerWeight(t.Data());
    fHistos->FillTH1("hEventCentrality" + t, fEventCentrality, weight);
    fHistos->FillTH1("hVertexZ" + t, fVertex[2], weight);
    
    for(auto trgclust : fTriggerClusters){
      fHistos->FillTH1("hTrgClustCounter" + t, static_cast<double>(trgclust), weight);
    }

    // Multiplicity correlation (no correction for downscaling)
    if(fDoFillMultiplicityHistograms){
      double data[6] = {v0amult, trackletmult, emctrackletmult, itsclustermult, emcclustermult, emccellocc};
      fHistos->FillTHnSparse("hMultiplicityCorrelation" + t, data);
    }
  }
}

std::vector<AliEMCALTriggerPatchInfo *> AliAnalysisTaskEmcalClustersRef::CorrelateToTrigger(Double_t etaclust, Double_t phiclust, const TList &triggerPatches) const {
  std::vector<AliEMCALTriggerPatchInfo *> foundpatches;
  for(auto patchIter : triggerPatches){
    Double_t boundaries[4];
    auto testpatch = static_cast<AliEMCALTriggerPatchInfo *>(patchIter);
    GetPatchBoundaries(*testpatch, boundaries);
    Double_t etamin = TMath::Min(boundaries[0], boundaries[1]),
        etamax = TMath::Max(boundaries[0], boundaries[1]),
        phimin = TMath::Min(boundaries[2], boundaries[3]),
        phimax = TMath::Max(boundaries[2], boundaries[3]);
    if(etaclust > etamin && etaclust < etamax && phiclust > phimin && phiclust < phimax){
      foundpatches.push_back(testpatch);
      break;
    }
  }
  return foundpatches;
}

void AliAnalysisTaskEmcalClustersRef::GetPatchBoundaries(AliEMCALTriggerPatchInfo &patch, Double_t *boundaries) const {
  boundaries[0] = patch.GetEtaMin();
  boundaries[1] = patch.GetEtaMax();
  boundaries[2] = patch.GetPhiMin();
  boundaries[3] = patch.GetPhiMax();
}

int AliAnalysisTaskEmcalClustersRef::CountEmcalClusters(double ecut){
	int nclusters = 0;
	for(auto clust : GetClusterContainer(fNameClusterContainer.Data())->all()){
	  if(!clust->IsEMCAL()) continue;
	  if(clust->GetIsExotic()) continue;
	  if(clust->E() > ecut) nclusters++;
	  nclusters++;
	}
	return nclusters;
}

int AliAnalysisTaskEmcalClustersRef::CountTracklets(double etamin, double etamax, double phimin, double phimax){
  int ntracklets = 0;
  AliVMultiplicity *mult = fInputEvent->GetMultiplicity();
  for(int itl = 0; itl < mult->GetNumberOfTracklets(); itl++){
    double eta = mult->GetEta(itl), phi = mult->GetPhi(itl);
    if(!(eta > etamin && eta < etamax)) continue;
    if(!(phi > phimin && phi < phimax)) continue;
    ntracklets++;
  }
  return ntracklets;
}

int AliAnalysisTaskEmcalClustersRef::GetEMCALCellOccupancy(double ecut){
  std::set<int> cellIDs;
  AliVCaloCells *emccells = fInputEvent->GetEMCALCells();
  for(short icell = 0; icell < emccells->GetNumberOfCells(); icell++){
    if(emccells->GetAmplitude(icell) > ecut){
      int cellID = emccells->GetCellNumber(icell);
      if(cellIDs.find(cellID) == cellIDs.end()) cellIDs.insert(cellID);
    }
  }
  return cellIDs.size();
}

void AliAnalysisTaskEmcalClustersRef::AddRequiredTriggerOverlap(const char *trigger){
  if(fRequiredOverlaps.FindObject(trigger)) return;
  fRequiredOverlaps.Add(new TObjString(trigger));
}

void AliAnalysisTaskEmcalClustersRef::AddExcludedTriggerOverlap(const char *trigger){
  if(fExcludedOverlaps.FindObject(trigger)) return;
  fExcludedOverlaps.Add(new TObjString(trigger));
}

AliAnalysisTaskEmcalClustersRef *AliAnalysisTaskEmcalClustersRef::AddTaskEmcalClustersRef(const TString &nclusters, const TString &suffix){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString clusName(nclusters == "usedefault" ? AliEmcalAnalysisFactory::ClusterContainerNameFactory(mgr->GetInputEventHandler()->InheritsFrom("AliAODInputHandler")) : nclusters);

  TString taskname = "emcalClusterQA_" + suffix;

  auto task = new AliAnalysisTaskEmcalClustersRef(taskname.Data());
  task->AddClusterContainer(clusName.Data());
  task->SetClusterContainer(clusName.Data());
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ClusterQA_" + TString(suffix);
  TString containername = "ClusterResults_" + TString(suffix);
  printf("Outfile: %s, container: %s\n", outfile.Data(), containername.Data());

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(containername.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}

AliAnalysisTaskEmcalClustersRef *AliAnalysisTaskEmcalClustersRef::AddTaskEmcalClustersRefDefault(const TString &nClusters){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  auto task = new AliAnalysisTaskEmcalClustersRef("emcalClusterQA");
  mgr->AddTask(task);

  // Adding cluster container
  TString clusName(nClusters == "usedefault" ? AliEmcalAnalysisFactory::ClusterContainerNameFactory(mgr->GetInputEventHandler()->InheritsFrom("AliAODInputHandler")) : nClusters);
  task->AddClusterContainer(clusName.Data());
  task->SetClusterContainer(clusName.Data());

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ClusterQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("ClusterResults", AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}


/**
 * Create new energy binning
 */
AliAnalysisTaskEmcalClustersRef::EnergyBinning::EnergyBinning():
  TCustomBinning()
{
  this->SetMinimum(0.);
  this->AddStep(1, 0.05);
  this->AddStep(2, 0.1);
  this->AddStep(4, 0.2);
  this->AddStep(7, 0.5);
  this->AddStep(16, 1);
  this->AddStep(32, 2);
  this->AddStep(40, 4);
  this->AddStep(50, 5);
  this->AddStep(100, 10);
  this->AddStep(200, 20);
}
