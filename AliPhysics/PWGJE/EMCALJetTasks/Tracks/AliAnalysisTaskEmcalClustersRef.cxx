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
#include "AliEmcalTriggerOfflineSelection.h"
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

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef)
/// \endcond

namespace EMCalTriggerPtAnalysis {

AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef() :
    AliAnalysisTaskEmcalTriggerBase(),
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
  TLinearBinning smbinning(21, -0.5, 20.5), etabinning(100, -0.7, 0.7), phibinning(200, 0., TMath::TwoPi()), timebinning(1000, -500e-9, 500e-9), ncellbinning(101, -0.5, 100.5);
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
  const TBinning *clustallbinning[6] = {&smbinning, &centralitybinning, &energybinning, &etabinning, &phibinning, &trgclustbinning};
  AliDebugStream(1) << "Using exclusive triggers: " << (fUseExclusiveTriggers ? "yes" : "no") << std::endl;
  for(auto trg : GetSupportedTriggers(fUseExclusiveTriggers)){
    AliDebugStream(1) << "Creating histograms for trigger " << trg << std::endl;
    fHistos->CreateTH1("hTrgClustCounter" + trg, "Event counter in trigger cluster " + trg, trgclustbinning);
    fHistos->CreateTH1("hEventCentrality" + trg, "Event centrality for trigger class " + trg, 103, -2., 101., optionstring);
    fHistos->CreateTH1("hVertexZ" + trg, "z-position of the primary vertex for trigger class " + trg, 200, -40., 40., optionstring);
    if(this->fDoFillMultiplicityHistograms) fHistos->CreateTHnSparse("hMultiplicityCorrelation" + trg, "Multiplicity correlation for trigger" + trg, 6, multbinning);
    fHistos->CreateTHnSparse("hClusterTHnSparseAll" + trg, "Cluster THnSparse (all) for trigger" + trg, 6, clustallbinning);
    fHistos->CreateTHnSparse("hClusterTHnSparseMax" + trg, "Cluster THnSparse (max) for trigger" + trg, 6, clustallbinning);
    if(fUseFiredTriggers) fHistos->CreateTHnSparse("hClusterTHnSparseFired" + trg, "Cluster THnSparse (firing) for trigger" + trg, 6, clustallbinning);
    fHistos->CreateTH2("hTimeEnergy" + trg, "Cluster time vs. energy for trigger class " + trg, timebinning, energybinning, optionstring);
    fHistos->CreateTH2("hNCellEnergy" + trg, "Cluster number of cells vs energy for trigger class " + trg, ncellbinning, energybinning, optionstring);
    if(fUseFiredTriggers){
      fHistos->CreateTH2("hCorrClusterEPatchADC" + trg, "Correlation between cluster E and patch ADC for trigger " + trg, energybinning, adcbinning);
      fHistos->CreateTH2("hCorrClusterEPatchE" + trg, "Correlation between cluster E and patch E for trigger " + trg, energybinning, energybinning);
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
    Bool_t isCENT(false), isCENTNOTRD(false), isCALO(false), isCALOFAST(false);
    for(auto trg : PWG::EMCAL::Triggerinfo::DecodeTriggerString(fInputEvent->GetFiredTriggerClasses().Data())){
      auto trgclust = trg.Triggercluster();
      if(trgclust == "CENT" && !isCENT){
        isCENT = true;
        fTriggerClusters.emplace_back(kTrgClusterCENT);
      }
      if(trgclust == "CENTNOTRD" && ! isCENTNOTRD){
        isCENTNOTRD = true;
        fTriggerClusters.emplace_back(kTrgClusterCENTNOTRD);
      }
      if(trgclust == "CALO" && !isCALO){
        isCALO = true;
        fTriggerClusters.emplace_back(kTrgClusterCALO);
      }
      if(trgclust == "CALOFAST" && ! isCALOFAST) {
        isCALOFAST = true;
        fTriggerClusters.emplace_back();
      } 
    }
    // Mixed clusters
    if(isCENT || isCENTNOTRD) {
      if(isCENT) {
        if(isCENTNOTRD) fTriggerClusters.emplace_back(kTrgClusterCENTBOTH);
        else fTriggerClusters.emplace_back(kTrgClusterOnlyCENT);
      } else fTriggerClusters.emplace_back(kTrgClusterOnlyCENTNOTRD);
    }
    if(isCALO || isCALOFAST){
      if(isCALO) {
        if(isCALOFAST) fTriggerClusters.emplace_back(kTrgClusterCALOBOTH);
        else fTriggerClusters.emplace_back(kTrgClusterOnlyCALO);
      } else fTriggerClusters.emplace_back(kTrgClusterOnlyCALOFAST);
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
          switch(decision->GetSelectionCuts()->GetSelectionMethod()){
            case PWG::EMCAL::AliEmcalTriggerSelectionCuts::kADC: energycomp = 0; break;
            case PWG::EMCAL::AliEmcalTriggerSelectionCuts::kEnergyOffline: energycomp = 1; break;
            case PWG::EMCAL::AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared: energycomp = 2; break;
          }
        }
      }
    }
  }

  auto supportedTriggers = GetSupportedTriggers(fUseExclusiveTriggers);
  Double_t energy, et, eta, phi;
  const TList *selpatches(nullptr);
  AliVCluster *maxcluster = nullptr; 
  for(auto clust : GetClusterContainer(fNameClusterContainer.Data())->all()){
    //AliVCluster *clust = static_cast<AliVCluster *>(*clustIter);
    if(!clust->IsEMCAL()) continue;
    if(clust->GetIsExotic()) continue;
    if(!fClusterTimeRange.IsInRange(clust->GetTOF())) continue;
    if(!maxcluster || clust->E() > maxcluster->E()) maxcluster = clust;

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
  double maxpoint[6] = {-1, fEventCentrality, 0, -1. -1., 0};
  if(maxcluster) {
    maxpoint[2] = maxcluster->E();
    TLorentzVector maxvector;
    maxcluster->GetMomentum(maxvector, fVertex);
    maxpoint[3] = maxvector.Eta();
    maxpoint[4] = maxvector.Phi();
    if(maxpoint[4] < 0) maxpoint[4] += TMath::TwoPi();
    Int_t supermoduleID = -1;
    fGeom->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
    maxpoint[0] = supermoduleID;
  }
  for(const auto & trg : fSelectedTriggers){
    if(std::find(supportedTriggers.begin(), supportedTriggers.end(), trg) == supportedTriggers.end()) continue;
    auto weight = GetTriggerWeight(trg);
    for(auto trgclust : fTriggerClusters) {
      maxpoint[5] = trgclust;
      fHistos->FillTHnSparse("hClusterTHnSparseMax" + trg, maxpoint, weight);
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
  Double_t weight = GetTriggerWeight(triggerclass);
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << triggerclass << std::endl;

  fGeom->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
  double point[6] = {static_cast<double>(supermoduleID), fEventCentrality, energy, eta, phi, static_cast<double>(trgcluster)};
  fHistos->FillTHnSparse("hClusterTHnSparseAll" + triggerclass, point, weight);

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
    fHistos->FillTHnSparse("hClusterTHnSparseFired" + triggerclass, point, weight);
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
    Double_t weight = GetTriggerWeight(t);
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

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef(taskname.Data());
  task->AddClusterContainer(clusName);
  task->SetClusterContainer(clusName);
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

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef("emcalClusterQA");
  mgr->AddTask(task);

  // Adding cluster container
  TString clusName(nClusters == "usedefault" ? AliEmcalAnalysisFactory::ClusterContainerNameFactory(mgr->GetInputEventHandler()->InheritsFrom("AliAODInputHandler")) : nClusters);
  task->AddClusterContainer(clusName.Data());
  task->SetClusterContainer(clusName);

  // Set Energy thresholds for additional patch selection:
  // These are events with offline patches of a given type where the trigger reached already the plateau
  // These numers are determined as:
  // EMC7: 3.5 GeV
  // EG1:  14 GeV
  // EG2:  8 GeV
  // EJ1:  22 GeV
  // EJ2:  12 GeV
  task->SetOfflineTriggerSelection(
      EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TriggerSelectionFactory(5, 14, 8, 22, 12)
  );

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


} /* namespace EMCalTriggerPtAnalysis */
