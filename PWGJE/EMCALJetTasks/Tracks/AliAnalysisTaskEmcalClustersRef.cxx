/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <array>
#include <bitset>
#include <iostream>
#include <map>
#include <vector>

#include <TArrayD.h>
#include <TClonesArray.h>
#include <TGrid.h>
#include <THashList.h>
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TParameter.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliClusterContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalList.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"
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
    fNameClusterContainer(""),
    fCentralityRange(-999., 999.),
    fRequestCentrality(false),
    fEventCentrality(-1),
    fCentralityEstimator("V0M"),
    fEnergyDefinition(kDefaultEnergy)
{
}

AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef(const char *name) :
    AliAnalysisTaskEmcalTriggerBase(name),
    fNameClusterContainer(""),
    fCentralityRange(-999., 999.),
    fRequestCentrality(false),
    fEventCentrality(-1),
    fCentralityEstimator("V0M"),
    fEnergyDefinition(kDefaultEnergy)
{
}

AliAnalysisTaskEmcalClustersRef::~AliAnalysisTaskEmcalClustersRef() {
}


void AliAnalysisTaskEmcalClustersRef::CreateUserHistos(){

  EnergyBinning energybinning;
  TLinearBinning smbinning(21, -0.5, 20.5), etabinning(100, -0.7, 0.7);

  /*
   * Exclusive classes are defined as triggered events
   * in a class without lower threshold classes firing.
   * This is needed to make the triggers statistically
   * independent.
   */
  std::array<Double_t, 5> encuts = {1., 2., 5., 10., 20.};
  Int_t sectorsWithEMCAL[10] = {4, 5, 6, 7, 8, 9, 13, 14, 15, 16};
  for(auto trg : GetSupportedTriggers()){
    fHistos->CreateTH1(Form("hEventCount%s", trg.Data()), Form("Event count for trigger class %s", trg.Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hEventCentrality%s", trg.Data()), Form("Event centrality for trigger class %s", trg.Data()), 103, -2., 101.);
    fHistos->CreateTH1(Form("hVertexZ%s", trg.Data()), Form("z-position of the primary vertex for trigger class %s", trg.Data()), 200, -40., 40.);
    fHistos->CreateTH1(Form("hClusterEnergy%s", trg.Data()), Form("Cluster energy for trigger class %s", trg.Data()), energybinning);
    fHistos->CreateTH1(Form("hClusterET%s", trg.Data()), Form("Cluster transverse energy for trigger class %s", trg.Data()), energybinning);
    fHistos->CreateTH1(Form("hClusterEnergyFired%s", trg.Data()), Form("Cluster energy for trigger class %s, firing the trigger", trg.Data()), energybinning);
    fHistos->CreateTH1(Form("hClusterETFired%s", trg.Data()), Form("Cluster transverse energy for trigger class %s, firing the trigger", trg.Data()), energybinning);
    fHistos->CreateTH2(Form("hClusterEnergySM%s", trg.Data()), Form("Cluster energy versus supermodule for trigger class %s", trg.Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hClusterETSM%s", trg.Data()), Form("Cluster transverse energy versus supermodule for trigger class %s", trg.Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hClusterEnergyFiredSM%s", trg.Data()), Form("Cluster energy versus supermodule for trigger class %s, firing the trigger", trg.Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hClusterETFiredSM%s", trg.Data()), Form("Cluster transverse energy versus supermodule for trigger class %s, firing the trigger", trg.Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hEtaEnergy%s", trg.Data()), Form("Cluster energy vs. eta for trigger class %s", trg.Data()), etabinning, energybinning);
    fHistos->CreateTH2(Form("hEtaET%s", trg.Data()), Form("Cluster transverse energy vs. eta for trigger class %s", trg.Data()), etabinning, energybinning);
    fHistos->CreateTH2(Form("hEtaEnergyFired%s", trg.Data()), Form("Cluster energy vs. eta for trigger class %s, firing the trigger", trg.Data()), etabinning, energybinning);
    fHistos->CreateTH2(Form("hEtaETFired%s", trg.Data()), Form("Cluster transverse energy vs. eta for trigger class %s, firing the trigger", trg.Data()), etabinning, energybinning);
    for(int ism = 0; ism < 20; ism++){
      fHistos->CreateTH2(Form("hEtaEnergySM%d%s", ism, trg.Data()), Form("Cluster energy vs. eta in Supermodule %d for trigger %s", ism, trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETSM%d%s", ism, trg.Data()), Form("Cluster transverse energy vs. eta in Supermodule %d for trigger %s", ism, trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaEnergyFiredSM%d%s", ism, trg.Data()), Form("Cluster energy vs. eta in Supermodule %d for trigger %s, firing the trigger", ism, trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETFiredSM%d%s", ism, trg.Data()), Form("Cluster transverse energy vs. eta in Supermodule %d for trigger %s, firing the trigger", ism, trg.Data()), etabinning, energybinning);
    }
    for(int isec = 0; isec < 10; isec++){
      fHistos->CreateTH2(Form("hEtaEnergySec%d%s", sectorsWithEMCAL[isec], trg.Data()), Form("Cluster energy vs.eta in tracking sector %d for trigger %s", sectorsWithEMCAL[isec], trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETSec%d%s", sectorsWithEMCAL[isec], trg.Data()), Form("Cluster transverse energy vs.eta in tracking sector %d for trigger %s", sectorsWithEMCAL[isec], trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaEnergyFiredSec%d%s", sectorsWithEMCAL[isec], trg.Data()), Form("Cluster energy vs.eta in tracking sector %d for trigger %s, firing the trigger", sectorsWithEMCAL[isec], trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETFiredSec%d%s", sectorsWithEMCAL[isec], trg.Data()), Form("Cluster transverse energy vs.eta in tracking sector %d for trigger %s, firing the trigger", sectorsWithEMCAL[isec], trg.Data()), etabinning, energybinning);
    }
    for(auto ien : encuts){
      fHistos->CreateTH2(Form("hEtaPhi%dG%s", static_cast<int>(ien), trg.Data()), Form("cluster #eta-#phi map for clusters with energy larger than %f GeV/c for trigger class %s", ien, trg.Data()), 100, -0.7, 0.7, 200, 0, 2*TMath::Pi());
      fHistos->CreateTH2(Form("hEtaPhiFired%dG%s", static_cast<int>(ien), trg.Data()), Form("cluster #eta-#phi map for clusters fired the trigger with energy larger than %f GeV/c for trigger class %s", ien, trg.Data()), 200, -0.7, 0.7, 200, 0, 2*TMath::Pi());
    }
  }
}

bool AliAnalysisTaskEmcalClustersRef::IsUserEventSelected(){
  fEventCentrality = -1;
  if(fRequestCentrality){
    AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
    if(!mult){
      AliErrorStream() << GetName() << ": Centrality selection enabled but no centrality estimator found" << std::endl;
      return false;
    }
    if(mult->IsEventSelected()) return false;
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
  return true;
}

bool AliAnalysisTaskEmcalClustersRef::Run(){
  AliDebugStream(1) << GetName() << ": UserExec start" << std::endl;

  TList ej1patches, dj1patches, ej2patches, dj2patches, eg1patches, dg1patches, eg2patches, dg2patches;
  FindPatchesForTrigger("EJ2", fTriggerPatchInfo, ej2patches);
  FindPatchesForTrigger("DJ2", fTriggerPatchInfo, dj2patches);
  FindPatchesForTrigger("EJ1", fTriggerPatchInfo, ej1patches);
  FindPatchesForTrigger("DJ1", fTriggerPatchInfo, dj1patches);
  FindPatchesForTrigger("EG2", fTriggerPatchInfo, eg2patches);
  FindPatchesForTrigger("DG2", fTriggerPatchInfo, dg2patches);
  FindPatchesForTrigger("EG1", fTriggerPatchInfo, eg1patches);
  FindPatchesForTrigger("DG1", fTriggerPatchInfo, dg1patches);

  Double_t energy, et, eta, phi;
  TList *selpatches(nullptr);
  for(auto clust : GetClusterContainer(fNameClusterContainer.Data())->all()){
    //AliVCluster *clust = static_cast<AliVCluster *>(*clustIter);
    if(!clust->IsEMCAL()) continue;
    if(clust->GetIsExotic()) continue;

    // Distinguish energy definition
    switch(fEnergyDefinition){
    case kDefaultEnergy: energy = clust->E(); break;
    case kNonLinCorrEnergy: energy = clust->GetNonLinCorrEnergy(); break;
    case kHadCorrEnergy: energy = clust->GetHadCorrEnergy(); break;
    };

    TLorentzVector posvec;
    energy = clust->GetNonLinCorrEnergy();
    clust->GetMomentum(posvec, fVertex);
    posvec.SetE(energy);        // use energy definition as selected in the task
    et = posvec.Et();
    eta = posvec.Eta();
    phi = posvec.Phi();

    // fill histograms allEta
    for(const auto & trg : fSelectedTriggers){
      selpatches = nullptr;
      if(trg.Contains("EJ2")) selpatches = &ej2patches;
      if(trg.Contains("DJ2")) selpatches = &dj2patches;
      if(trg.Contains("EJ1")) selpatches = &ej1patches;
      if(trg.Contains("DJ1")) selpatches = &dj1patches;
      if(trg.Contains("EG2")) selpatches = &eg2patches;
      if(trg.Contains("DG2")) selpatches = &dg2patches;
      if(trg.Contains("EG1")) selpatches = &eg1patches;
      if(trg.Contains("DG1")) selpatches = &dg1patches;
      FillClusterHistograms(trg.Data(), energy, et, eta, phi, nullptr);
    }
  }
  return true;
}

void AliAnalysisTaskEmcalClustersRef::FillClusterHistograms(const TString &triggerclass, double energy, double transverseenergy, double eta, double phi, TList *fTriggerPatches){
  Bool_t hasTriggerPatch = fTriggerPatches  ? CorrelateToTrigger(eta, phi, fTriggerPatches) : kFALSE;
  Int_t supermoduleID = -1, sector = -1;
  Double_t weight = GetTriggerWeight(triggerclass);
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << triggerclass << std::endl;

  fGeom->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
  fHistos->FillTH1(Form("hClusterEnergy%s", triggerclass.Data()), energy, weight);
  fHistos->FillTH1(Form("hClusterET%s", triggerclass.Data()), transverseenergy, weight);
  fHistos->FillTH2(Form("hEtaEnergy%s", triggerclass.Data()), eta, energy, weight);
  fHistos->FillTH2(Form("hEtaET%s", triggerclass.Data()), eta, transverseenergy, weight);
  if(supermoduleID >= 0){
    fHistos->FillTH2(Form("hClusterEnergySM%s", triggerclass.Data()), supermoduleID, energy, weight);
    fHistos->FillTH2(Form("hClusterETSM%s", triggerclass.Data()), supermoduleID, transverseenergy, weight);
    fHistos->FillTH2(Form("hEtaEnergySM%d%s", supermoduleID, triggerclass.Data()), eta, energy, weight);
    fHistos->FillTH2(Form("hEtaETSM%d%s", supermoduleID, triggerclass.Data()), eta, transverseenergy, weight);
    if(supermoduleID < 12)
      sector = 4 + int(supermoduleID/2); // EMCAL
    else
      sector = 13 + int((supermoduleID-12)/2);  // DCAL
    fHistos->FillTH2(Form("hEtaEnergySec%d%s", sector, triggerclass.Data()), eta, energy, weight);
    fHistos->FillTH2(Form("hEtaETSec%d%s", sector, triggerclass.Data()), eta, transverseenergy, weight);
  }
  if(hasTriggerPatch){
    fHistos->FillTH1(Form("hClusterEnergyFired%s", triggerclass.Data()), energy, weight);
    fHistos->FillTH1(Form("hClusterETFired%s", triggerclass.Data()), energy, weight);
    fHistos->FillTH2(Form("hEtaEnergyFired%s", triggerclass.Data()), eta, energy, weight);
    fHistos->FillTH2(Form("hEtaETFired%s", triggerclass.Data()), eta, energy, weight);
    if(supermoduleID >= 0){
      fHistos->FillTH2(Form("hClusterEnergyFiredSM%s", triggerclass.Data()), supermoduleID, energy, weight);
      fHistos->FillTH2(Form("hClusterETFiredSM%s", triggerclass.Data()), supermoduleID, transverseenergy, weight);
      fHistos->FillTH2(Form("hEtaEnergyFiredSM%d%s", supermoduleID, triggerclass.Data()), eta, energy,weight);
      fHistos->FillTH2(Form("hEtaETFiredSM%d%s", supermoduleID, triggerclass.Data()), eta, transverseenergy, weight);
      fHistos->FillTH2(Form("hEtaEnergyFiredSec%d%s", sector, triggerclass.Data()), eta, energy, weight);
      fHistos->FillTH2(Form("hEtaETFiredSec%d%s", sector, triggerclass.Data()), eta, transverseenergy, weight);
    }
  }
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(int ien = 0; ien < 5; ien++){
    if(energy > encuts[ien]){
      fHistos->FillTH2(Form("hEtaPhi%dG%s", static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi, weight);
      if(hasTriggerPatch){
        fHistos->FillTH2(Form("hEtaPhiFired%dG%s", static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi, weight);
      }
    }
  }
}

void AliAnalysisTaskEmcalClustersRef::UserFillHistosAfterEventSelection(){
  for(const auto &t : fSelectedTriggers){
    Double_t weight = GetTriggerWeight(t);
    fHistos->FillTH1(Form("hEventCount%s", t.Data()), 1, weight);
    fHistos->FillTH1(Form("hEventCentrality%s", t.Data()), fEventCentrality, weight);
    fHistos->FillTH1(Form("hVertexZ%s", t.Data()), fVertex[2], weight);
  }
}

Bool_t AliAnalysisTaskEmcalClustersRef::CorrelateToTrigger(Double_t etaclust, Double_t phiclust, TList *fTriggerPatches) const {
  Bool_t hasfound = kFALSE;
  for(TIter patchIter = TIter(fTriggerPatches).Begin(); patchIter != TIter::End(); ++patchIter){
    Double_t boundaries[4];
    GetPatchBoundaries(*patchIter, boundaries);
    Double_t etamin = TMath::Min(boundaries[0], boundaries[1]),
        etamax = TMath::Max(boundaries[0], boundaries[1]),
        phimin = TMath::Min(boundaries[2], boundaries[3]),
        phimax = TMath::Max(boundaries[2], boundaries[3]);
    if(etaclust > etamin && etaclust < etamax && phiclust > phimin && phiclust < phimax){
      hasfound = kTRUE;
      break;
    }
  }
  return hasfound;
}

void AliAnalysisTaskEmcalClustersRef::FindPatchesForTrigger(TString triggerclass, const TClonesArray * fTriggerPatches, TList &foundtriggers) const {
  foundtriggers.Clear();
  if(!fTriggerPatches) return;
  AliEmcalTriggerOfflineSelection::EmcalTriggerClass myclass = AliEmcalTriggerOfflineSelection::kTrgEL0;
  if(triggerclass == "EG1") myclass = AliEmcalTriggerOfflineSelection::kTrgEG1;
  if(triggerclass == "EG2") myclass = AliEmcalTriggerOfflineSelection::kTrgEG2;
  if(triggerclass == "EJ1") myclass = AliEmcalTriggerOfflineSelection::kTrgEJ1;
  if(triggerclass == "EJ2") myclass = AliEmcalTriggerOfflineSelection::kTrgEJ2;
  if(triggerclass == "DMC7") myclass = AliEmcalTriggerOfflineSelection::kTrgDL0;
  if(triggerclass == "DG1") myclass = AliEmcalTriggerOfflineSelection::kTrgDG1;
  if(triggerclass == "DG2") myclass = AliEmcalTriggerOfflineSelection::kTrgDG2;
  if(triggerclass == "DJ1") myclass = AliEmcalTriggerOfflineSelection::kTrgDJ1;
  if(triggerclass == "DJ2") myclass = AliEmcalTriggerOfflineSelection::kTrgDJ2;
  for(TIter patchiter = TIter(fTriggerPatches).Begin(); patchiter != TIter::End(); ++patchiter){
    if(!IsOfflineSimplePatch(*patchiter)) continue;
    if(AliEmcalTriggerOfflineSelection::IsDCAL(myclass)){
      if(!SelectDCALPatch(*patchiter)) continue;
    } else {
      if(SelectDCALPatch(*patchiter)) continue;
    }
    if(AliEmcalTriggerOfflineSelection::IsSingleShower(myclass)){
      if(!SelectSingleShowerPatch(*patchiter)) continue;
    } else {
      if(!SelectJetPatch(*patchiter)) continue;
    }
    double threshold = fTriggerSelection ? fTriggerSelection->GetThresholdForTrigger(myclass) : -1;
    if(GetPatchEnergy(*patchiter) > threshold) foundtriggers.Add(*patchiter);
  }
}

void AliAnalysisTaskEmcalClustersRef::GetPatchBoundaries(TObject *o, Double_t *boundaries) const {
  AliEMCALTriggerPatchInfo *patch= dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  boundaries[0] = patch->GetEtaMin();
  boundaries[1] = patch->GetEtaMax();
  boundaries[2] = patch->GetPhiMin();
  boundaries[3] = patch->GetPhiMax();
}

bool AliAnalysisTaskEmcalClustersRef::IsOfflineSimplePatch(TObject *o) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->IsOfflineSimple();
}

bool AliAnalysisTaskEmcalClustersRef::SelectDCALPatch(TObject *o) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->GetRowStart() >= 64;
}

bool AliAnalysisTaskEmcalClustersRef::SelectSingleShowerPatch(TObject *o) const{
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->IsGammaLowSimple();
}

bool AliAnalysisTaskEmcalClustersRef::SelectJetPatch(TObject *o) const{
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  if(!patch->IsOfflineSimple()) return false;
  return patch->IsJetLowSimple();
}

double AliAnalysisTaskEmcalClustersRef::GetPatchEnergy(TObject *o) const {
  double energy = 0.;
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  energy = patch->GetPatchE();
  return energy;
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
