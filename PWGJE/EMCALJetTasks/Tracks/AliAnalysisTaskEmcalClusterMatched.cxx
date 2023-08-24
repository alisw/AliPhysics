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
#include <THashList.h>
#include <THistManager.h>
#include <TLinearBinning.h>

#include "AliAnalysisTaskEmcalClusterMatched.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliEMCALGeometry.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliTrackContainer.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliVTrack.h"

#include <iostream>

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalClusterMatched)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalClusterMatched::AliAnalysisTaskEmcalClusterMatched():
    AliAnalysisTaskEmcalTriggerBase(),
    fTrackSelectionGlobal(nullptr),
    fTrackSelectionTPConly(nullptr),
    fTimeCut(-50e6, 50e6),
    fEnableSumw2(false)
{

}

AliAnalysisTaskEmcalClusterMatched::AliAnalysisTaskEmcalClusterMatched(const char *name):
    AliAnalysisTaskEmcalTriggerBase(name),
    fTrackSelectionGlobal(nullptr),
    fTrackSelectionTPConly(nullptr),
    fTimeCut(-50e6, 50e6),
    fEnableSumw2(false)
{

}

AliAnalysisTaskEmcalClusterMatched::~AliAnalysisTaskEmcalClusterMatched() {
  if(fTrackSelectionGlobal) delete fTrackSelectionGlobal;
  if(fTrackSelectionTPConly) delete fTrackSelectionTPConly;
}

void AliAnalysisTaskEmcalClusterMatched::CreateUserObjects(){
  Bool_t isAOD = fInputHandler->IsA() == AliAODInputHandler::Class();
  if(!(fTrackSelectionGlobal || fTrackSelectionTPConly)) InitializeTrackSelections(isAOD);
  AddClusterContainer(fNameClusterContainer);
}

void AliAnalysisTaskEmcalClusterMatched::CreateUserHistos(){
  EnergyBinning enbinning;
  TLinearBinning smbinning(20, -0.5, 19.5), timebinning(1000, -0.5e-6, 0.5e-6);

  TString optionstring = fEnableSumw2 ? "s" : "";

  for(const auto &t :GetSupportedTriggers()){
    fHistos->CreateTH1("hEventCount" + t, "Event count for trigger " + t, 1, 0.5, 1.5, optionstring);
    fHistos->CreateTH1("hVertexZ" + t, "Position of the z-vertex for trigger" + t, 500, -50., 50., optionstring);
    fHistos->CreateTH2("hClusterEnergyTime" + t, "Cluster energy vs time for trigger " + t, enbinning, timebinning, optionstring);
    fHistos->CreateTH2("hClusterEnergyAllSM" + t, "Cluster energy vs supermodule for all clusters for trigger " + t, smbinning, enbinning, optionstring);
    fHistos->CreateTH2("hClusterEnergyMatchGlobalSM" + t, "Cluster energy vs supermodule for matched clusters (gl) for trigger " + t, smbinning, enbinning, optionstring);
    fHistos->CreateTH2("hClusterEnergyTimeMatchGlobal" + t, "Cluster energy vs time for matched clusters (gl) for trigger" + t, enbinning, timebinning, optionstring);
    fHistos->CreateTH2("hClusterEnergyMatchTPConlySM" + t, "Cluster energy vs supermodule for matched clusters (tpc) for trigger " + t, smbinning, enbinning, optionstring);
    fHistos->CreateTH2("hClusterEnergyTimeMatchTPConly" + t, "Cluster energy vs time for matched clusters (tpc) for trigger" + t, enbinning, timebinning, optionstring);
    fHistos->CreateTH2("hClusterEnergyMatchTPCexclSM" + t, "Cluster energy vs supermodule for matched clusters (tpc excl) for trigger " + t, smbinning, enbinning, optionstring);
    fHistos->CreateTH2("hClusterEnergyTimeMatchTPCexcl" + t, "Cluster energy vs time for matched clusters (tpc excl) for trigger " + t, enbinning, timebinning, optionstring);
    fHistos->CreateTH2("hClusterEnergyUnmatchSM" + t, "Cluster energy vs supermodule for unmatched clusters for trigger " + t, smbinning, enbinning, optionstring);
    fHistos->CreateTH2("hClusterEnergyTimeUnmatch" + t, "Cluster energy vs time for unmatched clusters for trigger " + t, enbinning, timebinning, optionstring);
  }
}

void AliAnalysisTaskEmcalClusterMatched::UserFillHistosAfterEventSelection() {
  for(const auto &t : fSelectedTriggers) {
    double weight = GetTriggerWeight(t.Data());
    fHistos->FillTH1("hEventCount" + t, 1, weight);
    fHistos->FillTH1("hVertexZ" + t, fVertex[2], weight);
  }
}

bool AliAnalysisTaskEmcalClusterMatched::Run(){
  for(auto clust : GetClusterContainer(fNameClusterContainer)->all()){
    if(!clust->IsEMCAL()) continue;
    if(clust->GetIsExotic()) continue;
    if(!fTimeCut.IsInRange(clust->GetTOF())) continue;

    double energy =  clust->GetNonLinCorrEnergy();
    TLorentzVector clustervec;
    clust->GetMomentum(clustervec,fVertex);
    int supermoduleID;
    fGeom->SuperModuleNumberFromEtaPhi(clustervec.Eta(), clustervec.Phi(), supermoduleID);
    for(const auto &t : fSelectedTriggers){
      double weight = GetTriggerWeight(t.Data());
      fHistos->FillTH2("hClusterEnergyAllSM" + t, supermoduleID, energy, weight);
      fHistos->FillTH2("hClusterEnergyTime" + t, energy, clust->GetTOF(), weight);
    }

    // check whether cluster was matched to global or TPC-only track
    // attention: tpc tracks always subset of global tracks, exclusive TPC tracks
    // handled separately
    int nglobal(0), ntpc(0);
    if(clust->GetNTracksMatched()){
      AliDebugStream(1) << "Cluster matched to " << clust->GetNTracksMatched() << " Tracks" << std::endl;
    }
    for(int itrk = 0; itrk < clust->GetNTracksMatched(); itrk++){
      AliVTrack *matched = static_cast<AliVTrack *>(clust->GetTrackMatched(itrk));
      if(!matched) {
        AliDebugStream(2) << "Track is null" << std::endl;
        continue;
      }
      if(fTrackSelectionGlobal->IsTrackAccepted(matched)) nglobal++;
      if(fTrackSelectionTPConly->IsTrackAccepted(matched)) ntpc++;
    }

    if(nglobal || ntpc){
      for(const auto &t : fSelectedTriggers){
        double weight = GetTriggerWeight(t.Data());
        if(nglobal){
          fHistos->FillTH2("hClusterEnergyMatchGlobalSM" + t, supermoduleID, energy, weight);
          fHistos->FillTH2("hClusterEnergyTimeMatchGlobal" + t, energy, clust->GetTOF(), weight);
        }
        if(ntpc){
          fHistos->FillTH2("hClusterEnergyMatchTPConlySM" + t, supermoduleID, energy, weight);
          fHistos->FillTH2("hClusterEnergyTimeMatchTPConly" + t, energy, clust->GetTOF(), weight);
        }
        if(ntpc && !nglobal){
          fHistos->FillTH2("hClusterEnergyMatchTPCexclSM" + t, supermoduleID, energy, weight);
          fHistos->FillTH2("hClusterEnergyTimeMatchTPCexcl" + t, energy, clust->GetTOF(), weight);
        }
      }
    } else {
      // Unmatched cluster
      for(const auto &t : fSelectedTriggers){
        double weight = GetTriggerWeight(t.Data());
        fHistos->FillTH2("hClusterEnergyUnmatchSM" + t, supermoduleID, energy, weight);
        fHistos->FillTH2("hClusterEnergyTimeUnmatch" + t, energy, clust->GetTOF(), weight);
      }
    }
  }
  return true;
}

void AliAnalysisTaskEmcalClusterMatched::InitializeTrackSelections(bool isAOD){
  if(isAOD){
    AliEmcalTrackSelectionAOD *globalAOD = new AliEmcalTrackSelectionAOD;
    globalAOD->AddFilterBit(AliAODTrack::kTrkGlobal);
    fTrackSelectionGlobal = globalAOD;

    AliEmcalTrackSelectionAOD *tpcOnlyAOD = new AliEmcalTrackSelectionAOD;
    tpcOnlyAOD->AddFilterBit(AliAODTrack::kTrkTPCOnly);
    fTrackSelectionTPConly = tpcOnlyAOD;
  } else {
    AliEmcalTrackSelectionESD *globalESD = new AliEmcalTrackSelectionESD;
    globalESD->AddTrackCuts(AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1));
    fTrackSelectionGlobal = globalESD;

    AliEmcalTrackSelectionESD *tpcOnlyESD = new AliEmcalTrackSelectionESD;
    tpcOnlyESD->AddTrackCuts(AliESDtrackCuts::GetStandardTPCOnlyTrackCuts());
    fTrackSelectionTPConly = tpcOnlyESD;
  }
}

AliAnalysisTaskEmcalClusterMatched *AliAnalysisTaskEmcalClusterMatched::AddTaskEmcalClusterMatched(const char *suffix){
  TString suffixstring = suffix;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisTaskEmcalClusterMatched *clustertask = new AliAnalysisTaskEmcalClusterMatched("EmcalMatchedClusterAnalysis" + suffixstring);
  mgr->AddTask(clustertask);

  TString outputfile = mgr->GetCommonFileName();
  outputfile += (":MatchedClusterResults" + suffixstring);
  mgr->ConnectInput(clustertask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(clustertask, 1, mgr->CreateContainer("MatchedClusterResultHistos" + suffixstring, AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile));
  return clustertask;
}

/**
 * Create new energy binning
 */
AliAnalysisTaskEmcalClusterMatched::EnergyBinning::EnergyBinning():
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
